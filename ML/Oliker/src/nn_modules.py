import math
import pickle
from functools import partial
from typing import Any, Callable, List, Optional, Sequence

import numpy as np
import torch
from torch import nn, Tensor
from torch.nn import functional as F

from torchvision.ops.misc import Conv2dNormActivation, Permute
from torchvision.ops.stochastic_depth import StochasticDepth

from torchvision.models._utils import _ovewrite_named_param
from tqdm import trange

from src.utils_scaling import *
from src.nn_blocks import CNBlockConfig, CNBlock, LayerNorm2d


class GenModule(nn.Module):
    def __init__(
                self,
                output_linear_coeffs: List[float]=[1, 0],
                output_dim: int = 16,
                **kwargs
                ):
        super().__init__()

        self.output_dim = output_dim
        self.output_linear_coeffs = output_linear_coeffs

    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, (nn.Conv2d, nn.Linear)):
                nn.init.trunc_normal_(m.weight, std=0.02)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)

    def batch_eval(self, inputs, batch_size):
        num_imgs = inputs.shape[0]
        num_batches = math.ceil(num_imgs / batch_size)
        progress = trange(num_batches)

        outputs = torch.empty_like(inputs)

        start_ind = 0
        for _ in progress:
            end_ind = start_ind + batch_size
            end_ind = end_ind if end_ind < num_imgs else num_imgs
            cur_inds = range(start_ind, end_ind)

            outputs[cur_inds] = self(inputs[cur_inds])
            start_ind = end_ind

        return outputs

    def batch_scaled_eval(self, inputs, batch_size):
        outputs = self.batch_eval(inputs, batch_size)
        return linear_inverse_scaling(outputs, self.output_linear_coeffs)

    def scaled_eval(self, x):
        res = linear_inverse_scaling(self(x), self.output_linear_coeffs)
        res[x <= 0] = -np.inf
        return res

    def _pack_params(self):
        net_params = {'output_dim': self.output_dim,
                      'output_linear_coeffs': self.output_linear_coeffs
                      }
        return net_params

    def save(self, file_name):
        net_params = self._pack_params()

        with open(file_name, 'wb') as file:
            pickle.dump(net_params, file)
            pickle.dump(self.state_dict(), file)

class ConvNeXt(GenModule):

    def __init__(
        self,
        block_setting: List[CNBlockConfig],
        stochastic_depth_prob: float = 0.5,
        layer_scale: float = 1e-6,
        block: Optional[Callable[..., nn.Module]] = None,
        norm_layer: Optional[Callable[..., nn.Module]] = None,
        **kwargs
    ) -> None:
        super().__init__(**kwargs)

        self.block_setting = block_setting

        if not block_setting:
            raise ValueError("The block_setting should not be empty")
        elif not (isinstance(block_setting, Sequence) and all([isinstance(s, CNBlockConfig) for s in block_setting])):
            raise TypeError("The block_setting should be List[CNBlockConfig]")

        if block is None:
            block = CNBlock

        if norm_layer is None:
            norm_layer = partial(LayerNorm2d, eps=1e-6)

        layers: List[nn.Module] = []

        # Stem
        firstconv_output_channels = block_setting[0].input_channels
        self.init_layer = Conv2dNormActivation(
                                1,
                                firstconv_output_channels,
                                kernel_size=5,
                                stride=1,
                                padding=None,
                                norm_layer=norm_layer,
                                activation_layer=None,
                                bias=True,
                            )

        total_stage_blocks = sum(cnf.num_layers for cnf in block_setting)
        stage_block_id = 0
        for cnf in block_setting:
            # Bottlenecks
            stage: List[nn.Module] = []
            for _ in range(cnf.num_layers):
                # adjust stochastic depth probability based on the depth of the stage block
                sd_prob = stochastic_depth_prob * stage_block_id / (total_stage_blocks - 1.0)
                stage.append(block(cnf.input_channels, layer_scale, sd_prob))
                stage_block_id += 1
            layers.append(nn.Sequential(*stage))
            if cnf.out_channels is not None:
                # Downsampling
                layers.append(
                    nn.Sequential(
                        norm_layer(cnf.input_channels),
                        nn.Conv2d(cnf.input_channels, cnf.out_channels, kernel_size=2, stride=2),
                    )
                )

        self.features = nn.Sequential(*layers)
        self.avgpool = nn.AdaptiveAvgPool2d(1)

        lastblock = block_setting[-1]
        lastconv_output_channels = (
            lastblock.out_channels if lastblock.out_channels is not None else lastblock.input_channels
        )
        self.generator = nn.Sequential(
            norm_layer(lastconv_output_channels), nn.Flatten(1), nn.Linear(lastconv_output_channels, self.output_dim * self.output_dim)
        )

        self._init_weights()

    def _forward_impl(self, x: Tensor) -> Tensor:
        x = self.init_layer(x)
        x = self.features(x)
        x = self.avgpool(x)
        x = self.generator(x)
        return x.view(-1, 1, self.output_dim, self.output_dim)

    def forward(self, x: Tensor) -> Tensor:
        return self._forward_impl(x)

    def _pack_params(self):
        net_params = super()._pack_params()
        net_params.update({'block_setting': self.block_setting})
        return net_params

    @staticmethod
    def from_file(file_name):
        with open(file_name, 'rb') as file:
            net_params = pickle.load(file)
            state_dict = pickle.load(file)

        torch_dtype = list(state_dict.items())[0][1].dtype

        m = ConvNeXt(block_setting=net_params['block_setting'], output_dim=net_params['output_dim'], output_linear_coeffs=net_params['output_linear_coeffs'])
        m.to(dtype=torch_dtype)
        m.load_state_dict(state_dict)
        return m


class OlikerUNet(GenModule):

    def __init__(self, net_depth=2, net_width=2, base_degree=1, **kwargs):
        super().__init__(**kwargs)

        self.net_depth = net_depth
        self.net_width = net_width
        self.base_degree = base_degree

        self._init_layers()

    def forward(self, input):
        x = input

        for i in range(self.net_width):
            saved_vals = []
            for j in range(self.net_depth):
                c = self.down_layers[i][j][0](x)
                x = self.down_layers[i][j][1](c)
                saved_vals.append(c)

            x = self.bot_layers[i](x)

            for j in range(self.net_depth):
                x = self.up_layers[i][j][0](x)
                x = torch.cat((x, saved_vals[-(j + 1)]), 1)
                x = self.up_layers[i][j][1](x)

            x = self.end_layers[i](x)

            # # Segments connection
            # if i > 0:
            #     x = torch.cat((x, input), 1)
            #     x = self.middle_layers[i-1](x)

        return x

    @staticmethod
    def from_file(file_name):
        with open(file_name, 'rb') as file:
            net_params = pickle.load(file)
            state_dict = pickle.load(file)

        torch_dtype = list(state_dict.items())[0][1].dtype

        m = OlikerUNet(output_linear_coeffs=net_params['output_linear_coeffs'],
                       output_dim=net_params['output_dim'],
                       net_depth=net_params['net_depth'],
                       net_width=net_params['net_width'],
                       base_degree=net_params['base_degree'])
        m.to(dtype=torch_dtype)
        m.load_state_dict(state_dict)

        return m

    def _init_layers(self):
        self.up_layers = nn.ModuleList([])
        self.down_layers = nn.ModuleList([])
        self.bot_layers = nn.ModuleList([])
        self.end_layers = nn.ModuleList([])
        self.middle_layers = nn.ModuleList([])

        end_layer_degree = self.net_depth + self.base_degree + 3
        for i in range(self.net_width):
            cur_up_layers = nn.ModuleList([])
            cur_down_layers = nn.ModuleList([])
            nbi = 2 ** (self.net_depth + self.base_degree)

            for j in range(self.net_depth):
                ndil = 1 if j == 0 else 2 ** (j + self.base_degree)
                ndol = 2 ** (j + self.base_degree + 1)
                nuil = 2 ** (self.net_depth + self.base_degree + 1) if j == 0 else 2 ** (
                        self.net_depth + self.base_degree + 2 - j)
                nuol = 2 ** (self.net_depth + self.base_degree - j)

                d1 = nn.Sequential(
                    nn.Conv2d(ndil, ndol, (3, 3), padding='same'),
                    nn.ReLU(),
                    nn.Dropout(p=0.1),
                    nn.Conv2d(ndol, ndol, (3, 3), padding='same'),
                    nn.ReLU(),
                )
                d2 = nn.Sequential(
                    nn.BatchNorm2d(ndol),
                    nn.ReLU(),
                    nn.MaxPool2d((2, 2))
                )

                u1 = nn.ConvTranspose2d(nuil, nuol, (2, 2), stride=(2, 2))
                u2 = nn.Sequential(
                    nn.BatchNorm2d(2 * nuol),
                    nn.ReLU()
                )

                cur_down_layers.append(nn.ModuleList([d1, d2]))
                cur_up_layers.append(nn.ModuleList([u1, u2]))

            self.up_layers.append(cur_up_layers)
            self.down_layers.append(cur_down_layers)

            self.bot_layers.append(nn.Sequential(
                nn.Conv2d(nbi, 2 * nbi, (3, 3), padding='same'),
                nn.ReLU(),
                nn.BatchNorm2d(2 * nbi),
                nn.ReLU(),
                nn.Dropout(p=0.3),
                nn.Conv2d(2 * nbi, 2 * nbi, (3, 3), padding='same')
            ))

            self.end_layers.append(
                nn.Sequential(nn.Conv2d(2 ** (self.base_degree + 2), 2 ** end_layer_degree, (1, 1)),
                              nn.ReLU(),
                              nn.Conv2d(2 ** end_layer_degree, 1, (1, 1)),
                              nn.ReLU()
                              )
            )

            self.middle_layers.append(nn.Sequential(
                nn.Conv2d(2, 1, (1, 1)),
                nn.ReLU()
            ))

    def _pack_params(self):
        net_params = super()._pack_params()
        net_params.update({'net_depth': self.net_depth,
                           'net_width': self.net_width,
                           'base_degree': self.base_degree})
        return net_params


class UNext(GenModule):

    def __init__(self,
                 block_setting: List[CNBlockConfig],
                 stochastic_depth_prob: float = 0,
                 layer_scale: float = 1e-6,
                 **kwargs):
        super().__init__(**kwargs)

        self.block_setting = block_setting
        self.stochastic_depth_prob = stochastic_depth_prob
        self.layer_scale = layer_scale

        self._init_layers()

    def forward(self, input):
        x = input

        x = self.init_layer(x)

        saved_vals = []
        for j in range(self.net_depth-1):
            x = self.down_layers[j](x)
            saved_vals.append(x)

        x = self.bot_layer(x)

        for j in range(self.net_depth-1):
            x = torch.cat((x, saved_vals[-(j + 1)]), 1)
            x = self.up_layers[j](x)

        x = self.end_layers[0](x)
        x = torch.cat((x, input), 1)
        x = self.end_layers[1](x)

        return x

    @staticmethod
    def from_file(file_name):
        with open(file_name, 'rb') as file:
            net_params = pickle.load(file)
            state_dict = pickle.load(file)

        torch_dtype = list(state_dict.items())[0][1].dtype

        m = UNext(output_linear_coeffs=net_params['output_linear_coeffs'],
                  output_dim=net_params['output_dim'],
                  block_setting=net_params['block_setting'],
                  stochastic_depth_prob=net_params['stochastic_depth_prob'],
                  layer_scale=net_params['layer_scale'])
        m.to(dtype=torch_dtype)
        m.load_state_dict(state_dict)

        return m

    def _init_layers(self):
        self.net_depth = math.floor(len(self.block_setting)/2) + 1

        self.up_layers = nn.ModuleList([])
        self.down_layers = nn.ModuleList([])
        self.bot_layer = []
        self.end_layers = nn.ModuleList([])
        self.init_layer = []

        block = CNBlock
        norm_layer = partial(LayerNorm2d, eps=1e-6)

        # Stem
        firstconv_output_channels = self.block_setting[0].input_channels
        self.init_layer = Conv2dNormActivation(
                                                1,
                                                firstconv_output_channels,
                                                kernel_size=5,
                                                stride=1,
                                                padding=None,
                                                norm_layer=norm_layer,
                                                activation_layer=None,
                                                bias=True,
                                                )

        total_stage_blocks = sum(cnf.num_layers for cnf in self.block_setting)
        stage_block_id = 0
        stage_id = 1
        for cnf in self.block_setting:
            # Bottlenecks
            stage: List[nn.Module] = []
            for _ in range(cnf.num_layers):
                # adjust stochastic depth probability based on the depth of the stage block
                sd_prob = self.stochastic_depth_prob * stage_block_id / (total_stage_blocks - 1.0)
                stage.append(block(cnf.input_channels, self.layer_scale, sd_prob))
                stage_block_id += 1

            if stage_id < self.net_depth:
                # Downsampling
                stage.append(
                    nn.Sequential(
                        norm_layer(cnf.input_channels),
                        nn.Conv2d(cnf.input_channels, cnf.out_channels, kernel_size=2, stride=2),
                    )
                )
                self.down_layers.append(nn.Sequential(*stage))
            elif stage_id == self.net_depth:
                # Bottom Layer Upsampling
                stage.append(
                    nn.Sequential(
                        norm_layer(cnf.input_channels),
                        # nn.ConvTranspose2d(cnf.input_channels, cnf.out_channels, (2, 2), stride=(2, 2))
                    )
                )
                self.bot_layer = nn.Sequential(*stage)
            else:
                # Upsampling
                stage.append(
                    nn.Sequential(
                        norm_layer(cnf.input_channels),
                        nn.ConvTranspose2d(cnf.input_channels, cnf.out_channels, (2, 2), stride=(2, 2))
                    )
                )
                self.up_layers.append(nn.Sequential(*stage))

            stage_id += 1

        lastconv_output_channels = self.block_setting[-1].out_channels
        self.end_layers.append(nn.Sequential(
                                    norm_layer(lastconv_output_channels),
                                    nn.Conv2d(lastconv_output_channels, 1, kernel_size=3, padding='same', bias=True)
                                ))
        self.end_layers.append(nn.Conv2d(2, 1, kernel_size=3, padding=1, bias=True))


        self._init_weights()


    def _pack_params(self):
        net_params = super()._pack_params()
        net_params.update({'block_setting': self.block_setting,
                           'stochastic_depth_prob': self.stochastic_depth_prob,
                           'layer_scale': self.layer_scale})
        return net_params


def convnext_base(output_dim=16, output_linear_coeffs=[1,0], weights=None, progress: bool = True, **kwargs: Any) -> ConvNeXt:
    block_setting = [
        CNBlockConfig(64, 128, 3),
        CNBlockConfig(128, 256, 3),
        CNBlockConfig(256, 512, 27),
        CNBlockConfig(512, None, 3),
    ]
    stochastic_depth_prob = kwargs.pop("stochastic_depth_prob", 0.5)

    model = ConvNeXt(block_setting,output_dim=output_dim, stochastic_depth_prob=stochastic_depth_prob, output_linear_coeffs=output_linear_coeffs, **kwargs)

    if weights is not None:
        model.load_state_dict(weights.get_state_dict(progress=progress, check_hash=True))

    return model


def unext_base(output_dim=16, output_linear_coeffs=[1,0], weights=None, progress: bool = True, **kwargs: Any) -> UNext:
    block_base = 64
    block_setting = [
        CNBlockConfig(block_base, block_base * 2, 3),
        CNBlockConfig(block_base * 2, block_base * 4, 3),
        CNBlockConfig(block_base * 4, block_base * 4, 9),
        CNBlockConfig(block_base * 8, block_base * 2, 3),
        CNBlockConfig(block_base * 4, block_base, 3)
    ]
    stochastic_depth_prob = kwargs.pop("stochastic_depth_prob", 0.5)

    model = UNext(block_setting,output_dim=output_dim, stochastic_depth_prob=stochastic_depth_prob, output_linear_coeffs=output_linear_coeffs, **kwargs)

    if weights is not None:
        model.load_state_dict(weights.get_state_dict(progress=progress, check_hash=True))

    return model