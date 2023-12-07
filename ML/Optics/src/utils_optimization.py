import math
import time

import torch
from tqdm import tqdm

from src.utils_energy_distribution import *


def numpy2torch_batch(np_batch, device):
    torch_dtype = torch.float32 if np_batch.dtype == np.single else torch.double
    return torch.tensor(np.expand_dims(np_batch, axis=1), dtype=torch_dtype, device=device)


def torch2numpy_batch(torch_batch):
    return torch_batch.squeeze().cpu().detach().numpy()


def run_oliker_optimization(quad_primitive, req_distrs, rays,
                            loss_func=rrmse_loss_batch, step=1, num_iters=5, end_step=1e-3,
                            quadrics=1, fixed_quads=0):
    dtype = rays.dtype
    batch_size = req_distrs.shape[0]
    # num_saved_iters = 2
    # iters_to_save = math.ceil(num_iters / (num_saved_iters + 1))

    k = math.log(step / end_step) / (num_iters - 1)
    a = step * math.pow(step / end_step, 1 / (num_iters - 1))

    if np.isscalar(quadrics):
        quadrics = np.ones(req_distrs.shape, dtype=dtype) * quadrics
    elif req_distrs.shape != quadrics.shape:
        raise Exception("quadrics shape must be the same as req_distrs, or it should be scalar value")

    quad_primitive.to_Oliker(quadrics, req_distrs)
    # if quad_primitive.use_closest_surface:
    #     quadrics[req_distrs <= 0] = np.inf
    # else:
    #     quadrics[req_distrs <= 0] = -np.inf

    is_fixed_quad = False
    if not np.isscalar(fixed_quads):
        is_fixed_quad = True
        fixed_quads_values = quadrics[fixed_quads]
        # optimisable_quads = req_distrs > 0

    loss = np.zeros((batch_size, num_iters), dtype=dtype)
    saved_quads = []
    time_start = time.time()

    for i in tqdm(range(num_iters)):
        trace_results, _ = quad_primitive.trace_rays_batch(quadrics, rays)

        loss[:, i] = loss_func(req_distrs, trace_results)

        cur_step = a * math.exp(-k * i)
        quadrics = quad_primitive.update_quadrics_batch(quadrics, cur_step, req_distrs, trace_results)

        if is_fixed_quad:
            # quads_diffs = (quadrics[fixed_quads] - fixed_quads_values)/64
            # for k in range(batch_size):
            #     quadrics[k][optimisable_quads[k]] += -quads_diffs[k]
            quadrics[fixed_quads] = fixed_quads_values

        # if ((i + 1) % iters_to_save == 0) & ((i + 1) < num_iters):
        #     saved_quads.append((i, quadrics.copy()))

    time_end = time.time()
    time_total = time_end - time_start

    print("Total time: " + str(time_total / num_iters) + " s/iter")
    return quadrics, loss, time_total, saved_quads
