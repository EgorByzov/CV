import pickle
import numpy as np
import torch
from numba import njit, prange


class FocusLensQuadricPrimitive:
    def __init__(self, exit_points,
                 plane_z,
                 n_ref=1.493,
                 use_closest_surface=True):
        self.n_ref = n_ref
        self.use_closest_surface = use_closest_surface
        self.exit_points = exit_points
        self.plane_z = plane_z

    def trace_rays(self, quadric_values, rays):
        return _trace_focus_quads_plane_to_points(self.n_ref, self.plane_z, self.exit_points, self.use_closest_surface,
                                                  quadric_values, rays)

    def update_quadrics(self, quadric_values, cur_step, req_distrs, trace_results):
        en_diff = trace_results - req_distrs
        if not self.use_closest_surface:
            en_diff = -en_diff

        quadric_values = quadric_values + cur_step * en_diff
        return quadric_values

    def trace_rays_batch(self, quadric_values, rays):
        return _trace_focus_quads_plane_to_points_batch(self.n_ref, self.plane_z, self.exit_points,
                                                        self.use_closest_surface,
                                                        quadric_values, rays)

    def update_quadrics_batch(self, quadric_values, cur_step, req_distrs, trace_results):
        return _update_quads_batch(self.use_closest_surface, quadric_values, cur_step, req_distrs, trace_results)

    def to_Oliker(self, quadric_values, req_distrs):
        if self.use_closest_surface:
            quadric_values[req_distrs <= 0] = np.inf
        else:
            quadric_values[req_distrs <= 0] = -np.inf

    def to_CNN(self, quadric_values, req_distrs):
        quadric_values[req_distrs <= 0] = 0

    def save(self, file_name):
        state_dict = {'n_ref': self.n_ref,
                      'use_closest_surface': self.use_closest_surface,
                      'exit_points': self.exit_points,
                      'plane_z': self.plane_z}

        with open(file_name, 'wb') as file:
            pickle.dump(state_dict, file)

    @staticmethod
    def load(file_name):
        with open(file_name, 'rb') as file:
            state_dict = pickle.load(file)

        o = FocusLensQuadricPrimitive(state_dict['exit_points'],
                                      state_dict['plane_z'],
                                      n_ref=state_dict['n_ref'],
                                      use_closest_surface=state_dict['use_closest_surface'])
        return o


####################################### UTIL FUNCS ####################################################

@njit(parallel=True)
def _trace_focus_quads_plane_to_points(n_ref, exit_z, exit_points, use_closest_surface, quadric_values, rays):
    def_shape = quadric_values.shape

    num_rays = rays.shape[0]
    a = n_ref - 1
    b2 = a * (n_ref + 1)

    x_exit = exit_points[:, :, 0]
    y_exit = exit_points[:, :, 1]

    res_t = np.zeros((1, num_rays), dtype=quadric_values.dtype)
    res_illum = np.zeros(def_shape, dtype=quadric_values.dtype)

    for i in prange(num_rays):
        cur_ray = rays[i, :]
        d2 = (cur_ray[0] - x_exit) ** 2 + (cur_ray[1] - y_exit) ** 2
        b1 = a * (exit_z + n_ref * quadric_values)
        b3 = a ** 2 * (exit_z - quadric_values) ** 2

        t = (-np.sqrt(b2 * d2 + b3) + b1) / b2

        if use_closest_surface:
            extremum_ind = t.argmin()
        else:
            extremum_ind = t.argmax()

        extremum_val = t.flat[extremum_ind]
        res_illum.flat[extremum_ind] += cur_ray[-1]
        res_t[0, i] = extremum_val

    return res_illum, res_t


@njit(parallel=True)
def _trace_focus_quads_plane_to_points_batch(n_ref, exit_z, exit_points, use_closest_surface, quadric_values,
                                             rays):
    def_shape = quadric_values.shape
    batch_size = def_shape[0]

    num_rays = rays.shape[0]
    a = n_ref - 1
    b2 = a * (n_ref + 1)

    x_exit = exit_points[:, :, 0]
    y_exit = exit_points[:, :, 1]

    res_t = np.zeros((batch_size, num_rays), dtype=quadric_values.dtype)
    res_illum = np.zeros(def_shape, dtype=quadric_values.dtype)

    for k in prange(batch_size):
        for i in prange(num_rays):
            cur_ray = rays[i, :]
            d2 = (cur_ray[0] - x_exit) ** 2 + (cur_ray[1] - y_exit) ** 2

            cur_quads = quadric_values[k]

            b1 = a * (exit_z + n_ref * cur_quads)
            b3 = a ** 2 * (exit_z - cur_quads) ** 2

            t = (-np.sqrt(b2 * d2 + b3) + b1) / b2

            if use_closest_surface:
                extremum_ind = t.argmin()
            else:
                extremum_ind = t.argmax()

            extremum_val = t.flat[extremum_ind]
            res_illum[k].flat[extremum_ind] += cur_ray[-1]
            res_t[k, i] = extremum_val

    return res_illum, res_t


@njit(parallel=True)
def _update_quads_batch(use_closest_surface, quadric_values, cur_step, req_distrs, trace_results):
    batch_size = quadric_values.shape[0]

    for k in prange(batch_size):
        en_diff = trace_results[k] - req_distrs[k]
        if not use_closest_surface:
            en_diff = -en_diff

        quadric_values[k] = quadric_values[k] + cur_step * en_diff

    return quadric_values
