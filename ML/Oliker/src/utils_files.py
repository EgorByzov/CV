import os
import pickle

import numpy as np


def save_quads(filepath, quads, req_distrs, rays, params, loss=[], time_total=0, rrmse=[], optim_params=[],
               tracing_results=[]):
    design_logs = {'loss': loss,
                   'time_total': time_total,
                   'rrmse': rrmse,
                   'optim_params': optim_params,
                   'tracing_results': tracing_results}

    project = {'params': params,
               'rays': rays,
               'quadrics': quads,
               'req_distrs': req_distrs,
               'design_logs': design_logs}

    with open(filepath, 'wb') as file:
        pickle.dump(project, file)


def load_quads(filepath):
    quadrics = []
    req_distrs = []
    rays = []
    params = []
    design_logs = []

    with open(filepath, 'rb') as file:
        project = pickle.load(file)

        quadrics = project['quadrics']
        req_distrs = project['req_distrs']
        rays = project['rays']
        params = project['params']
        design_logs = project['design_logs']

    return params, quadrics, req_distrs, rays, design_logs


def unite_quads(new_file_name, file_names=[], folder_name=[]):
    if os.path.isdir(folder_name):
        with os.scandir(folder_name) as it:
            for entry in it:
                if entry.is_file:
                    file_names.append(entry.path)

    params = []
    for i in range(len(file_names)):
        params_cur, quadrics_cur, req_distrs_cur, rays_cur, design_logs_cur = load_quads(entry)

        if i == 0:
            params = params_cur
            rays = rays_cur
            quads = quadrics_cur
            req_distrs = req_distrs_cur

            loss = design_logs_cur['loss']
            rrmse = design_logs_cur['rrmse']
            optim_params = design_logs_cur['optim_params']
            time_total = design_logs_cur['time_total']
            tracing_results = design_logs_cur['tracing_results']
        else:
            quads = np.concatenate((quads, quadrics_cur), axis=0)
            req_distrs = np.concatenate((req_distrs, req_distrs_cur), axis=0)

            loss = np.concatenate((loss, design_logs_cur['loss']), axis=0)
            rrmse = np.concatenate((rrmse, design_logs_cur['rrmse']), axis=0)
            optim_params = np.concatenate((optim_params, design_logs_cur['optim_params']), axis=0)
            time_total += design_logs_cur['time_total']
            tracing_results = np.concatenate((tracing_results, design_logs_cur['tracing_results']), axis=0)

    if params:
        save_quads(new_file_name, quads, req_distrs, rays, params,
                   loss=loss, time_total=time_total, rrmse=rrmse,
                   optim_params=optim_params, tracing_results=tracing_results)

        print(f'Total number of samples = {rrmse.shape[0]}')
