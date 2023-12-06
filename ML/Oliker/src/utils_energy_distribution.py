from numba import *
import numpy as np

@njit
def rrmse_loss(target_distr, current_distr):
    result = (target_distr - current_distr) ** 2
    result = np.sqrt(np.mean(result)) * 100 / np.mean(target_distr)
    return result


@njit(parallel=True)
def rrmse_loss_batch(target_distr, current_distr):
    batch_size = target_distr.shape[0]
    result = np.zeros((batch_size))

    for k in prange(batch_size):
        result[k] = rrmse_loss(target_distr[k], current_distr[k])

    return result


def masked_rrmse_loss(target_distr, current_distr):
    cur_mask = target_distr > 0
    result = (target_distr - current_distr) ** 2
    result = result[cur_mask]
    target_distr = target_distr[cur_mask]
    result = np.sqrt(np.mean(result)) * 100 / np.mean(target_distr)
    return result


def masked_rrmse_loss_batch(target_distr, current_distr):
    batch_size = target_distr.shape[0]
    result = np.zeros((batch_size))

    for k in range(batch_size):
        result[k] = masked_rrmse_loss(target_distr[k], current_distr[k])

    return result


def random_distribution(side_num_pnts, np_min, np_max, total_flux=1, batch_size=1, dtype=np.double):
    rng = np.random.default_rng()
    are_pnts_sparced = rng.choice([False, True], batch_size)

    distributions = np.zeros((batch_size, side_num_pnts, side_num_pnts), dtype=dtype)
    inds = np.arange(side_num_pnts ** 2)

    for i in range(batch_size):
        num_points = rng.choice(np_max - np_min + 1) + np_min
        en_value = total_flux / num_points

        if are_pnts_sparced[i]:
            rng.shuffle(inds)
            cur_inds = inds[0:num_points]
            (rows, cols) = np.unravel_index(cur_inds, (side_num_pnts, side_num_pnts))
            distributions[i, rows, cols] = en_value
        else:
            start_ind = rng.choice(inds)
            (start_row, start_col) = np.unravel_index(start_ind, (side_num_pnts, side_num_pnts))

            field_size = int(np.ceil(np.sqrt(num_points)))

            start_row -= int(np.floor(field_size / 2))
            start_col -= int(np.floor(field_size / 2))

            start_row = 0 if start_row < 0 else start_row
            start_col = 0 if start_col < 0 else start_col

            start_row = side_num_pnts - field_size if start_row + field_size > side_num_pnts else start_row
            start_col = side_num_pnts - field_size if start_col + field_size > side_num_pnts else start_col

            cur_inds = np.arange(field_size ** 2)
            rng.shuffle(cur_inds)
            cur_inds = cur_inds[0:num_points]
            (rows, cols) = np.unravel_index(cur_inds, (field_size, field_size))
            distributions[i, rows + start_row, cols + start_col] = en_value

    return distributions
