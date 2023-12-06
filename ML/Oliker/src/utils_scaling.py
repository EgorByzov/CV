def linear_scaling(values, new_interval):
    old_interval = [values.min(), values.max()]

    nA = (new_interval[1] - new_interval[0]) / (old_interval[1] - old_interval[0])
    nB = new_interval[0] - old_interval[0] * nA

    linear_coeffs = [nA, nB]

    values = linear_coeffs[0] * values + linear_coeffs[1]

    return values, linear_coeffs, old_interval


def linear_inverse_scaling(values, linear_coeffs):
    return (values - linear_coeffs[1]) / linear_coeffs[0]
