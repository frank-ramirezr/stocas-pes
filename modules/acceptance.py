import math as m

h2j = 4.35974 * 10 ** (-18)
kb = 1.380649 * 10 ** (-23)


def boltz(actual_temp, e0, e_iter):

    actual_temp = float(actual_temp)
    e0 = float(e0)
    e_iter = float(e_iter)

    boltzmann = m.exp(-h2j * (e_iter - e0) / (kb * actual_temp))

    return boltzmann


def tsallis_accept(actual_temp, e0, e_iter, q=1.2):

    actual_temp = float(actual_temp)
    e0 = float(e0)
    e_iter = float(e_iter)

    delta_e = e_iter - e0
    # Calculate the standard dimensionless energy scaling using your SI factors
    x = (h2j * delta_e) / (kb * actual_temp)

    # Core Tsallis condition
    base = 1.0 - (1.0 - q) * x

    if base <= 0:
        # If the base drops below zero, probability of acceptance is zero
        return 0.0

    p_tsallis = base ** (1.0 / (1.0 - q))
    return p_tsallis