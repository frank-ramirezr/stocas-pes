def linear(t0, tf, n_temps, iteration):

    t0 = float(t0)
    tf = float(tf)
    n_temps = float(n_temps)
    iteration = float(iteration)

    actual_temp = t0 - (t0 - tf) * (iteration / n_temps)
    return actual_temp

