import math as m

def linear(t0, tf, n_temps, iteration):

    t0 = float(t0)
    tf = float(tf)
    n_temps = float(n_temps)
    iteration = float(iteration)

    actual_temp = t0 - (t0 - tf) * (iteration / n_temps)
    
    return actual_temp

def exponential(t0, tf, n_temps, iteration):

    t0 = float(t0)
    tf = float(tf)
    n_temps = float(n_temps)
    iteration = float(iteration)

    actual_temp = t0 * (tf / t0) ** (iteration / n_temps)
    
    return actual_temp

def decay(t0, tf, n_temps, iteration):

    t0 = float(t0)
    tf = float(tf)
    n_temps = float(n_temps)
    iteration = float(iteration)

    beta = (n_temps + 1) / n_temps

    actual_temp = t0 - (t0 - tf) * beta * iteration / (iteration + 1.0)

    return actual_temp