import random as rn
import math as m

precision = 12

def trasl_matrix(max_trasl):

    phi = rn.uniform(a=0.0, b=2.0 * m.pi)
    costheta = rn.uniform(a=-1.0, b=1.0)
    theta = m.acos(costheta)

    r = max_trasl * rn.uniform(a=0.0, b=1.0)

    dx = r * m.sin(theta) * m.cos(phi)
    dy = r * m.sin(theta) * m.sin(phi)
    dz = r * m.cos(theta)

    return dx, dy, dz