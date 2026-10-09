import random as rn
import numpy as np
import math as m


def rot_matrix(max_rot=0.1):
    zero = 0.0
    one = 1.0

    alpha = rn.uniform(a=-max_rot, b=max_rot)
    beta = rn.uniform(a=-max_rot, b=max_rot)
    gamma = rn.uniform(a=-max_rot, b=max_rot)

    matx = [
        [one, zero, zero],
        [zero, m.cos(alpha), -m.sin(alpha)],
        [zero, m.sin(alpha), m.cos(alpha)]
    ]

    maty = [
        [m.cos(beta), zero, m.sin(beta)],
        [zero, one, zero],
        [-m.sin(beta), zero, m.cos(beta)]
    ]

    matz = [
        [m.cos(gamma), -m.sin(gamma), zero],
        [m.sin(gamma), m.cos(gamma), zero],
        [zero, zero, one]
    ]

    mat_xy = np.matmul(matx, maty)
    mat_xyz = np.matmul(mat_xy, matz)

    rot11, rot12, rot13 = np.asarray(mat_xyz[0], dtype='float')
    rot21, rot22, rot23 = np.asarray(mat_xyz[1], dtype='float')
    rot31, rot32, rot33 = np.asarray(mat_xyz[2], dtype='float')

    return rot11, rot12, rot13, rot21, rot22, rot23, rot31, rot32, rot33

