import os as os
import math as m
from modules import rot_constants

def rot(rot_deviation):

    cdir = os.getcwd()
    scratch_dir = f'{cdir}/scratch'
    rot_comp = f'{scratch_dir}/rot_comp.dat'

    nlines = 0

    min_max_dev = float('inf')
    best_devs = (0.0, 0.0, 0.0)

    ra2, rb2, rc2 = rot_constants.inertia()

    with open(rot_comp, 'r') as rot_data:

        rconst = []
        nlines = 0
        for line in rot_data:

            nlines = nlines + 1
            rconst.append(line.split())

    for i in range(0, nlines):

        ra1, rb1, rc1 = float(rconst[i][0]), float(rconst[i][1]), float(rconst[i][2])

        ra_dev = 100.0 * m.fabs((ra1 - ra2) / ra1)
        rb_dev = 100.0 * m.fabs((rb1 - rb2) / rb1)
        rc_dev = 100.0 * m.fabs((rc1 - rc2) / rc1)

        if (
            ra_dev < rot_deviation
            and rb_dev < rot_deviation
            and rc_dev < rot_deviation
        ):
            return False, ra_dev, rb_dev, rc_dev

        # Track the closest structure for reporting purposes
        max_dev_curr = max(ra_dev, rb_dev, rc_dev)
        if max_dev_curr < min_max_dev:
            min_max_dev = max_dev_curr
            best_devs = (ra_dev, rb_dev, rc_dev)

# If no match was found across all stored configurations, it is a new unique structure
    return True, best_devs[0], best_devs[1], best_devs[2]