import math
import numpy as np
import logging

logger = logging.getLogger(__name__)


def angle3Dchain(liste):
    logger.debug(f"liste {liste}")
    X0 = np.asarray(liste[-2])
    X1 = np.asarray(liste[-1])
    [x, y, z] = X1-X0
    
    logger.debug(f"X0:{X0}; X1:{X1}; [x,y,z]: {[x,y,z]}")

    length = (x*x + y*y + z*z)**0.5
    phi = math.asin(z/length)

    if y >= 0 and (x*x+y*y != 0):
        theta = math.acos(x/(x*x+y*y)**0.5)
    elif y >= 0 and (x*x+y*y == 0):
        theta = 0.0
    else:
        theta = 2*math.pi-math.acos(x/(x*x+y*y)**0.5)
    return theta, phi
