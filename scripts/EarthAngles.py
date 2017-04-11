import numpy as np
from numpy import cos, sin, pi


def calcGamma(alpha, theta, t):
    w = 2*pi/24.0;
    nhat = np.array([cos(theta)*cos(w*t), cos(theta)*sin(w*t), sin(theta)])
    v_chi = np.array([-sin(alpha), 0, -cos(alpha)])
    return np.arccos(np.sum(np.multiply(nhat, v_chi)))
    
def getLatitude(loc):
    angle = 0
    if (loc == "LNGS"):
        angle = 42.5
    elif (loc == "CJPL"):
        angle = 28.15
    elif (loc == "SUPL"):
        angle = -37.1
    elif (loc == "EQUATOR"):
        angle = 0.0
    elif (loc == "SOUTH POLE"):
        angle = -90.0
    elif (loc == "NORTH POLE"):
        angle = 90.0
    elif (loc == "INO"):
        angle = 9.7
    else:
        angle = 0.0
    
    return angle*pi/180.0
    