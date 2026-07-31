import numpy as np



def hat(w):
    wx, wy,wz = w
    return np.array([
        [0.0, -wz, wy],
        [wz, 0.0, -wx],
        [-wy, wx, 0.0]
    ])


def vee(H):
    wx,wy,wz = -H[1,2], H[0,2], -H[0,1]
    return np.array([wx,wy,wz])



