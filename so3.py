import numpy as np



def hat(w):
    wx, wy,wz = w
    return np.array([
        [0.0, -wz, wy],
        [wz, 0.0, -wx],
        [-wy, wx, 0.0]
    ])


