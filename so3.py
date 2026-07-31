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


def exp_so3(w):
    theta = np.linalg.norm(w)
    K  = hat(w)
    K2 = K @ K
    I = np.eye(3)
    if theta < 1e-4:
        t2 = theta * theta  

        A = 1 - t2 / 6 + t2 * t2 / 120
        B = 0.5 - t2 / 24 + t2 * t2 / 720
    else:
        A = np.sin(theta) /theta
        B = (1 - np.cos(theta)) / theta**2

    return I + A * K + B * K2


