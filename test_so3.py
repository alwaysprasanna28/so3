import numpy as np
from so3 import hat,vee


class TestHat:

    def test_skew_symmetric(self):
        w = np.array([1.2, -0.4, 3.6])

        H = hat(w)

        assert np.allclose(H, -H.T)


class TestVee:
    def test_vee_hat_inverse(self):
        w = np.array([2.0,4.0,-5])

        assert np.allclose(
            vee(hat(w)), w
        )

    def test_hat_vee_inverse(self):
        H = np.array([
                [0.0, -3.0, 2.0],
                [3.0, 0.0, -1.0],
                [-2.0, 1.0, 0.0],
            ])
    
        assert np.allclose(
                hat(vee(H)), H
            )
    

