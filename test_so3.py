import numpy as np
from so3 import hat


class TestHat:

    def test_skew_symmetric(self):
        w = np.array([1.2, -0.4, 3.6])

        H = hat(w)

        assert np.allclose(H, -H.T)
        


