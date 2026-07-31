import numpy as np
from so3 import hat,vee, exp_so3, log_so3
from scipy.spatial.transform import  Rotation

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


class TestExpSO3:
    def test_zero_rotation(self):
        R = exp_so3(np.zeros(3))

        assert np.allclose(R, np.eye(3))

    def test_orthogonal(self):
        w = np.array([0.3, -0.5, 1.2])
        R = exp_so3(w)

        assert np.allclose(
            R @ R.T,
            np.eye(3),
            atol = 1e-10
        )

    def test_determinant(self):
        w = np.array([0.3,0.7,-0.345])
        R = exp_so3(w)
        assert np.allclose(
            np.linalg.det(R),
            1.0,
        )

    def test_against_scipy(self):
        rng = np.random.default_rng(0)

        for _ in range(100):
            w = rng.normal(size=3)

            R_expected = Rotation.from_rotvec(w).as_matrix()
            R_actual  = exp_so3(w)   

            assert np.allclose(
                R_actual,
                R_expected,
                atol = 1e-10
            )
            



class TestLogSO3:

    def test_identity(self):
        w = np.array([0.3,0.4,-0.345])
        R = exp_so3(w)
        w2 = log_so3(R)

        assert np.allclose(
            w2,
            w,
            atol = 1e-10
        )