# SO(3) Module

A minimal implementation of the Lie group **SO(3)** in Python, implemented from first principles rather than copied formulas.

## Features

- `hat(w)` : $\mathbb{R}^3 \rightarrow \mathfrak{so}(3)$
- `vee(K)` : $\mathfrak{so}(3) \rightarrow \mathbb{R}^3$
- `exp_so3(w)` using Rodrigues' formula
- `log_so3(R)` using the matrix logarithm
- Small-angle Taylor approximations
- PyTest test suite
- Validation against SciPy

---

## Repository Structure

```text
so3/
├── so3.py
├── test_so3.py
├── requirements.txt
└── README.md
```

---

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

or with `uv`:

```bash
uv venv
source .venv/bin/activate
uv pip install -r requirements.txt
```

---

## Mathematical Background

SO(3) is the set of rotation matrices:

$$SO(3) = \{ R \in \mathbb{R}^{3 \times 3} \mid RR^T = I, \ \det(R) = 1 \}$$

Its Lie algebra is the set of skew-symmetric matrices:

$$\mathfrak{so}(3) = \{ K \in \mathbb{R}^{3 \times 3} \mid K^T = -K \}$$

---

## Hat Operator

The hat operator converts a vector into a skew-symmetric matrix. Given

$$\omega = \begin{bmatrix} \omega_x \\\ \omega_y \\\ \omega_z \end{bmatrix}$$

define

$$\hat{\omega} = \begin{bmatrix} 0 & -\omega_z & \omega_y \\\ \omega_z & 0 & -\omega_x \\\ -\omega_y & \omega_x & 0 \end{bmatrix}$$

It satisfies

$$\hat{\omega} \, x = \omega \times x$$

The inverse map is

$$\mathrm{vee}(\hat{\omega}) = \omega$$

**Tests:**

```python
assert np.allclose(vee(hat(w)), w)
assert np.allclose(hat(w), -hat(w).T)
```

---

## Exponential Map

We begin with the matrix exponential

$$R = \exp(K) = \sum_{k=0}^{\infty} \frac{K^k}{k!}, \qquad K = \hat{\omega}, \quad \theta = \|\omega\|$$

**Step 1.** Using the vector identity $\omega \times (\omega \times x) = \omega(\omega \cdot x) - x(\omega \cdot \omega)$, we obtain

$$K^2 = \omega \omega^T - \theta^2 I$$

**Step 2.** Then

$$K^3 = K \left( \omega \omega^T - \theta^2 I \right)$$

Since $K\omega = \omega \times \omega = 0$, it follows that

$$K^3 = -\theta^2 K$$

Therefore $K^4 = -\theta^2 K^2$, $K^5 = \theta^4 K$, and every higher power reduces to either $K$ or $K^2$.

**Step 3.** Collect odd powers:

$$\frac{K}{1!} + \frac{K^3}{3!} + \frac{K^5}{5!} + \cdots = \frac{\sin\theta}{\theta} K$$

Collect even powers:

$$\frac{K^2}{2!} + \frac{K^4}{4!} + \frac{K^6}{6!} + \cdots = \frac{1 - \cos\theta}{\theta^2} K^2$$

Therefore:

$$R = I + \frac{\sin\theta}{\theta} K + \frac{1 - \cos\theta}{\theta^2} K^2$$

This is **Rodrigues' rotation formula**.

**Tests:**

- Orthogonality: $RR^T = I$
- Determinant: $\det(R) = 1$
- Compare against `Rotation.from_rotvec(w).as_matrix()`

---

## Logarithm Map

The trace satisfies

$$\operatorname{tr}(R) = 1 + 2\cos\theta$$

Therefore

$$\theta = \arccos\left( \frac{\operatorname{tr}(R) - 1}{2} \right)$$

The antisymmetric part gives

$$R - R^T = \frac{2\sin\theta}{\theta} K$$

Hence

$$K = \frac{\theta}{2\sin\theta} \left( R - R^T \right)$$

Finally

$$\omega = \frac{\theta}{2\sin\theta} \, \mathrm{vee}\left( R - R^T \right)$$

**Tests:**

```python
assert np.allclose(log_so3(exp_so3(w)), w)
assert np.allclose(exp_so3(log_so3(R)), R)
```

for rotations satisfying $\|\omega\| < \pi$.

---

## Small-Angle Approximation

Direct evaluation becomes numerically unstable near zero. Instead use

$$\frac{\sin\theta}{\theta} \approx 1 - \frac{\theta^2}{6} + \frac{\theta^4}{120}$$

and

$$\frac{1 - \cos\theta}{\theta^2} \approx \frac{1}{2} - \frac{\theta^2}{24} + \frac{\theta^4}{720}$$

Similarly, $\dfrac{\theta}{2\sin\theta}$ requires its own Taylor expansion inside `log_so3()`.

---

## Validation

The implementation is validated by:

- `hat` / `vee` inverse tests
- Orthogonality
- Determinant
- Round-trip tests
- SciPy comparison
- Small-angle continuity