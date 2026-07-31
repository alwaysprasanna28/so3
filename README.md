# SO(3) Module

A lightweight Python implementation of the Lie group **SO(3)** and its Lie algebra **𝔰𝔬(3)**, implemented from first principles.

The goal of this project is **to derive the mathematics before writing the code**. Rather than copying Rodrigues' formula or the logarithm map from a textbook, each implementation follows directly from the matrix exponential and the algebraic properties of skew-symmetric matrices.

---

## Features

- `hat(w)` — map from ℝ³ → 𝔰𝔬(3)
- `vee(K)` — inverse map 𝔰𝔬(3) → ℝ³
- `exp_so3(w)` — exponential map (Rodrigues' formula)
- `log_so3(R)` — logarithm map
- Numerical small-angle handling
- Pytest unit tests
- Validation against SciPy

---

## Repository Structure

```text
so3/
├── so3.py             # hat, vee, exp_so3, log_so3
├── test_so3.py        # pytest tests
├── requirements.txt
└── README.md
```

Each function was implemented independently and committed only after all associated tests passed.

---

# Mathematical Background

SO(3) is the Lie group of 3D rotation matrices,

\[
SO(3)=
\{R\in\mathbb{R}^{3\times3}\mid
RR^T=I,\;
\det(R)=1
\}
\]

Its Lie algebra is the set of skew-symmetric matrices,

\[
\mathfrak{so}(3)=
\{
K\in\mathbb{R}^{3\times3}
\mid
K^T=-K
\}
\]

The exponential map connects the Lie algebra to the Lie group,

\[
R=\exp(K).
\]

---

# 1. Hat and Vee Operators

The **hat operator** converts a vector into its skew-symmetric matrix,

\[
\omega=
\begin{bmatrix}
\omega_x\\
\omega_y\\
\omega_z
\end{bmatrix}
\]

\[
[\omega]_\times=
\hat\omega=
\begin{bmatrix}
0 & -\omega_z & \omega_y\\
\omega_z & 0 & -\omega_x\\
-\omega_y & \omega_x & 0
\end{bmatrix}.
\]

It satisfies

\[
\hat\omega x
=
\omega\times x
\]

for every vector \(x\).

The inverse mapping is the **vee operator**

\[
\mathrm{vee}(\hat\omega)=\omega.
\]

---

## Tests

- `vee(hat(w)) == w`
- `hat(w) == -hat(w).T`

---

# 2. Exponential Map

We begin from the matrix exponential,

\[
R
=
\exp(K)
=
\sum_{k=0}^{\infty}
\frac{K^k}{k!},
\]

where

\[
K=[\omega]_\times,
\qquad
\theta=\|\omega\|.
\]

Computing infinitely many powers directly is impractical.

The key observation is that powers of a skew-symmetric matrix are **not independent**.

---

## Computing \(K^2\)

Using

\[
\omega\times(\omega\times x)
=
\omega(\omega\cdot x)
-
x(\omega\cdot\omega),
\]

we obtain

\[
K^2
=
\omega\omega^T
-
\theta^2I.
\]

---

## Computing \(K^3\)

\[
K^3
=
K(\omega\omega^T-\theta^2I)
=
(K\omega)\omega^T
-
\theta^2K.
\]

Since

\[
K\omega
=
\omega\times\omega
=
0,
\]

it follows that

\[
K^3
=
-\theta^2K.
\]

This identity implies

\[
K^4=-\theta^2K^2,
\]

\[
K^5=\theta^4K,
\]

\[
K^6=\theta^4K^2,
\]

and so on.

Every power collapses back to either \(K\) or \(K^2\).

---

## Splitting the exponential

Odd powers become

\[
K
\left(
1
-
\frac{\theta^2}{3!}
+
\frac{\theta^4}{5!}
-
\cdots
\right)
=
\frac{\sin\theta}{\theta}K.
\]

Even powers become

\[
K^2
\left(
\frac1{2!}
-
\frac{\theta^2}{4!}
+
\frac{\theta^4}{6!}
-
\cdots
\right)
=
\frac{1-\cos\theta}{\theta^2}K^2.
\]

---

## Rodrigues' Formula

Combining both series,

\[
\boxed{
R
=
I
+
\frac{\sin\theta}{\theta}K
+
\frac{1-\cos\theta}{\theta^2}K^2
}
\]

---

## Tests

- \(RR^T=I\)
- \(\det(R)=1\)
- Compare against

```python
Rotation.from_rotvec(w).as_matrix()
```

for many random rotation vectors.

---

# 3. Logarithm Map

The logarithm map inverts Rodrigues' formula.

---

## Recovering the angle

Since

\[
\operatorname{tr}(K)=0,
\]

and

\[
\operatorname{tr}(K^2)
=
-2\theta^2,
\]

Rodrigues' formula gives

\[
\operatorname{tr}(R)
=
1
+
2\cos\theta.
\]

Therefore,

\[
\boxed{
\theta
=
\arccos
\left(
\frac{\operatorname{tr}(R)-1}{2}
\right)
}
\]

---

## Recovering the axis

Transpose Rodrigues' formula,

\[
R^T
=
I
-
\frac{\sin\theta}{\theta}K
+
\frac{1-\cos\theta}{\theta^2}K^2.
\]

Subtract,

\[
R-R^T
=
\frac{2\sin\theta}{\theta}K.
\]

Hence,

\[
K
=
\frac{\theta}{2\sin\theta}
(R-R^T).
\]

Finally,

\[
\boxed{
\omega
=
\frac{\theta}{2\sin\theta}
\operatorname{vee}(R-R^T)
}
\]

---

## Tests

Round-trip:

```python
log_so3(exp_so3(w))
```

should recover

```python
w
```

for

\[
\|\omega\|<\pi.
\]

Likewise,

```python
exp_so3(log_so3(R))
```

should recover the original rotation matrix.

The restriction

\[
\|\omega\|<\pi
\]

exists because rotations by

\[
(\theta,\hat u)
\]

and

\[
(-\theta,-\hat u)
\]

represent the same element of SO(3) once the rotation angle exceeds π. Consequently, the logarithm is no longer single-valued.

---

# Small-Angle Numerical Stability

The exponential map contains

\[
\frac{\sin\theta}{\theta},
\qquad
\frac{1-\cos\theta}{\theta^2}.
\]

Analytically,

\[
\frac{\sin\theta}{\theta}\rightarrow1,
\]

\[
\frac{1-\cos\theta}{\theta^2}\rightarrow\frac12.
\]

The problem is numerical.

For

\[
\theta\approx10^{-8},
\]

double-precision arithmetic rounds

```python
np.cos(theta)
```

to exactly

```python
1.0
```

making

```python
1 - np.cos(theta)
```

equal to zero.

This is catastrophic cancellation.

---

## Taylor Approximation

For sufficiently small angles,

\[
\frac{\sin\theta}{\theta}
\approx
1
-
\frac{\theta^2}{6}
+
\frac{\theta^4}{120},
\]

\[
\frac{1-\cos\theta}{\theta^2}
\approx
\frac12
-
\frac{\theta^2}{24}
+
\frac{\theta^4}{720}.
\]

The implementation switches to these series below a threshold (e.g. \(10^{-4}\)).

Similarly,

\[
\frac{\theta}{2\sin\theta}
\]

in the logarithm map has a removable singularity at zero and is handled with its own Taylor expansion.

A separate numerical issue occurs near

\[
\theta=\pi,
\]

where the recovered axis becomes ill-conditioned because

\[
\sin\theta\rightarrow0.
\]

Handling this robustly typically requires specialized axis extraction methods (e.g. Shepperd's algorithm), which are outside the scope of this implementation.

---

# Validation

The implementation is validated by:

- Analytical identities
- Orthogonality checks
- Determinant checks
- Round-trip consistency
- Comparison against SciPy's `Rotation` implementation
- Small-angle continuity tests across

```
1e-10
1e-8
1e-6
1e-4
1e-2
```

ensuring no discontinuity at the Taylor-series threshold.

---

# References

- Ethan Eade, *Lie Groups for 2D and 3D Transformations*
- Timothy D. Barfoot, *State Estimation for Robotics*
- Joan Solà, *Quaternion Kinematics for the Error-State Kalman Filter*
- Sophus Lie group formulations for SO(3)