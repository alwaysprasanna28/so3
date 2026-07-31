# SO(3) Module — Structure & Math

## Repo structure

```
so3/
├── so3.py             # hat, vee, exp_so3, log_so3
├── test_so3.py         # pytest — one test class per function
├── requirements.txt     # numpy, scipy, pytest
└── README.md            # your derivation notes, small-angle findings
```

Commit after each function has a passing test, not in one batch at the end.

---

## 1. `hat` and `vee`

The map between ℝ³ and 𝔰𝔬(3) (3×3 skew-symmetric matrices), defined so that
`hat(w) @ x == np.cross(w, x)` for all x:

$$
[\omega]_\times = \hat\omega =
\begin{bmatrix}
0 & -\omega_3 & \omega_2 \\
\omega_3 & 0 & -\omega_1 \\
-\omega_2 & \omega_1 & 0
\end{bmatrix}
$$

`vee` is just the inverse: pull the three independent entries back out.

**Test:** `vee(hat(w)) == w`, and `hat(w) == -hat(w).T`.

---

## 2. `exp_so3(w) -> R` — derive, don't paste

Start from the definition of the matrix exponential, with $K = [\omega]_\times$ and $\theta = \|\omega\|$:

$$
R = \exp(K) = \sum_{k=0}^{\infty} \frac{K^k}{k!}
$$

This is an infinite sum of 3×3 matrices — useless as-is. The way out is a **closure identity**: powers of $K$ don't grow into new independent matrices, they cycle between just $K$ and $K^2$.

**Step 1 — find $K^2$.**
Using the vector identity $\omega \times (\omega \times x) = \omega(\omega \cdot x) - x(\omega \cdot \omega)$, applied to arbitrary $x$:

$$
K^2 = \omega\omega^T - \theta^2 I
$$

**Step 2 — find $K^3$.**

$$
K^3 = K \cdot K^2 = K(\omega\omega^T - \theta^2 I) = (K\omega)\omega^T - \theta^2 K
$$

But $K\omega = \omega \times \omega = 0$, so:

$$
K^3 = -\theta^2 K
$$

This is the whole trick. Every power beyond the second collapses back onto $K$ or $K^2$ with a sign flip and a factor of $\theta^2$:

$$
K^1, K^2, K^3=-\theta^2K,\ K^4=-\theta^2K^2,\ K^5=\theta^4K,\ K^6=\theta^4K^2,\ \dots
$$

**Step 3 — split the exponential sum into odd and even powers of $K$.**

Odd terms (all proportional to $K$):
$$
\frac{K}{1!} + \frac{K^3}{3!} + \frac{K^5}{5!} + \cdots
= K\left(1 - \frac{\theta^2}{3!} + \frac{\theta^4}{5!} - \cdots\right)
= K \cdot \frac{1}{\theta}\left(\theta - \frac{\theta^3}{3!} + \frac{\theta^5}{5!} - \cdots\right)
= \frac{\sin\theta}{\theta}\,K
$$

Even terms (all proportional to $K^2$):
$$
\frac{K^2}{2!} + \frac{K^4}{4!} + \frac{K^6}{6!} + \cdots
= K^2\left(\frac{1}{2!} - \frac{\theta^2}{4!} + \frac{\theta^4}{6!} - \cdots\right)
= \frac{1-\cos\theta}{\theta^2}\,K^2
$$

(Recognize the bracketed series as the Taylor series for $\sin\theta/\theta$ and $(1-\cos\theta)/\theta^2$ respectively — that's the whole derivation, just term matching after the $K^3=-\theta^2K$ substitution.)

**Result — Rodrigues' rotation formula:**

$$
\boxed{R = I + \frac{\sin\theta}{\theta}K + \frac{1-\cos\theta}{\theta^2}K^2}
$$

**Tests:**
- $R R^T = I$, $\det R = +1$, for many random $\omega$
- Cross-check every sample against `Rotation.from_rotvec(w).as_matrix()`

---

## 3. `log_so3(R) -> w` — the inverse problem

**Angle**, from the trace. Since $K$ is traceless and $\mathrm{tr}(K^2) = \mathrm{tr}(\omega\omega^T) - 3\theta^2 = -2\theta^2$:

$$
\mathrm{tr}(R) = 3 + 0 + \frac{1-\cos\theta}{\theta^2}(-2\theta^2) = 1 + 2\cos\theta
\quad\Rightarrow\quad
\theta = \arccos\!\left(\frac{\mathrm{tr}(R)-1}{2}\right)
$$

**Axis**, from the antisymmetric part of $R$. $K$ is the only odd-in-$K$ term in Rodrigues' formula, so transposing flips its sign while $I$ and $K^2$ (symmetric) survive unchanged:

$$
R^T = I - \frac{\sin\theta}{\theta}K + \frac{1-\cos\theta}{\theta^2}K^2
\quad\Rightarrow\quad
R - R^T = \frac{2\sin\theta}{\theta}K
$$

Solve for $K$, then `vee` it and scale by $\theta$:

$$
\omega = \frac{\theta}{2\sin\theta}\,\mathrm{vee}(R - R^T)
$$

**Tests:**
- `log_so3(exp_so3(w)) ≈ w` for $\|\omega\| < \pi$ — think about *why* this bound exists: rotations by $\theta$ and by $-\theta$ about the flipped axis are the same $R$ once $\theta$ passes $\pi$, so the inverse stops being single-valued (SO(3) is the closed ball of radius $\pi$ with antipodal boundary points identified — that's the actual geometric content here, not just a numerical caveat).
- `exp_so3(log_so3(R)) ≈ R` for random valid rotations.

---

## 4. Small-angle branch — where it actually breaks

Look at the two coefficients as $\theta \to 0$:

$$
\frac{\sin\theta}{\theta} \to 1, \qquad \frac{1-\cos\theta}{\theta^2} \to \frac12
$$

Both limits are perfectly finite analytically. The failure is purely a float64 problem, and the two terms fail *differently* — work out which before you write the guard:

- $\sin\theta/\theta$ at $\theta=10^{-8}$: numerator and denominator are both tiny but well-scaled: this one is comparatively benign.
- $(1-\cos\theta)/\theta^2$: $\cos(10^{-8})$ rounds to exactly `1.0` in float64 (since $\theta^2/2 \approx 5\times10^{-17}$ is below machine epsilon relative to 1), so the numerator $1-\cos\theta$ becomes an exact `0.0` — total loss of precision from **catastrophic cancellation**, then divided by an already-tiny $\theta^2$.

So it's specifically the *cosine* term that collapses first, and it collapses much earlier (around $\theta \sim 10^{-8}$, roughly $\sqrt{\epsilon_{\text{machine}}}$) than you'd guess from eyeballing "small number near zero."

**Fix**: below some threshold (e.g. $\theta < 10^{-4}$ for safety margin), replace both coefficients with their Taylor series directly, never computing $\sin$, $\cos$, or the division at all:

$$
\frac{\sin\theta}{\theta} \approx 1 - \frac{\theta^2}{6} + \frac{\theta^4}{120}, \qquad
\frac{1-\cos\theta}{\theta^2} \approx \frac12 - \frac{\theta^2}{24} + \frac{\theta^4}{720}
$$

`log_so3` has its own version of this: $\theta/(2\sin\theta)$ has a removable singularity at $\theta=0$ (both hit exact zero), needing the same kind of Taylor guard, and a *separate* numerical issue near $\theta=\pi$ where $\sin\theta \to 0$ again but the axis becomes genuinely ill-conditioned (not just a float issue — small perturbations in $R$ swing the recovered axis a lot). Different regime, different fix (trace-based / Shepperd's-method axis extraction) — worth noting in your README even if you don't implement it yet.

**Test:** sweep $\|\omega\|$ across $10^{-10}, 10^{-8}, 10^{-6}, 10^{-4}, 10^{-2}$ and confirm `R` and the round-trip error are continuous across your threshold — no jump at the switchover.
