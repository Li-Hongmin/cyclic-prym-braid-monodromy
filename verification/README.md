# Braid Monodromy Rigidity on Symmetric Cyclic Prym Slices: Exact Symbolic and Finite Checks

This supplement provides a standalone verification script (`verify_slice_monodromy.py`) and its saved execution log (`verify_slice_monodromy.out`) for the computable algebraic identities, monodromy matrices, ranks, characteristic polynomials, invariant Hermitian form minors, Ramanujan-sum Galois trace averages, length-three cyclotomic separation invariants, method-boundary counterexamples, and determinant-character identities in the manuscript.

## Requirements

- Python 3.10+ (standard library `math`, `cmath`, `fractions`) and `sympy`.
- No external project files, private data, or network access are required.

## Reproduction

From inside the `supplement/` directory, run:

```bash
python3 verify_slice_monodromy.py
```

To verify file integrity against `CONTENTS.sha256`:

```bash
shasum -a 256 -c CONTENTS.sha256
```

## Summary of Checks

1. **Check 1 (Section 2, Lemma 2.1 & Theorem A(1))**: Symbolic verification of the characteristic polynomials of the $4 \times 4$ reduced Gassner puncture generators $g_1(t), g_0(t), g_{-1}(t), g_\infty(t)$ (with $g_\infty(t) = (g_1(t) g_0(t) g_{-1}(t))^{-1}$), the rank-$(2,2)$ bisection identities $(g_0(t) - I_4)(g_0(t) - w I_4) = 0$ and $(g_\infty(t) - (b_2 b_4)^{-1} I_4)(g_\infty(t) - b_\infty I_4) = 0$, and finite verification of $\operatorname{rank}(g_p(t) - c I_4) \ge 2$ ($p \in \{1, 0, -1, \infty\}$) and $\operatorname{rank}(g_1(t)^k - c I_4) \ge 2$ for small conductors $m = 15N$ ($N \in \{7, 11, 13\}$).
2. **Check 2 (Sections 2.3 & 3, Proposition 2.2, Lemma 3.1, Theorem 3.2)**: Symbolic verification of skew-Hermitian symmetry $H(\mathbf{b})^* = -H(\mathbf{b})$, monodromy invariance $g^* J(t) g = J(t)$, the closed-form leading principal minors $\det H(\mathbf{b})_{1..k} = (P_{k+1} - 1)/\prod_{j=1}^{k+1}(b_j - 1)$, the diagonalization $C(t)^{-1} g_1(t) C(t)$, the $12$ off-diagonal entries and $2 \times 2$ minors of $M(t) = C(t)^{-1} g_0(t) C(t)$, and finite verification of $P_{2..5} \neq 1$, $|\operatorname{Spec}(\operatorname{Ad} g_1(t))| = 7$, and the signature distributions across all $23$ representative conductors $N$ in Table 9.1.
3. **Check 3 (Section 4, Lemma 4.1 & Theorem B(1))**: Symbolic verification of the trace and determinant formulas for $W(t) = (g_1(t) g_0(t))^2$, the degree-$8$ and degree-$4$ Chebyshev–Laurent polynomials $M_{\mathrm{prim}}(u)$ and $M_{\mathrm{comp}}(u)$, exact rational Ramanujan-sum Galois means $\mathcal{M}_{\mathrm{prim}}(N)$ and $\mathcal{M}_{\mathrm{comp}}(N)$ across Table 9.1 and all $N \in [7, 300]$ with $\gcd(N, 30) = 1$, and the $N = 7$ companion cyclotomic certificate $|\operatorname{tr}(W(t))|^2 > 16$.
4. **Check 4 (Section 5, Lemma 5.1 & Theorem B(2))**: Symbolic verification of the linearity in $w$ of $\operatorname{tr}(g_1 g_0)$ and $\operatorname{tr}(g_1^2 g_0)$, the Laurent polynomial factorizations in Step 2 of Theorem 5.2, and numerical `cmath` verification of complete separation of all non-trivial primitive cosets $[t] \in \mathcal{T}_{\mathrm{prim}} / \pm 1$ by $W_1 = g_1^2 g_0$ across all $23$ representative conductors $N$ in Table 9.1 on both the primary and companion families.
5. **Check 5 (Section 6, Theorem C / Theorem 6.1)**: Symbolic verification that $\operatorname{tr}(g_1 g_{-1})$ is independent of $b_3$ and satisfies exact conjugate symmetry under $b_1 b_5 = b_2 b_4 = \eta_+$, symbolic verification of the four binomial factorizations (6.6)–(6.8), and numerical `cmath` verification of length-$2$ collisions $I_v(g_1 g_0) = I_{-1}(g_1 g_0)$, $I_v(g_1 g_{-1}) = I_{\pm 1}(g_1 g_{-1})$ and length-$3$ separation $I_v(g_1^2 g_0) \notin \{I_1(g_1^2 g_0), I_{-1}(g_1^2 g_0)\}$ for all $8$ pairs $(a, m')$ in Table 9.2.
6. **Check 6 (Section 7 & Section 8.2, Theorem D(1)–(2) & Corollary E(3))**: Symbolic verification of the pure-tensor Laurent factorizations $F_{--}, F_{+-}, F_{-+}$, symbolic resultant verification of the degree-$6$ minimal polynomial of $|\operatorname{tr}(W(92))|^2$ at $(a, b, N) = (5, 3, 7)$, finite verification of all Theorem D(1) triples in Table 9.2, explicit reproduction of all four boundary counterexamples (b.1)–(b.4) of Theorem D(2) (including exact $\mathfrak{sl}_4$ eigenvalue multiplicities $7, 4, 4$ at $(a, b, N) = (5, 7, 3)$), and exact integer verification of the Corollary E(3) cross-orbit determinant identity $d_{8g} = d_{13g}, d_{8g} + d_{-13g} = 6$ for all $g \in (\mathbf{Z}/105\mathbf{Z})^\times$ on the primary $m = 105, N = 7$ family.
7. **Check 7 (Sections 5–7 & 9, Exact Cyclotomic Ring Certificates)**: 100% floating-point-free exact integer polynomial remainder verification of complete coordinate separation $\tau_t(g_1^2 g_0) \neq \pm \tau_{\pm 1}(g_1^2 g_0)$ in $\mathbf{Z}[z]/(\Phi_N(z))$ across all $23$ conductors $N$ in Table 9.1 (both Primary and Companion families) and in the tensor-product cyclotomic ring $\mathbf{Z}[x_a, x_b, z]/(\Phi_a(x_a), \Phi_b(x_b), \Phi_N(z))$ across all $12$ three-factor triples $(a, b, N)$ in Table 9.2 (together with exact non-vanishing of the individual projected tensor factors).

These symbolic, exact-ring, and finite computations verify the explicit formulas and representative tables in the manuscript; the general theorems for all admissible conductors are established by the mathematical proofs in Sections 2–8 of the paper.
