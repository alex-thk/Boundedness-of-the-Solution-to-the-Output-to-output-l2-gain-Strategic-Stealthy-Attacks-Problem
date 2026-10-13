# On the Boundedness of the Solution to the Output-to-output ℓ2-gain Strategic Stealthy Attacks Problem

MATLAB code accompanying the IEEE TAC paper, which characterizes, via an infeasibility certificate, when the output-to-output ℓ2-gain of a stealthy attack is unbounded.

## What the code reproduces

| Script | Paper |
|---|---|
| `quadruple_tank_example.m` | **Example 1** quadruple-tank process - system matrices in Eq. (21): construction of the certificate $Y$, the solver-based certificate, infeasibility of Problem 3, and the data of **Fig. 2** (cumulative energies of $y_p$ and $y_r$). |
| `toy_example_perturbation_and_epsilon.m` | **Example 2** shared unstable zero - perturbation study of **Fig. 3** and convergence of $\gamma_3(ε_k)$ to $\gamma_2^\star$ (**Remark 1**). |

`psi_map.m` builds the bijection mapping in Eq. (1), `get_Uij.m` builds the corresponding $U_{ij}$ matrix for given indiced and `get_F3_block` builds $[F_{\psi(i,j)}]_3$. 

## Requirements

- MATLAB with the Control System Toolbox
- [YALMIP](https://yalmip.github.io/)
- [MOSEK](https://www.mosek.com/) 

## How to run

1. Add this folder (and YALMIP/MOSEK) to the MATLAB path.
2. Run `quadruple_tank_example.m` and `toy_example_perturbation_and_epsilon.m`.
3. The generated CSV files are saved in `data/`.

In `toy_example_perturbation_and_epsilon.m`, the random seed is fixed *i.e.* (`rng(1)`), so Fig. 3 is reproducible.

## Paper

Submitted for review at IEEE Transactions on Automatic Control
