# Numerical method

A short description of what the code computes and where. The master's thesis behind this project derives each piece in detail; the references are listed in the [README](../README.md#references).

## Equations

The 2-D compressible Euler equations in conservation form,

$$
\frac{\partial \mathbf{U}}{\partial t} + \frac{\partial \mathbf{F}(\mathbf{U})}{\partial x} + \frac{\partial \mathbf{G}(\mathbf{U})}{\partial y} = 0,
\qquad
\mathbf{U} = \begin{pmatrix} \rho \\ \rho u \\ \rho v \\ \rho E \end{pmatrix},
$$

closed by the ideal-gas law $p = (\gamma - 1)\left(\rho E - \tfrac12 \rho (u^2 + v^2)\right)$ with $\gamma = 1.4$. The flux through a face with unit normal $\mathbf{n} = (n_x, n_y)$ is

$$
\mathbf{F}_n(\mathbf{U}) = \begin{pmatrix} \rho u_n \\ \rho u\, u_n + p\, n_x \\ \rho v\, u_n + p\, n_y \\ (\rho E + p)\, u_n \end{pmatrix},
\qquad u_n = u n_x + v n_y .
$$

Code: `src/core/physics_eulerphysics.*` (`pressure`, `soundSpeed`, `physicalFluxNormal`).

## Finite-volume discretisation

The mesh is a uniform Cartesian grid of $N_x \times N_y$ cells (`src/core/mesh_structured.*`). Integrating over cell $C$ with area $\Omega_C$ gives the semi-discrete equation

$$
\Omega_C \frac{d\mathbf{U}_C}{dt} = \mathbf{R}_C = -\sum_{f \in \partial C} \hat{\mathbf{F}}(\mathbf{U}_f^-, \mathbf{U}_f^+; \mathbf{n}_f)\, |f| ,
$$

where $\hat{\mathbf{F}}$ is a numerical flux and $\mathbf{U}_f^\mp$ are the states on either side of face $f$.

The residual is assembled in **one loop over faces** (`EulerSolver::computeResiduals`). Each internal face stores its left and right cell and a unit normal pointing from left to right; its flux is subtracted from the left cell and added to the right cell. Every flux therefore leaves one cell and enters its neighbour exactly, so the scheme conserves mass, momentum and energy by construction. Boundary faces carry the outward normal and use a ghost state (see below).

## Numerical fluxes

**Rusanov** (local Lax–Friedrichs):

$$
\hat{\mathbf{F}} = \tfrac12\left(\mathbf{F}_n(\mathbf{U}_L) + \mathbf{F}_n(\mathbf{U}_R)\right) - \tfrac12 S_{\max}\,(\mathbf{U}_R - \mathbf{U}_L),
\qquad S_{\max} = \max\left(|u_{n,L}| + a_L,\ |u_{n,R}| + a_R\right).
$$

Simple and robust, but it smears contact discontinuities because it uses a single wave speed.

**HLLC** (Toro §10.4–10.5) restores the contact wave. The outer wave speeds use the pressure-based estimate

$$
S_L = u_{n,L} - a_L q_L,\quad S_R = u_{n,R} + a_R q_R,\qquad
q_K = \begin{cases} 1 & p^* \le p_K \\ \sqrt{1 + \frac{\gamma+1}{2\gamma}\left(\frac{p^*}{p_K} - 1\right)} & p^* > p_K \end{cases}
$$

with the PVRS estimate $p^* = \max\left(0,\ \tfrac12(p_L + p_R) - \tfrac12(u_{n,R} - u_{n,L})\,\bar\rho\,\bar a\right)$. The contact speed is

$$
S^* = \frac{p_R - p_L + \rho_L u_{n,L}(S_L - u_{n,L}) - \rho_R u_{n,R}(S_R - u_{n,R})}{\rho_L (S_L - u_{n,L}) - \rho_R (S_R - u_{n,R})},
$$

and the flux is $\mathbf{F}_n(\mathbf{U}_L)$, $\mathbf{F}_n(\mathbf{U}_L) + S_L(\mathbf{U}^*_L - \mathbf{U}_L)$, $\mathbf{F}_n(\mathbf{U}_R) + S_R(\mathbf{U}^*_R - \mathbf{U}_R)$ or $\mathbf{F}_n(\mathbf{U}_R)$ depending on where $0$ lies among $S_L \le S^* \le S_R$. The star states keep the tangential velocity of their side. HLLC resolves a stationary contact exactly, which the Shu–Osher results show.

## Reconstruction

**Piecewise constant**: the face states are the cell averages. First order in space.

**MUSCL** with the MinMod limiter (`musclReconstruct` in `physics_solver.cpp`). For a face between cells $L$ and $R$ along one grid direction, with neighbours $LL$ and $RR$,

$$
\mathbf{U}_f^- = \mathbf{U}_L + \tfrac12\,\mathrm{minmod}(\mathbf{U}_L - \mathbf{U}_{LL},\ \mathbf{U}_R - \mathbf{U}_L),
\qquad
\mathbf{U}_f^+ = \mathbf{U}_R - \tfrac12\,\mathrm{minmod}(\mathbf{U}_R - \mathbf{U}_L,\ \mathbf{U}_{RR} - \mathbf{U}_R),
$$

where $\mathrm{minmod}(a, b)$ is zero when $a$ and $b$ have opposite signs and otherwise the one with the smaller magnitude. The limiter is applied to each conserved variable separately. Cells next to a boundary use a zero slope, and boundary faces are always first order. MUSCL is second order on smooth flow when combined with SSP-RK2.

## Time integration

The time step is global:

$$
\Delta t = \mathrm{CFL}\cdot \min_C \frac{\Omega_C}{\sum_{f \in \partial C} \left(|u_n| + a\right)_C |f|} .
$$

The sum runs over **all four faces** of each cell, so this CFL number is about **twice the classical Courant number**: in 1-D, $\mathrm{CFL} = 1$ gives $\Delta t \approx \tfrac12 \Delta x / (|u| + a)$, and the linear stability limit is about $\mathrm{CFL} = 2$.

- **Forward Euler**: $\mathbf{U}^{n+1} = \mathbf{U}^n + \frac{\Delta t}{\Omega} \mathbf{R}(\mathbf{U}^n)$.
- **SSP-RK2** (Heun, Shu–Osher form): $\mathbf{U}^* = \mathbf{U}^n + \frac{\Delta t}{\Omega}\mathbf{R}(\mathbf{U}^n)$, then $\mathbf{U}^{n+1} = \tfrac12\left(\mathbf{U}^n + \mathbf{U}^* + \frac{\Delta t}{\Omega}\mathbf{R}(\mathbf{U}^*)\right)$. A convex combination of Forward Euler steps, so it adds no oscillations that Forward Euler would not produce.

When a test case has an end time, the last step is shortened to land on it exactly.

## Boundary conditions

Each boundary face builds a ghost state $\mathbf{U}_g$ from the interior state $\mathbf{U}_i$ and passes $(\mathbf{U}_i, \mathbf{U}_g)$ to the same numerical flux as an internal face (`src/core/physics_boundarycondition.h`).

| Type | Ghost state |
|---|---|
| Slip wall, symmetry | $\rho$, $p$ and tangential velocity of the interior; normal velocity reversed. Zero mass flux. |
| Supersonic inlet, far-field | The prescribed free-stream state. |
| Supersonic outlet, zero gradient | $\mathbf{U}_g = \mathbf{U}_i$ (waves leave without reflection). |
| Subsonic inlet | Outgoing Riemann invariant $R^+ = u_n + \frac{2a}{\gamma-1}$ from the interior; incoming invariant $R^-$, entropy $p/\rho^\gamma$ and tangential velocity from the free stream. |
| Subsonic outlet | Back pressure $p_\text{back}$ prescribed; entropy, tangential velocity and $R^+$ from the interior. |

The subsonic inlet fixes the free stream's static entropy and $R^-$, not its total pressure, so a back pressure below the free-stream pressure accelerates the flow (see the channel study).

## Robustness

- After every update, and after the first stage of SSP-RK2, the solver checks that every cell has positive, finite density and pressure (`EulerSolver::checkStates`). If not, it stops with an error that names the cell, step and time. The GUI shows the message, and the CLI exits with code 2.
- MUSCL face states keep a positive density (MinMod bounds them by the neighbouring cell values) but can have $p \le 0$. Both fluxes floor the pressure at $10^{-10}$ before square roots and divisions, so such a face still gives a finite flux.

## Exact Riemann solver

The reference for Sod (`src/core/sod_exact.*`, and `analysis/euler_analysis/exact.py` in Python) solves $f_L(p^*) + f_R(p^*) + u_R - u_L = 0$ by Newton iteration from the PVRS guess (Toro ch. 4) and samples the self-similar solution in $\xi = (x - x_0)/t$: shock or rarefaction on each side, and the contact at $u^*$. For Sod it gives $p^* = 0.30313$ and $u^* = 0.92745$.
