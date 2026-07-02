# Numerical experiments for *Toward a Higher-Dimensional Poincaré–Hopf Theory for Filippov Vector Fields*

This repository contains the MATLAB code that generates the figures and reproduces
the numerical checks referenced in the paper

> **Toward a higher-dimensional Poincaré–Hopf theory for Filippov vector fields**
> Sri Venkata Durga Sudarsan Madhyannapu.

The code is written in **plain MATLAB** and needs **no toolboxes** (tested on MATLAB Online).

## How to run

```matlab
cd matlab
run_all
```

Each script prints clearly labelled output values and saves its figure(s) to `figures/`
as both `.pdf` (vector) and `.png`.

## What each script does (and where it lives in the paper)

| Script | Checks | Paper reference | Status |
|---|---|---|---|
| `s1_wedge_identity.m` | `‖u∧v‖² + ⟨u,v⟩² = ‖u‖²‖v‖²` | eq. (wedge), Def. (angmag) | **proved** — sanity check |
| `s2_planar_index_cmn.m` | planar local index is an **integer** and equals the Sotomayor–Teixeira regularized index | Def. (planarindex), Thm (cmn) | **proved** (Casimiro–Martins–Novaes) |
| `s3_gauss_degree_S2.m` | Gauss-map degree normalization on `S²` (identity → +1, `z↦z²` → +2, antipodal → −1) | eq. (W), `vol(Sⁿ⁻¹)` normalization | **proved** — sanity check |
| `s4_candidate_index_n3.m` | candidate index in `n = 3` vs regularized degree, on specific linear examples | eq. (interp), (J), (index) | **numerical evidence only** |
| `s5_delta_independence.m` | candidate index is (numerically) independent of the tube radius `δ` | Conjecture (deltaindep) | **numerical evidence only** |
| `s6_proximity_condition.m` | `‖a−b‖ + δ(‖a′‖+‖b′‖) < μ` guarantees a nonvanishing straight-line homotopy | eq. (proximity), Prop. (condhom) | **proved** — illustration |

### ⚠️ Scientific-honesty note

For dimension `n ≥ 3`, the paper **does not claim** that the candidate local index is an
integer or that it is invariant under regularization — these are stated as **open
conjectures**. Scripts `s4` and `s5` therefore provide **numerical evidence on particular
examples only**; they are *not* proofs and must not be described as such in the paper.

## Reproducing the comparison in the paper

Run `run_all`, then compare the printed values against the paper:

- **S2** should show `I_index = I_reg` (an integer) for every generic example.
- **S3** should return `+1`, `−1`, `+2` (up to ~1e-2 discretization error).
- **S4** should show `candidate ≈ regularized`, both close to an integer.

(Small deviations of order `1e-2`–`1e-3` are the finite-grid discretization error of the
surface-degree integral; increase `Nu, Nv` in `deg_sphere.m` to reduce them.)

## Figures produced

`fig_s1_wedge`, `fig_s2_planar`, `fig_s3_degree`, `fig_s4_tube`, `fig_s5_delta`, `fig_s6_proximity`
(each as `.pdf` and `.png`). Drop-in LaTeX `figure` blocks are in `latex/figure_snippets.tex`.

## Repository layout

```
.
├── README.md
├── LICENSE
├── CITATION.cff
├── matlab/
│   ├── run_all.m
│   ├── s1_wedge_identity.m
│   ├── s2_planar_index_cmn.m
│   ├── s3_gauss_degree_S2.m
│   ├── s4_candidate_index_n3.m
│   ├── s5_delta_independence.m
│   ├── s6_proximity_condition.m
│   ├── deg_sphere.m          % numerical Brouwer degree on S^2
│   ├── assembled_field.m     % caps + equatorial-tube interpolation (eq. interp)
│   ├── reg_field.m           % Sotomayor–Teixeira regularization (eq. reg)
│   └── exportfig.m           % figure export helper
├── latex/
│   ├── data_availability.tex % the Data-availability statement (already in the paper)
│   └── figure_snippets.tex   % drop-in \begin{figure} blocks
└── figures/                  % output (created on first run)
```

## Citing

See `CITATION.cff`. Please also cite Casimiro–Martins–Novaes (*Commun. Pure Appl. Anal.*
**23** (2024), no. 11, 1770–1796) for the planar theory that this work builds on.
