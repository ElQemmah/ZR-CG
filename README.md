# ZR-CG: Zonotopic Robust Data-Driven Command Governor

MATLAB implementation of the **Zonotopic Robust Data-Driven Command Governor (ZR-CG)** for unknown pre-controlled linear systems subject to zonotopically bounded disturbances.
This repository contains the code reproducing the numerical examples of the paper, including the data-driven identification of the matrix zonotope, its interval-hull over-approximation, the tractable zonotopic reachability propagation, and the online convex QP.

---

### Dependencies
- MATLAB (R2021b or later recommended).
- [CORA 2022 Toolbox](https://tumcps.github.io/CORA/) — zonotope, matrix zonotope and interval matrix classes.
- [CoGoV Toolbox](https://github.com/vinz-uts/CoGoV/) — useful for the main ZR-CG class.
- [YALMIP](https://yalmip.github.io/) — modelling interface for the online QP and the offline LMI.
- A QP/SDP solver, e.g. [Gurobi](https://www.gurobi.com) (used in the paper).


## Features
- **Identification**: computation of the matrix zonotope of all closed-loop models consistent with noisy state/command data (Lemma 1) and of its interval hull.
- **Stability certificate**: offline tiered verification of the robust Schur stability of the state-transition block of the interval hull (majorant test and vertex-based LMI).
- **Reachability**: tractable forward propagation of zonotopic reachable sets through the interval-hull/zonotope product (Lemma 2), with linear generator growth.
- **Supervision**: online Command Governor as a convex QP, with epigraph reformulation of the robust constraints and terminal constraint tightening.
- **Conservatism analysis**: support-function gaps and volume ratios against the exact matrix-zonotope propagation and against a per-step order-reduced pipeline.

## Note
The exact matrix-zonotope propagation implemented in `exact_paper_enclosure.m` is provided for comparison purposes only. As discussed in the paper, the number of generators grows exponentially with the prediction step, so the exact reachable sets can be computed, and plotted, only for the first prediction steps; beyond them the computation becomes strictly intractable. The proposed interval-hull pipeline in `compute_set_evolution.m` exhibits instead a linear growth.


---


### Running the examples
```matlab
addpath(genpath('util'));

% Illustrative example: Figures 1 and 2
run('illustrative_example/main_illustrative.m');

% Mass-spring-damper example: Figures 3-5 and Table II
run('mass_spring_damper/main_ZRCG.m');
```

## Repository Structure

```
README.md
license.txt
util/
    RobustZRCommandGovernor.m      % ZR-CG class: online QP (Lemma 2 propagation,
                                   % epigraph reformulation, terminal constraint)
    compute_set_evolution.m        % forward propagation of the reachable sets
                                   % (interval hull + Lemma 2 enclosure)
    exact_paper_enclosure.m        % exact matrix-zonotope/zonotope product
                                   % (standard zonotopic outer-approximation)
    analyze_conservatism.m         % support-function gaps and volume ratios
                                   % between the exact and the proposed sets
illustrative_example/
    main_illustrative.m            % second-order example: identification of M_cl,
                                   % interval hull, projections (Figure 1) and
                                   % reachable set evolution (Figure 2)
mass_spring_damper/
    main_ZRCG.m                    % full supervision example: data collection,
                                   % identification, stability certificate,
                                   % terminal margins, closed-loop simulation
                                   % (Figures 3-5, Table II)
```

---

## Citation
If you use this software, please cite the associated paper:

```
@article{ZRCG2026,
  title   = {Computationally Tractable Robust Data-Driven Command Governors via Zonotopic Reachability},
  author  = {El Qemmah, Ayman and Casavola, Alessandro and Tedesco, Francesco},
  journal = {IEEE Transactions on Automatic Control},
  year    = {2026},
  doi     = {DOI TBA}
}
```

---

## License
This software is licensed under the custom LaSa, DIMES license (see headers in source files).

---
