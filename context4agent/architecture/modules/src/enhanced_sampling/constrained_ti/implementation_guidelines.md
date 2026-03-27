# constrained_ti — Implementation Guidelines

## Role

Constrained TI convergence diagnostics, free energy integration, and constant-potential correction. Diagnoses sampling sufficiency via 4-step pipeline, integrates dA/dxi with error propagation, and optionally applies Norskov correction for constant-charge -> constant-potential conversion.

## File Organization

| File | Purpose |
|------|---------|
| `config.py` | Threshold constants + output filename constants |
| `models.py` | Frozen dataclasses (all result types) + exception hierarchy |
| `workflow.py` | Orchestrator: analyze_single_point, analyze_standalone, analyze_ti, CSV exports |
| `integration.py` | Trapezoid weights, SEM targets, free energy integration, optimal allocation |
| `io.py` | Directory discovery (ti_target_*/xi_*, reverse support) + batch parsing |
| `plot.py` | 2x2 diagnostics PNG + free energy profile PNG |
| `correction.py` | Norskov constant-potential correction + corrected CSV/PNG |
| `analysis/` | Four independent diagnostic engines (see analysis/ mirror) |

## 4-Step Diagnostic Pipeline

1. **Running Average** — cumulative mean drift D < 3.0 * SEM
2. **ACF** — Sokal (1997) self-consistent cutoff -> tau_corr -> N_eff >= 50
3. **Block Average** — F&P (1989) pow2 blocks -> delta_SEM plateau detection
4. **Geweke** — front 10% vs rear 50% z-test (|z| < 1.96)

SEM selection: F&P plateau (primary) -> ACF fallback. Cross-check warning at >15% disagreement.

## Integration

- Non-uniform trapezoid: `compute_trapezoid_weights(xi)` supports ascending/descending xi
- Forces: dA/dxi = -lambda_mean (negated in workflow)
- Error: sigma_A = sqrt(sum(w^2 * sem^2))
- Units: lambda in a.u., only final delta_A converted to eV

## Constant-Potential Correction

`correction_eV = delta_sigma[e/A^2] * delta_phi[V] * A[A^2] / 2`

- sigma from per-point Bader ensemble averages (trajectory_surface_charge)
- phi from calibration mapper (mapper.predict)
- No error propagation on correction term
- Missing bader/ -> WARN + skip correction

## Module-Internal Constraints

- `workflow` -> `plot` (one-way); `plot` must NOT import `workflow`
- `io` must NOT import `integration`
- Analysis engines must NOT cross-import
- No reverse dependency: other packages do not depend on constrained_ti
