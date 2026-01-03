# Notebook Execution Time Analysis

This document contains the execution time analysis for all notebooks based on a successful GitHub Actions run (run ID: 20676974568).

## Summary

- **Total notebooks processed**: 51
- **Total execution time**: 4778.03 seconds (~79.6 minutes)
- **Average execution time**: 93.69 seconds per notebook
- **Threshold for long-running notebooks**: 150 seconds (2.5 minutes)

## Long-Running Notebooks (Moved to `notebooks/long_running/`)

These notebooks take more than 150 seconds to execute and have been moved to a separate folder for optimization:

| Notebook | Time (s) | Time (min) |
|----------|----------|------------|
| N-002-fit-benchmark.jl | 359.18 | 5.99 |
| N-008-regge_NN_pluto.jl | 294.78 | 4.91 |
| N-093-Xtoppbar.jl | 289.87 | 4.83 |
| N-078-CL-KDE.jl | 280.05 | 4.67 |
| N-065-LcXX2Lcpipi.jl | 211.41 | 3.52 |
| N-094-cut-on-dalitz.jl | 195.22 | 3.25 |
| N-077-mixed-event-sample.jl | 181.11 | 3.02 |
| N-048-resonancesondalitzplot.jl | 165.17 | 2.75 |

**Total time for long-running notebooks**: 1976.79 seconds (~32.9 minutes) - 41.4% of total time

## All Notebooks (Sorted by Execution Time)

| Notebook | Time (s) | Time (min) | Status |
|----------|----------|------------|--------|
| N-002-fit-benchmark.jl | 359.18 | 5.99 | ⚠️ Long-running |
| N-008-regge_NN_pluto.jl | 294.78 | 4.91 | ⚠️ Long-running |
| N-093-Xtoppbar.jl | 289.87 | 4.83 | ⚠️ Long-running |
| N-078-CL-KDE.jl | 280.05 | 4.67 | ⚠️ Long-running |
| N-065-LcXX2Lcpipi.jl | 211.41 | 3.52 | ⚠️ Long-running |
| N-094-cut-on-dalitz.jl | 195.22 | 3.25 | ⚠️ Long-running |
| N-077-mixed-event-sample.jl | 181.11 | 3.02 | ⚠️ Long-running |
| N-048-resonancesondalitzplot.jl | 165.17 | 2.75 | ⚠️ Long-running |
| N-056-Chew-Mandelstam.jl | 139.74 | 2.33 | ✅ Normal |
| N-062-poisson-CL.jl | 133.46 | 2.23 | ✅ Normal |
| N-016-XicK_phasespace_suppression.jl | 125.33 | 2.09 | ✅ Normal |
| N-051-Xibpipi-feeddown.jl | 124.55 | 2.08 | ✅ Normal |
| N-059-PDG-API.jl | 122.36 | 2.04 | ✅ Normal |
| N-049-singlechannelkmatrix.jl | 116.07 | 1.93 | ✅ Normal |
| N-017-UnitarizedBackground.jl | 101.37 | 1.69 | ✅ Normal |
| N-066-virtual-and-bound-states.jl | 101.21 | 1.69 | ✅ Normal |
| N-060-double-regge.jl | 94.88 | 1.58 | ✅ Normal |
| N-063-poisson-profile.jl | 94.07 | 1.57 | ✅ Normal |
| N-054-B2jpsietaK.jl | 90.92 | 1.52 | ✅ Normal |
| N-054-Mmatrix-and-Kmatrix.jl | 89.51 | 1.49 | ✅ Normal |
| N-061-complex-random-variance.jl | 88.62 | 1.48 | ✅ Normal |
| N-039-logsqrt_complex_structure.jl | 88.00 | 1.47 | ✅ Normal |
| N-038-Spuriouspole_a1.jl | 86.44 | 1.44 | ✅ Normal |
| N-057-a12piKK_toys.jl | 75.08 | 1.25 | ✅ Normal |
| N-001-Adams_regge_integral.jl | 74.68 | 1.24 | ✅ Normal |
| N-012-Pythia_FHist_plots.jl | 73.20 | 1.22 | ✅ Normal |
| N-046-LcLcK_momentumdistribution.jl | 69.15 | 1.15 | ✅ Normal |
| N-003-advanced_bw.jl | 69.14 | 1.15 | ✅ Normal |
| N-029-introduction2complexplane.jl | 55.12 | 0.92 | ✅ Normal |
| N-013-MixedModelExample.jl | 54.58 | 0.91 | ✅ Normal |
| N-035-numconv_prototype.jl | 53.72 | 0.90 | ✅ Normal |
| L-014-MixedModel-chi2.jl | 50.77 | 0.85 | ✅ Normal |
| N-041-derivativeBWonArgand.jl | 48.12 | 0.80 | ✅ Normal |
| L-052-Xib2Xibpipi.jl | 44.27 | 0.74 | ✅ Normal |
| N-018-LcLcK.jl | 42.58 | 0.71 | ✅ Normal |
| N-050-pipiDeck.jl | 42.58 | 0.71 | ✅ Normal |
| N-042-detDsymbolics.jl | 42.43 | 0.71 | ✅ Normal |
| L-044-DsDsState.jl | 42.27 | 0.70 | ✅ Normal |
| L-047-effectiverangeof2bloop.jl | 37.49 | 0.62 | ✅ Normal |
| N-058-dispesion-blatt-weisskopf.jl | 34.13 | 0.57 | ✅ Normal |
| L-055-Kmatrix-sliders.jl | 33.67 | 0.56 | ✅ Normal |
| N-022-OcXX2XicK.jl | 31.27 | 0.52 | ✅ Normal |
| N-045-a1dalitz.jl | 29.47 | 0.49 | ✅ Normal |
| N-028-pole_search_optim.jl | 29.08 | 0.48 | ✅ Normal |
| N-026-fourier_massspectrum.jl | 28.96 | 0.48 | ✅ Normal |
| L-024-X2Dx0D0Triangle.jl | 27.31 | 0.46 | ✅ Normal |
| L-025-psi-triangle.jl | 27.30 | 0.46 | ✅ Normal |
| L-032-quiverplot.jl | 27.16 | 0.45 | ✅ Normal |
| N-033-LcXX2Lcpipi_angle_crosscheck.jl | 22.12 | 0.37 | ✅ Normal |
| N-070-TrackLorentzGroup.jl | 21.95 | 0.37 | ✅ Normal |
| N-030-merging_poles.jl | 17.08 | 0.28 | ✅ Normal |

## Recommendations

1. **Optimize long-running notebooks**: The 8 notebooks in `notebooks/long_running/` account for 41.4% of total execution time. Consider:
   - Reducing computation complexity
   - Using cached results where possible
   - Breaking down into smaller notebooks
   - Using approximate methods for expensive calculations

2. **CI/CD Strategy**: 
   - Consider running long-running notebooks separately or on a schedule
   - Use parallel execution for independent notebooks
   - Cache results between runs when possible

3. **Monitoring**: Re-run the timing analysis periodically to track improvements.

## How to Re-run Analysis

```bash
# Get a successful run ID
gh run list --workflow=pages.yaml --limit 10

# Extract timing data
python3 scripts/extract_notebook_times.py <run_id>
```

