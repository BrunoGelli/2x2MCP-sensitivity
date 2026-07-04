# MCP rough sensitivity study

This repository contains a rough expected-yield sensitivity pipeline for millicharged particles (MCPs) crossing a liquid-argon pixel detector. The main notebook is:

```text
mcp_rough_sensitivity_v10_powerlaw_hybrid_commit_ready.ipynb
```

The current goal is not to produce a final statistical exclusion. The goal is to build a transparent, inspectable sensitivity estimate that separates production/geometry from detector response, and to understand how MCP charge, deposited-energy thresholds, and exposure affect the expected double-blip signal yield.

## Five-point summary

1. **Production and geometry are factorized into an acceptance table.**  
   The notebook starts from a table `A(m)` giving the expected number of MCP crossings through the detector per `epsilon^2` for the chosen exposure and geometry.

2. **Detector response is treated separately from production.**  
   For each crossing MCP track, the notebook estimates the probability of producing at least two above-threshold pixel blips.

3. **The deposited-energy histogram is used as a calibrated shape model, not as an absolute efficiency.**  
   The raw histogram normalization is ambiguous because zeros/underflow or missing pixel opportunities may not be fully represented. Therefore the histogram tail is normalized to the trusted reference point `p_pix(epsilon = 0.01, threshold = 1/3 MIP) = 2.37e-5`.

4. **The final histogram-tail response uses empirical survival where possible and a power-law extrapolation only where needed.**  
   A simple exponential tail was too punitive and produced unphysical-looking flat sensitivity contours. The power-law tail gives behavior that is smoother and more consistent with the q2 baseline and ArgoNeuT-like scaling.

5. **The final sensitivity curves are expected-yield contours.**  
   The notebook plots where `N_signal >= 3` under a zero-background rough-reach assumption. This is useful for model comparison and pipeline development, but it is not yet a full confidence-level limit.

## Physics logic

### 1. Production plus geometry

The calculation begins with a mass-dependent crossing normalization:

```text
N_cross(m, epsilon) = A(m) * epsilon^2
```

Here:

- `m` is the MCP mass;
- `epsilon` is the MCP electric charge in units of the electron charge;
- `A(m)` is the expected number of MCPs crossing the detector per `epsilon^2` for the chosen exposure and geometry.

The current notebook embeds the working acceptance table directly. This table should be treated as an input from the production/geometric-acceptance simulation, not something derived inside the sensitivity notebook.

Because `A(m)` falls over many orders of magnitude, the notebook interpolates it in log-log space using a shape-preserving PCHIP interpolation. This is important: ordinary cubic splines can create artificial bumps or fake improvements near steep acceptance drops. The notebook includes a dedicated plot of the acceptance table and its interpolation so this can be checked visually before looking at the final sensitivity curves.

### 2. MCP ionization and threshold scaling

For ionization-dominated energy deposition, the mean deposited energy scales approximately as:

```text
E_dep(epsilon) proportional to epsilon^2
```

Instead of rerunning the detector response at every charge, the notebook starts from a deposited-energy histogram at a reference charge:

```text
epsilon_ref = 0.01
```

A detector threshold `T` at charge `epsilon` is mapped into an equivalent threshold at the reference charge:

```text
E_equiv = T * (epsilon_ref / epsilon)^2
```

Then the probability of a pixel crossing threshold is controlled by the survival probability of the reference histogram above `E_equiv`.

### 3. Why the histogram is calibrated

The most tempting approach would be:

```text
p_pix = tail_counts_above_threshold / total_pixel_opportunities
```

However, this is dangerous unless the histogram is guaranteed to include every traversed pixel opportunity, including zero-deposit pixels and underflow. The current deposited-energy histogram is excellent for studying the shape of the tail, but its absolute denominator is not trusted enough to define the final detection efficiency.

Therefore the notebook uses the histogram as a shape model and calibrates it to the trusted q2 baseline:

```text
p_pix(epsilon_ref = 0.01, threshold = 1/3 MIP) = 2.37e-5
```

The calibrated histogram-tail model is:

```text
p_pix(epsilon, T) = p_pix_ref * S(E_equiv) / S(T_ref)
```

where:

- `S(E)` is the survival probability of the reference deposited-energy distribution above energy `E`;
- `T_ref` is the 1/3 MIP threshold at the reference charge;
- `p_pix_ref = 2.37e-5` is the trusted per-pixel reference probability.

This guarantees that the calibrated histogram model and the q2 baseline agree at the reference point.

### 4. Survival-tail model

The high-energy tail of the deposited-energy histogram controls the sensitivity at low charge. A simple exponential survival model,

```text
S(E) ~ exp(-constant * E)
```

was found to be too punitive, because after the charge rescaling it behaves roughly like:

```text
p_pix(epsilon) ~ exp(-constant / epsilon^2)
```

That produces flat, exposure-insensitive sensitivity curves that do not match the expected ArgoNeuT-like behavior.

The current recommended model uses the working v6/v9 behavior:

- use the empirical histogram survival where the histogram has enough support;
- fit the reliable region of the tail, currently about `0.25--0.60 MeV`;
- use a power-law extrapolation only outside the reliable histogram region;
- show exponential and stretched-exponential fits only as diagnostics.

The power-law extrapolation is not claimed to be final detector physics. It is a controlled, inspectable model that avoids the pathological behavior of a pure exponential tail while preserving the measured histogram shape where available.

### 5. Track-level double-blip probability

The signal topology is not a single above-threshold pixel. One MCP track crosses many pixels, and the analysis asks for at least two above-threshold blips along the track.

The current simplified model assumes a fixed average path length:

```text
n_pixels_crossed = 300
```

If the per-pixel above-threshold probability is `p_pix`, the probability of at least two above-threshold pixels along the track is:

```text
P_track_ge2 = 1 - (1 - p_pix)^N - N * p_pix * (1 - p_pix)^(N - 1)
```

where `N = n_pixels_crossed`.

This is the binomial probability for two or more successes in `N` trials. A future improvement should replace the fixed `N = 300` with a simulated path-length distribution.

### 6. Expected signal yield

The final expected signal yield is:

```text
N_signal(m, epsilon, T) = A(m) * epsilon^2 * P_track_ge2(epsilon, T)
```

The rough sensitivity contour is defined by:

```text
N_signal >= 3
```

This is a zero-background expected-yield reach criterion. It is useful for comparing threshold and exposure scenarios, but it is not yet a statistically rigorous exclusion limit.

## Programming logic

The notebook is organized as a linear pipeline.

### 1. Global configuration

The first code section defines:

- `epsilon_ref = 0.01`;
- `min_signal_events = 3.0`;
- `n_pixels_crossed = 300`;
- `mip_energy_mev = 0.9`;
- `threshold_ref_mip_fraction = 1/3`;
- `p_pix_ref_at_threshold_ref = 2.37e-5`;
- tail-fit range and plot ranges.

These constants are intentionally centralized so the assumptions are easy to audit.

### 2. Acceptance table and interpolation

The notebook embeds two arrays:

```python
mass_data
acc_over_eps2_base
```

These define `A(m)`. The function `build_loglog_interpolator` creates the mass interpolation used by the scan. If SciPy is available, the notebook uses `PchipInterpolator`; otherwise it falls back to log-log linear interpolation.

The acceptance diagnostic plot is one of the most important sanity checks. If the spline is not smooth and monotonic-looking where expected, the final sensitivity contours should not be trusted.

### 3. Histogram loading

The notebook looks for deposited-energy histograms with these names:

```text
mcp_edep_histogram_eps0p01.txt
mcp_edep_histogram_eps0p01 (2).txt
```

The expected file format is two columns:

```text
edep_mev   counts
```

Alternatively, arrays can be pasted directly into:

```python
edep_hist_mev
edep_hist_counts
```

The raw histogram denominator is used only for diagnostics. It is not used to set the final absolute detection efficiency.

### 4. Survival model

The notebook builds several survival functions from the reference histogram:

- empirical survival;
- exponential-in-energy fit;
- power-law fit;
- stretched-exponential fit, if SciPy is available.

Only the hybrid empirical-plus-power-law model is used for the final sensitivity curves. The other fits are kept in the survival-tail plot so the model choice remains visible.

### 5. Detector-response functions

The main response functions are:

```python
p_pix_q2_reference(eps)
p_pix_histogram_tail_calibrated_power_law(eps, threshold_mip_fraction)
p_track_ge2_from_p_pix(p_pix, n_pixels_crossed)
p_track_ge2_q2_reference(eps)
p_track_ge2_hist_cal_power_law(eps, threshold_mip_fraction)
```

The q2 reference is used as a check and comparison model. The calibrated histogram-tail function is the recommended detector-response model for the final scenarios.

### 6. Scenario definition

The final notebook uses exactly four scenarios:

```text
1/3 MIP, current exposure
1/3 MIP, 20x exposure
1/4 MIP, current exposure
1/4 MIP, 20x exposure
```

Each scenario changes only the threshold and/or exposure scale. The same production table and detector-response logic are used for all four.

### 7. Final sensitivity scan

For each mass and charge on the scan grid, the notebook evaluates:

```text
N_signal = exposure_scale * A(m) * epsilon^2 * P_track_ge2(epsilon, threshold)
```

Then it extracts the contour where `N_signal = 3`. External comparison curves, such as the ArgoNeuT two-blip limit, are loaded and plotted as overlays.

## Validation checks included in the notebook

The notebook includes several checks that should be run before interpreting the final plot.

### Acceptance check

The acceptance table and interpolation are plotted together. This checks that the production/geometric input is smooth and that the interpolator is not introducing artifacts.

### Reference-point q2 check

The calibrated histogram model must reproduce the trusted q2 baseline at the reference point:

```text
epsilon = 0.01
threshold = 1/3 MIP
```

The notebook prints this explicitly. The expected result is:

```text
q2 p_pix(eps_ref)                     = 2.37e-5
hist-cal p_pix(eps_ref, 1/3 MIP)      = 2.37e-5
q2 P_track_ge2(eps_ref)               = hist-cal P_track_ge2(eps_ref, 1/3 MIP)
```

If this check fails, the sensitivity plot should not be used.

### Track-probability plot

The track-level double-blip probability is plotted as a function of `epsilon`. This plot shows whether the histogram-tail model behaves reasonably compared with the q2 baseline.

### Survival-tail plot

The survival-tail plot shows the empirical histogram survival and the three candidate fit shapes. This is the main diagnostic for understanding whether the extrapolation is too aggressive or too gentle.

## Inputs

Recommended files to commit:

```text
mcp_rough_sensitivity_v10_powerlaw_hybrid_commit_ready.ipynb
mcp_edep_histogram_eps0p01.txt
ArgoNeuT - 2 blip.txt
README.md
```

Optional supporting files:

```text
mcp_rough_sensitivity_v10_powerlaw_hybrid_commit_ready_executed.ipynb
```

The executed notebook is useful for quick review, but the clean notebook is the preferred source file for development.

## Outputs

The notebook produces:

1. production/geometric acceptance interpolation plot;
2. deposited-energy histogram plot;
3. survival-tail plot with fit comparison;
4. track-level double-blip probability plot;
5. final sensitivity plot with the four threshold/exposure scenarios and external overlays.

## Interpreting the final plot

Lower curves correspond to better sensitivity in charge. The `20x exposure` curves should improve relative to current exposure. The `1/4 MIP` curves should improve relative to `1/3 MIP`, but the improvement depends on the histogram-tail shape rather than a simple q2 scaling.

The q2 baseline and the calibrated histogram-tail model are forced to agree only at the reference point:

```text
epsilon = 0.01, threshold = 1/3 MIP
```

They are not expected to agree at all charges, because the histogram-tail model contains additional threshold-shape information.

## Known limitations

This is still a rough sensitivity model. The main limitations are:

1. **Zero-background criterion**  
   The contour uses `N_signal >= 3`, not a full statistical treatment with backgrounds and systematics.

2. **Fixed path length**  
   The model uses `n_pixels_crossed = 300`. A future version should average over the simulated path-length distribution.

3. **Histogram-shape dependence**  
   The detector response depends on the assumed survival-tail extrapolation. The power-law model is currently preferred because it behaves well and avoids the pathological exponential suppression, but it should be validated with additional toy MC charge points.

4. **Histogram normalization is not used absolutely**  
   This is intentional. The histogram may not represent every pixel opportunity, so the model is calibrated to the trusted q2 reference point.

5. **Production table provenance matters**  
   The embedded `A(m)` table must match the intended geometry, POT, meson production model, and mass range. If the production code changes, this table should be regenerated and the acceptance spline plot should be checked again.

6. **External overlays need unit checks**  
   The ArgoNeuT overlay file may label mass units differently from the numerical scale used in the plot. Always verify whether the external curve should be interpreted as GeV or MeV before comparing.

## Recommended next steps

Before treating this as more than a rough reach estimate, the most important improvements are:

1. rerun the detector toy MC at one or more additional charges, for example `epsilon = 0.005` and/or `epsilon = 0.02`;
2. directly compare the measured charge scaling of the per-pixel threshold probability against the calibrated histogram-tail prediction;
3. store per-track outputs: number of pixels crossed, number above threshold, and per-pixel deposited energies;
4. replace the fixed `n_pixels_crossed = 300` with the simulated path-length distribution;
5. add a background estimate and replace the `N_signal >= 3` criterion with a proper statistical limit.

## Current status

The current recommended notebook is `v10`. It should be considered a clean, commit-ready rough sensitivity pipeline with a transparent calibrated histogram-tail detector response. The physics architecture is now reasonable and inspectable, but the detector-response model still needs validation from additional toy MC samples before being used as a final result.
