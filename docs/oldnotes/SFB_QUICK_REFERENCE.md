# SFB Turbulence Quick Reference

## Essential Parameters

```
<turb_driving>
basis_type = 1      # REQUIRED: 1 for SFB, 0 for Cartesian
lmax       = 15     # Max angular momentum (10-20 typical)
nmax       = 10     # Max radial modes (5-10 typical)
r0_turb    = 0.4    # Turbulence radius (code units)
```

## Mode Selection
- Total modes before filtering: ~(lmax+1)² × nmax
- After k-filtering: typically 50-500 modes
- Memory: ~5 × N_modes × N_cells_active

## Power Spectrum Options

### Kolmogorov (default)
```
spect_form = 2      # Power law
expo       = 1.667  # 5/3
```

### Burgers
```
spect_form = 2      # Power law
expo       = 2.0    # Shock-dominated
```

### Peaked spectrum
```
spect_form = 1      # Parabolic
kpeak      = 25.13  # Peak at k = 4×2π
```

## Wavenumber Range
```
nlow  = 2           # k_min ~ 2×(2π/L)
nhigh = 8           # k_max ~ 8×(2π/L)
```
Note: Actual k-range for SFB modes:
- k_ln = x_ln / r0_turb
- x_ln = n-th root of spherical Bessel function j_l

## Time Evolution
```
dedt           = 0.1   # Energy injection rate
tcorr          = 0.5   # Correlation time (O-U process)
dt_turb_update = 0.01  # Update interval
```

## Typical Configurations

### Small-scale turbulence (star formation)
```
lmax = 20, nmax = 10, r0_turb = 0.1
nlow = 10, nhigh = 40
```

### Large-scale stirring (galaxy clusters)
```
lmax = 10, nmax = 5, r0_turb = 0.8
nlow = 1, nhigh = 5
```

### High-resolution study
```
lmax = 30, nmax = 15, r0_turb = 0.3
Grid: 512³ or larger
```

## Performance Tips

1. **Mode count**: N_modes ∝ lmax² × nmax × (k_range_fraction)
2. **Memory**: ∝ N_modes × (r0_turb/L)³ × N³
3. **CPU time**: Mode setup is one-time cost
4. **Parallel scaling**: Excellent (no communication)

## Common Issues

| Problem | Solution |
|---------|----------|
| Too few modes | Increase lmax or nmax |
| Poor k-coverage | Adjust nlow/nhigh |
| High memory use | Reduce lmax or increase r0_turb |
| Divergence at r=0 | Fixed in implementation |

## Analysis Commands

```bash
# Single file
python3 analyze_sfb_output.py output.athdf --r0_turb 0.4

# Time series
python3 analyze_sfb_output.py "*.athdf" --history file.hst

# Power spectrum
python3 analyze_sfb_output.py output.athdf --power-spectrum
```