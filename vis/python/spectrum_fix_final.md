For the **2D scalar-field path itself**, the core spectrum logic is now basically sound:

* centered shells via `floor(kmag + 0.5)` are a good choice;
* `shell_sum` is the exact discrete shell-integrated power;
* `2πk * shell_mean` is a reasonable continuum-style isotropic estimator;
* dividing rebinned sums by bin width is the right thing if you want to preserve slopes.

So I do **not** see another hidden low-(k) bug of the same kind.

What I do still see is this:

## 1. The continuum spectrum will be wrong near Nyquist unless `kmax` is conservative

This is the biggest remaining mathematical issue.

Your plotted spectrum is

```python
continuum = 2.0 * np.pi * kvals * shell_mean
```

in 2D. That assumes a **full circular annulus** exists at radius (k). But on an FFT grid, Fourier space is a square, not a disk. Near Nyquist, the annulus is clipped by the square boundary.

So:

* `shell_sum_raw` is still valid all the way out;
* `spec_raw` / `spec_continuum_raw` are **not** valid arbitrarily close to the grid cutoff.

This is only okay if `base._spectrum_kmax_from_field(...)` already limits you to the inscribed-circle range.

For 2D, the safe continuum cutoff is roughly

```python
safe_kmax_2d = int(np.floor(min(nx, ny) / 2.0 - 0.5))
```

For 3D full-shell continuum spectra, similarly:

```python
safe_kmax_3d = int(np.floor(min(nx, ny, nz) / 2.0 - 0.5))
```

For your 3D slice-averaged 2D spectra, the same idea applies per slice plane; the common safe cutoff is again tied to the smallest grid dimension if spacing is uniform.

If `base._spectrum_kmax_from_field` is not doing this, then the high-(k) end of `spec_raw` and `spec_rebinned` is biased high.

## 2. The 3D slice-average path is not a true 3D isotropic spectrum

This is a conceptual issue, not a coding accident.

In 3D you are doing 2D FFTs on slices and averaging them:

```python
plane = np.take(array, indices=index, axis=axis)
kvals, shell_sum, shell_count, _, _ = _theta_spectrum(plane)
```

That gives a **2D slice spectrum**, not the full 3D shell spectrum of the 3D field.

So if someone later interprets the 3D outputs as “the 3D scalar spectrum,” that would be wrong in amplitude and definition. It can still be a useful diagnostic, and inertial-range slopes may still be informative, but it is not interchangeable with a full 3D FFT spectrum.

Your warning string helps, but I would make that distinction very explicit in filenames / metadata / plotting labels too.

## 3. `_average_raw_shell_stats()` throws away valid data on non-cubic grids

This is a real structural bug if the 3D domain is not cubic.

You do:

```python
common_length = min(len(shell_sum) for _, shell_sum, _ in spectra)
```

and then truncate every slice spectrum to that minimum length.

That means:

* larger planes lose their valid higher-(k) content;
* the combined spectrum is forcibly clipped to the smallest slice orientation.

If all your production grids are cubic, this is harmless. If not, it is wrong.

A better combination is to accumulate to the **maximum** length and let `shell_count` be zero where a given slice has no support:

```python
def _combine_raw_shell_stats(
    spectra: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    max_length = max(len(shell_sum) for _, shell_sum, _ in spectra)
    raw_sum = np.zeros(max_length, dtype=np.float64)
    raw_count = np.zeros(max_length, dtype=np.float64)

    for _, shell_sum, shell_count in spectra:
        m = len(shell_sum)
        raw_sum[:m] += np.asarray(shell_sum, dtype=np.float64)
        raw_count[:m] += np.asarray(shell_count, dtype=np.float64)

    raw_k = np.arange(1, max_length + 1, dtype=np.float64)
    return raw_k, raw_sum, raw_count
```

Then apply a physically justified `kmax` afterward.

One nuance: that combines all slice modes together. If you instead want **equal weight per slice**, you should average per-slice continuum spectra, not raw sums/counts.

## 4. Your `k` is still in grid-mode units, not physical wavenumber

You define

```python
kx = np.fft.fftfreq(nx) * nx
```

so `k` is an integer Fourier mode number.

That is fine for slopes. But if you want physical wavenumber, you still need the box length:

[
k_\mathrm{phys} = \frac{2\pi}{L} , n
]

So if later you compare to forcing scales, Batchelor scales, Kolmogorov scales, or any dimensional theory, this code still needs a domain-length conversion.

## 5. The velocity “spectrum” is power, not kinetic energy, unless that is intentional

Here you use

```python
energy = (|uhat_x|^2 + |uhat_y|^2 + |uhat_z|^2) / N^2
```

If you want a **kinetic-energy spectrum** in the conventional sense, there is usually a factor of `0.5`.

So either:

* rename it conceptually as a velocity-power spectrum, or
* multiply by `0.5`.

This does not affect the slope, only the normalization.

## 6. Naming / metadata are a bit misleading

A few fields will confuse future-you:

* `theta_spectrum_raw` is **not** the raw discrete shell sum. It stores the continuum estimator.
* the true raw discrete spectrum is `theta_shell_sum_raw`.
* `theta_spectrum_raw` and `theta_spectrum_continuum_raw` are duplicates.
* the file docstring says `"shell-integrated continuum spectra"`, which mixes two different ideas.

I would rename for clarity:

* `theta_shell_sum_raw`
* `theta_shell_mean_raw`
* `theta_spectrum_continuum_raw`

and drop or rename `theta_spectrum_raw`.

Also:

```python
method=np.asarray(str(meta.get("method", "")))
```

looks unrelated to the actual spectral method. The reliable field is `spectrum_method`, not `method`.

## 7. Minor: the 3D continuum factor is only asymptotically correct

If you ever use the full 3D branch again, you have

```python
4.0 * np.pi * k**2 * shell_mean
```

For half-integer shell bins, the exact shell volume is

[
\frac{4\pi}{3}\left[(k+\tfrac12)^3 - (k-\tfrac12)^3\right]
= 4\pi k^2 + \frac{\pi}{3}.
]

So `4πk²` is the large-(k) approximation, not the exact low-(k) shell measure.

This is minor, and your production path seems to avoid full 3D anyway, but it is a real detail.

## 8. The first few bins can still be jagged, and that is expected

Even after this fix, the lowest few (k) bins are still based on very few modes. So some wiggle there is normal.

That is no longer a sign that the shelling is broken.
