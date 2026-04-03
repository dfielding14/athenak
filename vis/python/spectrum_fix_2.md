Your hunch is good, but the `ceil` is only part of the story.

The main issue is this:

`np.ceil(np.sqrt(kx2*kx2 + ky2*ky2) - 1e-12)`

does **hard shell assignment** with shells

* shell 1: (0 < |k| \le 1)
* shell 2: (1 < |k| \le 2)
* shell 3: (2 < |k| \le 3)

So in 2D, your first shell contains only the four axis modes:
((\pm1,0)), ((0,\pm1)).

That is extremely sparse. A 4096×4096 grid still has only 4 modes in that shell. So low-(k) jaggedness is expected even on huge grids.

## What `ceil` is doing

It is not “wrong” mathematically, but it is an awkward convention for an isotropic spectrum because it labels shells by their **upper edge**, not their center.

For example:

* (|k|=1) goes to shell 1
* (|k|=\sqrt{2}\approx1.41) goes to shell 2
* (|k|=2) goes to shell 2

So your `k=1` shell is unusually under-populated, which makes the first point especially noisy.

A more standard choice is to use shells centered on integer radii:

[
k-\tfrac12 \le |k| < k+\tfrac12
]

That corresponds to:

```python
shell = np.floor(kmag + 0.5).astype(np.int64)
```

or equivalently `np.rint(kmag).astype(np.int64)`.

That change alone usually makes the lowest few shells look less pathological.

## A better shell definition

For 2D:

```python
def _spectrum_shells_2d(nx: int, ny: int) -> np.ndarray:
    kx = np.fft.fftfreq(nx) * nx
    ky = np.fft.fftfreq(ny) * ny
    ky2, kx2 = np.meshgrid(ky, kx, indexing="ij")
    kmag = np.hypot(kx2, ky2)
    return np.floor(kmag + 0.5).astype(np.int64)
```

For 3D:

```python
def _spectrum_shells_3d(nx: int, ny: int, nz: int) -> np.ndarray:
    kx = np.fft.fftfreq(nx) * nx
    ky = np.fft.fftfreq(ny) * ny
    kz = np.fft.fftfreq(nz) * nz
    kz3, ky3, kx3 = np.meshgrid(kz, ky, kx, indexing="ij")
    kmag = np.sqrt(kx3 * kx3 + ky3 * ky3 + kz3 * kz3)
    return np.floor(kmag + 0.5).astype(np.int64)
```

## What is really causing the jagged low-k spectrum

Even with better shelling, low-(k) in 2D is often jagged because the number of modes per annulus is small.

With your current `ceil` binning, the first few shell occupancies are roughly:

* (k=1): 4 modes
* (k=2): 8 modes
* (k=3): 16 modes
* (k=4): 20 modes

That is not enough for a smooth radial average.

So the honest answer is:

* **yes**, the `ceil` choice is contributing;
* **no**, replacing it will not completely remove the jaggedness.

## Two spectra are getting mixed together

Right now you compute:

```python
spectrum = np.bincount(shell.ravel(), weights=power.ravel())
```

That is the **total power in each shell**. This is fine if you want a shell-integrated spectrum whose sum gives the variance.

But shell sums in 2D inherit the irregular integer-lattice mode counts. That alone makes them jagged.

A very useful diagnostic is to also compute the number of modes per shell:

```python
counts = np.bincount(shell.ravel())
```

Then compare:

```python
shell_sum = spectrum[1:]
shell_mean = shell_sum / counts[1:]
```

* `shell_sum` = total power in each annulus
* `shell_mean` = average power per Fourier mode in that annulus

If `shell_mean` is much smoother, then the jaggedness is mostly shell occupancy, not physics.

For a smoother 2D isotropic estimate, a common choice is:

```python
E2d = 2.0 * np.pi * kvals * shell_mean
```

That uses the smooth continuum ring area (2\pi k) instead of the exact discrete mode count. It will usually look much less jagged at small (k).

## A practical rewrite

This keeps your current normalization and adds the shell counts:

```python
def _theta_spectrum(field: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    array = _canonical_field(field)

    theta_mean = float(np.mean(array))
    fluctuating = array - theta_mean

    shat = np.fft.fftn(fluctuating)
    power = np.abs(shat) ** 2 / array.size**2

    if array.ndim == 2:
        shell = _spectrum_shells_2d(array.shape[1], array.shape[0])
    else:
        shell = _spectrum_shells_3d(array.shape[2], array.shape[1], array.shape[0])

    shell_sum = np.bincount(shell.ravel(), weights=power.ravel())
    shell_count = np.bincount(shell.ravel())

    kvals = np.arange(1, len(shell_sum), dtype=np.float64)
    return kvals, shell_sum[1:], shell_count[1:], theta_mean
```

Then:

```python
kvals, shell_sum, shell_count, theta_mean = _theta_spectrum(field)

shell_mean = shell_sum / np.maximum(shell_count, 1)

# shell-integrated spectrum: sums to variance
E_shell = shell_sum

# smoother 2D continuum-style estimate
E_2d = 2.0 * np.pi * kvals * shell_mean
```

## One more thing: your log rebinning

There is a separate issue in `_log_rebin_spectrum`:

```python
rebinned.append(np.sum(spectrum[mask]))
```

This gives the **total power in each logarithmic bin**, not a spectral density. Since the bins get wider with (k), the plotted slope changes.

For a rebinned spectrum density, divide by the bin width:

```python
rebinned.append(np.sum(spectrum[mask]) / (hi - lo))
```

That will not fix the low-(k) jaggedness, but it will make the rebinned spectrum more interpretable.

## Best ways to reduce the jaggedness

In order of usefulness:

1. Change from `ceil(...)` to centered shells via `floor(kmag + 0.5)`.
2. Track `shell_count` and inspect `shell_sum / shell_count`.
3. For a smoother 2D isotropic estimate, use `2πk * shell_mean`.
4. Rebin with broader low-(k) bins, or bins containing a minimum number of modes.
5. Average spectra over multiple snapshots if you have them.

One subtle but important point: making the grid larger does **not** give you more modes in the very first shells. The `k=1` and `k=2` shells stay sparse no matter how large the domain is.

## The `-1e-12` is not the real problem

That tiny subtraction is just there to keep exact integers from rounding up accidentally. The issue is the use of `ceil`, not the epsilon.

The shortest diagnosis is:

* `ceil` makes your lowest shells badly centered and under-populated.
* low-(k) jaggedness in 2D is mostly a finite-mode-count effect.
* use centered shells and consider plotting `2πk ×` the shell-averaged power.

Paste a small example spectrum or a plot of `shell_count` versus `k`, and I can tell you very quickly whether the roughness is mostly binning noise or something physical.
