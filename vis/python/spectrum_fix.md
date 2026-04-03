We’ve stumbled into a classic problem in discrete Fourier spectral analysis often related to the **Gauss circle problem**.

Here is what is happening and how to fix it.

### The Problem: Mode Sparsity and the "Gauss Circle Problem"

When you take an FFT on a Cartesian grid, the available wavevectors $(k_x, k_y)$ exist only at discrete integer coordinates. To get a 1D power spectrum $E(k)$, you are grouping these discrete points into circular (or spherical) shells based on their radius $k = \sqrt{k_x^2 + k_y^2}$.

Currently, your code calculates the spectrum by **summing** the power of all modes that fall into a given integer shell using `np.bincount(..., weights=power)`.

The issue is that at small $k$, the number of discrete grid points that fall into an annulus fluctuates wildly from one bin to the next. For instance:

* A shell at $k=1$ might contain exactly 4 modes: $(\pm1, 0)$ and $(0, \pm1)$.
* A shell at $k=2$ might contain 4 modes.
* A shell at $k=3$ might contain 16 modes.

Because the mode count jumps up and down randomly at small $k$, your **summed power will be jagged**, perfectly mirroring the fluctuating mode counts rather than the actual physical spectrum of your field.

### The Solution: Average and Rescale

To get a smooth spectrum, you must detach your power calculation from the discrete grid's geometry. You can do this in two steps:

1. **Calculate the *mean* power per mode** in each shell (by dividing your summed power by the exact number of modes in that shell). This gives you an accurate, smooth measure of the 2D Power Spectral Density, $P(k)$.
2. **Multiply by the theoretical continuous area/volume** of the shell. In 2D, the shell circumference is $2\pi k$. In 3D, the spherical surface area is $4\pi k^2$.

By doing $E(k) = \text{mean}(|\hat{\theta}|^2) \times 2\pi k$, you completely bypass the discrete mode counting noise.

### The Code Fix

Here is how you can update your code. I also swapped your `np.ceil(k - 1e-12)` logic for `np.round(k)` inside the shell generators. Rounding ensures that your bins (e.g., bin 1 covers $k \in [0.5, 1.5)$) are perfectly centered on the integer `kvals` you output.

```python
import numpy as np

# A small constant to avoid magic numbers
SPECTRUM_NBINS = 100

def _spectrum_shells_2d(nx: int, ny: int) -> np.ndarray:
    kx = np.fft.fftfreq(nx) * nx
    ky = np.fft.fftfreq(ny) * ny
    ky2, kx2 = np.meshgrid(ky, kx, indexing="ij")
    # Using np.round centers the bins perfectly on integer k values
    return np.round(np.sqrt(kx2 * kx2 + ky2 * ky2)).astype(np.int64)

def _spectrum_shells_3d(nx: int, ny: int, nz: int) -> np.ndarray:
    kx = np.fft.fftfreq(nx) * nx
    ky = np.fft.fftfreq(ny) * ny
    kz = np.fft.fftfreq(nz) * nz
    kz3, ky3, kx3 = np.meshgrid(kz, ky, kx, indexing="ij")
    return np.round(np.sqrt(kx3 * kx3 + ky3 * ky3 + kz3 * kz3)).astype(np.int64)

def _theta_spectrum(field: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    array = _canonical_field(field)
    theta_mean = float(np.mean(array))
    fluctuating = array - theta_mean
    shat = np.fft.fftn(fluctuating)
    power = np.abs(shat) ** 2 / array.size**2

    if array.ndim == 2:
        shell = _spectrum_shells_2d(array.shape[1], array.shape[0])
        area_factor = 2.0 * np.pi
        exponent = 1.0
    else:
        shell = _spectrum_shells_3d(array.shape[2], array.shape[1], array.shape[0])
        area_factor = 4.0 * np.pi
        exponent = 2.0

    # 1. Sum the power in each discrete shell
    spectrum_sum = np.bincount(shell.ravel(), weights=power.ravel())

    # 2. Count the exact number of modes that fell into each shell
    mode_counts = np.bincount(shell.ravel())

    # 3. Calculate the mean power per mode (handling empty bins gracefully)
    valid = mode_counts > 0
    mean_power = np.zeros_like(spectrum_sum)
    mean_power[valid] = spectrum_sum[valid] / mode_counts[valid]

    # kvals start at 1 to match spectrum[1:]
    kvals = np.arange(1, len(spectrum_sum), dtype=np.float64)

    # 4. Extrapolate to the continuous analytical shell area/volume
    spectrum = mean_power[1:] * area_factor * (kvals ** exponent)

    return kvals, spectrum, theta_mean

```
