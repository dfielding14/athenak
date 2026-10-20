# GOTHAM dust rates against their sources

Audit of `DustSource` in `src/pgen/gotham.cpp`, 2026-10-19. **No code was
changed.** Two of the four processes are normalized differently from the
papers the code cites:

| Process | Code vs. source | Effect in the code |
|---|---|---|
| Coagulation | prefactor 2.7 Myr; Dubois et al. (2024) Eq. (20) has 0.27 Myr | 10× too slow |
| Thermal sputtering | Tsai & Mathews (1995) rate at n_H = 1e-3 cm^-3; it belongs at n_H = 4.5e-4 cm^-3 | 2.2× too slow |
| Shattering | matches Dubois et al. (2024) Eq. (18) | — |
| Accretion | follows `docs/dust_growth_time_note.tex` | not audited here |

The legacy generator `src/pgen/cgm_cooling_flow_full.cpp` carries the same
expressions.

## How the code applies a timescale

Each process removes `req = dt * rho_D / (tau * f(a))` per step, so the
effective mass timescale is `tau * f(a)`: the grain-size factor sits next to
`tau` in the denominator. Shattering confirms this reading. The code has
`tau_shatt_Myr = 540` and divides by `a_large` in µm, giving 540 × 0.1 =
54 Myr at a = 0.1 µm, exactly Dubois et al. Eq. (18):

```latex
t_{\rm sha,i}= 54 a_{\rm 0.1,L}  s_{3,i} \left(\frac{\mathcal{D}_{\rm L}}{0.01}\right)^{-1}\left(\frac{n}{1\,\rm H\,cm^{-3}}\right)^{-1}\left(\frac{\sigma_{\rm D,L}}{10\,\rm km\,s^{-1}}\right)^{-1}\, \rm Myr\, ,
```

## Coagulation: 10× too slow

Code (`gotham.cpp`, section 4 of `DustSource`; `f_coag = 0.5`,
`v_coag_ref_kms = 0.1`):

```cpp
Real tau_coag_Myr = 2.7;
tau_coag_Myr *= (rho_grain_cgs / 3.0);
tau_coag_Myr *= (1.0e3 / nH);
tau_coag_Myr *= (0.01 / D_small);
tau_coag_Myr *= (0.1 / v_coag_ref_kms);
tau_coag_Myr /= f_coag;
...
const Real req = dt_Myr * rho_Ds / (tau_coag_Myr * (a_small / 0.005));
```

This gives an effective timescale of

    t_coa(code) = 2.7 Myr (a/0.005 µm) (s/3 g cm^-3) / F (D_S/0.01)^-1 (n_H/10^3 cm^-3)^-1 (v/0.1 km s^-1)^-1

Dubois et al. 2024, A&A 687, A240 (arXiv:2402.18515v2), Sect. 3.6, Eq. (20).
This is verbatim from the arXiv source, `article.tex` line 527:

```latex
t_{\rm coa}= 0.27 \frac{a_{\rm 0.005} s_3}{F} \left(\frac{\mathcal{D}_{\rm S}}{0.01}\right)^{-1}\left(\frac{n}{10^3\,\rm H\, cm^{-3}}\right)^{-1} \left( \frac{\sigma_{\rm D,S}}{0.1\, {\rm km\, s^{-1}}} \right)^{-1}\, \rm Myr\, ,
```

- Here `a_0.005 = a/5 nm` and `s_3 = s/3 g cm^-3` (line 444), and F is a fudge
  factor with fiducial value F = 0.5 (line 532).
- Aoyama et al. 2017, MNRAS 466, 105 (arXiv:1609.07547), Sect. 2.3.2, Eq. (22)
  gives the same number: `ken.tex` lines 881-892 have
  `2.71 \times 10^{5} yr (D_S/0.01)^{-1} (v_co/0.1 km s^-1)^{-1} / f_dense`,
  evaluated at a = 5e-3 µm and s = 3 g cm^-3, for dense clouds at
  n_H = 1e3 cm^-3.

Every factor in the code matches Eq. (20): `a/0.005`, `s/3`, `/F`, `D_S`,
`n_H` and `v`. Only the prefactor differs: 2.7 Myr instead of 0.27 Myr. With
GOTHAM's inputs (`a_gr_s = 0.003`, `rho_gr = 3`) at n_H = 1e3 cm^-3 and
D_S = 0.01:

- Eq. (20) gives 0.324 Myr;
- the code gives 3.24 Myr.

**How it got this way.**
- Before 5b4d5c5c7 (2026-04-14, "Fix bugs in dust module"),
  `cgm_cooling_flow_full.cpp` had `t_coag = 0.27 * F` with size factor
  `a_s / 0.05`. That is 0.0135 Myr × (a/0.005 µm), 40× faster than Eq. (20):
  it multiplied by F instead of dividing (4×), and normalized to 0.05 µm
  instead of 0.005 µm (10×).
- 5b4d5c5c7 corrected F but also changed both 10× factors in the same
  direction (prefactor 0.27 → 2.7, size factor 0.05 → 0.005), leaving the
  code 10× too slow.
- Its new comment, `2.7 Myr ... × f_coag` with `a_small / 0.05`, matched
  neither the code nor Eq. (20).
- `gotham.cpp` inherited this in 4de5ecede.
- The comment now states what the code computes.

**Fix, not applied:** `Real tau_coag_Myr = 0.27;`. This makes coagulation 10×
faster.

**Other differences (modeling choices, not errors).**
- Dubois et al., following Aoyama et al., coagulate only the fraction F of
  gas with n > 0.1 H cm^-3, T < 1e4 K and an unresolved Jeans length. They
  evaluate it at a subgrid n = 1e3 H cm^-3 with σ = 0.1 km/s.
- GOTHAM coagulates every cell at its resolved n_H with v = 0.1 km/s, and
  still divides by F = 0.5.

## Thermal sputtering: 2.2× too slow

Code (`gotham.cpp`, section 1 of `DustSource`; `nH = X rho / m_u` with
X = 0.75):

```cpp
Real tau_sput_Myr = 1657.0 * (1.0e-3 / nH);
const Real T_ratio = T_K / 2.0e6;
tau_sput_Myr *= 1.0 + 1.0 / (T_ratio * T_ratio * sqrt(T_ratio));
...
const Real req = dt_Myr * 3.0 * rho_Ds / (tau_sput_Myr * a_small);
```

The grain-radius timescale this implies is

    a/|da/dt| (code) = 1657 Myr (a/1 µm) (1e-3 cm^-3 / n_H) [1 + (2e6 K/T)^2.5]

and the mass timescale is a third of that.

Tsai & Mathews 1995, ApJ 448, 84, Eqs. (14)-(15), as written by McKinnon et
al. 2017, MNRAS 468, 1505 (arXiv:1606.02714), `draft.tex` lines 243-262:

```latex
\frac{\diff a}{\diff t} = -(3.2 \times 10^{-18} \, \text{cm}^{4} \, \text{s}^{-1}) \left( \frac{\rho}{m_\text{p}} \right) \left[ \left( \frac{T_0}{T} \right)^\omega + 1 \right]^{-1},
\tau_\text{sp} = a \left| \frac{\diff a}{\diff t} \right|^{-1} \approx (0.17 \, \text{Gyr}) \left( \frac{a_{-1}}{\rho_{-27}} \right) \left[ \left( \frac{T_0}{T} \right)^\omega + 1 \right],
```

Here ω = 2.5, T_0 = 2e6 K and ρ_-27 = ρ / 1e-27 g cm^-3, "which corresponds
to an effective number density of n ≈ 6 × 10^-4 cm^-3" (that is ρ/m_p).

- **The 1657 Myr comes from this formula.** At ρ = 1e-27 g cm^-3 the Tsai &
  Mathews rate gives 1e-4 cm / (3.2e-18 cm^4 s^-1 × 5.98e-4 cm^-3) =
  1656 Myr per µm (0.17 Gyr per 0.1 µm).
- **The code attaches it to the wrong density.** It uses n_H = 1e-3 cm^-3,
  but with the code's own conversion n_H = X ρ / m_u, ρ = 1e-27 g cm^-3 is
  n_H = 0.75 × 1e-27 / 1.6605e-24 = 4.52e-4 cm^-3.
- **Hence the factor.** The code's timescale is 1e-3 / 4.52e-4 = 2.21× the
  source's at every n_H and T. For a = 0.1 µm, T >> 2e6 K and
  n_H = 1e-3 cm^-3 that is 165.7 Myr against 74.8 Myr.

**History.** Before 5b4d5c5c7 the prefactor was 102 Myr, which was 7.3× too
fast. 5b4d5c5c7 set it to the Tsai & Mathews value at ρ = 1e-27 g cm^-3 but
kept the 1e-3 cm^-3 reference density.

**Fix, not applied:** `1657.0 * (4.52e-4 / nH)`. This makes sputtering 2.2×
faster. Dubois et al. (Sect. 3.3) found a 10× longer sputtering time made no
noticeable difference in their ISM-focused galaxies. Whether 2.2× matters for
GOTHAM's dust in the hot CGM has not been checked.

## Reproduce the ratios

The expressions are transcribed from the code and the papers:

```python
X, m_u, m_p, Myr = 0.75, 1.66053907e-24, 1.67262192e-24, 3.15576e13

def t_coag_code(nH, Ds, a=0.003, s=3.0, F=0.5, v=0.1):
    return 2.7*(s/3.0)*(1.0e3/nH)*(0.01/Ds)*(0.1/v)/F*(a/0.005)
def t_coag_dubois(nH, Ds, a=0.003, s=3.0, F=0.5, v=0.1):         # Eq. (20)
    return 0.27*(a/0.005)*(s/3.0)/F*(Ds/0.01)**-1*(nH/1e3)**-1*(v/0.1)**-1
def t_shat_code(nH, Dl, a=0.1, s=3.0, v=10.0):
    return 540.0*(1.0/nH)*(s/3.0)*(0.01/Dl)*(10.0/v)*a
def t_shat_dubois(nH, Dl, a=0.1, s=3.0, v=10.0):                 # Eq. (18)
    return 54.0*(a/0.1)*(s/3.0)*(Dl/0.01)**-1*(nH/1.0)**-1*(v/10.0)**-1
def t_sput_code(nH, a):                                         # a/|da/dt|, T >> T_0
    return 1657.0*(1.0e-3/nH)*a
def t_sput_tm(nH, a):                                           # Tsai & Mathews Eq. (14)
    return (a*1e-4)/(3.2e-18*nH*m_u/(X*m_p))/Myr

for nH in (1e-3, 1.0, 1e3):
    print(nH, t_coag_code(nH, 0.01)/t_coag_dubois(nH, 0.01),     # 10.0
          t_shat_code(nH, 0.01)/t_shat_dubois(nH, 0.01),         # 1.0
          t_sput_code(nH, 0.1)/t_sput_tm(nH, 0.1))               # 2.215
```
