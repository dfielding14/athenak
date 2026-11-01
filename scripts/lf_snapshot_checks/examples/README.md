# Smoke data provenance

These two native active CGL $24^3$ outputs are from periodic matched CPU pilots at $t=3$, with initial $\beta_i=16$, $\rho=1$, $B_0=\sqrt{1/8}$, and $p_\parallel=p_\perp=1$. Both used the same forcing seed 314159, finite soft rate 64, backup rate $10^{10}$, and history cooling. Background collision rates were zero and 20; the second has $\nu_{\rm coll}L/v_{A0}=56.57$. The output fields are float32 from a double-precision executable. The reached $u_{\rm rms}/v_{A0}$ was 1.49 and 1.51. The executable SHA-256 was `46f4f35c3a4636d54d731e02c488a528201a13012a987ad8a58e7d2cbb0dfe05`. The full input is embedded in each binary.

The small reference arrays are fixed-$k_L$ and finite-reach signed fluxes and conservative tendencies extracted from the independent all-face pilot analysis (reaches $0.125L,0.25L,0.5L$ at ray step $\Delta x/4$). The validation points include periodic boundaries. These snapshots are not production evidence.
