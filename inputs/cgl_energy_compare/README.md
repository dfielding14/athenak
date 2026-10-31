# Energy comparison plots

Run the plotting script with the matching MHD and CGL histories.

```sh
python3 scripts/plot_cgl_energy_compare.py \
  --mhd-history /path/to/mhd/run.mhd.hst \
  --cgl-history /path/to/cgl/run.mhd.hst \
  --output-dir /path/to/plots
```

Keep each run's companion `run.user.hst` beside its `run.mhd.hst`. The supplied
`turb` inputs enable user history, whose `B^2` column is the volume integral
of the squared **cell-centered** magnetic field. The script uses half of this
quantity for magnetic energy and subtracts it and kinetic energy from total
energy to obtain thermal energy. The cell-centered field is also used by the
EOS when recovering thermal energy. The ordinary `1-ME`, `2-ME`, and `3-ME`
columns instead average squared face fields and cannot be used for this
thermal-energy subtraction.

Paper and cooling histories supply `therm_cgl` directly. For an adiabatic
passive run, `therm_mhd` is the thermal energy associated with the independent
gas total energy. Its passive CGL thermal energy belongs to a separate
reservoir. Direct thermal histories take precedence over magnetic subtraction.
When direct thermal and magnetic measurements are both available, the script
uses each directly, retaining magnetic energy even when it is tiny compared
with thermal energy.

The script matches samples by time and removes abandoned rows after a restart
or a duplicate final write, keeping the later branch. Each MHD sample needs
a user-history sample at the same recorded time. Extra user-history samples
are allowed. Within each file the column names must agree across restart headers.

An older history lacking a direct thermal sum or the cell-centered magnetic
sum requires matching saved snapshots to recover exact thermal energy. The
script reports the missing data. The CSV and magnetic-energy curves use the
same physical quantities as the thermal-energy calculation.
