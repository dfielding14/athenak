# Squire archived reference B4 experiment

Immutable archived source: `44a29a70f2869b44e97d725393d4b9a020a60adf`.

Archive command: `git --git-dir=/Users/dbf75/Work/Research/AthenaK/Athena++/AthenaCGL/.git_old_dev_version_history archive 44a29a70f2869b44e97d725393d4b9a020a60adf | tar -x -C /tmp/cgl-b4-compare-20261006/squire`.

Configure: `python3 configure.py --prob=b4_reference --coord=cartesian --eos=cgl --flux=hlle -b --cxx=clang++-apple`.

Build: `make -j2`. Commands and output are in configure.log, build.log, rebuild.log. Runs: `python3 run_cases.py`; full per-case commands/inputs/status are in summary.json and runs/*/.

Only custom pgen and read-only stage diagnostics were added. All archived physics files retain their exact bytes. Changed tracked source files: ['src/task_list/time_integrator.cpp']. The two stage calls dump after C2P and after collision/P2C, first five cycles only. Snapshot stage time is beginning-of-cycle time; cycle snapshots record end-of-cycle time.

## Actual magnetic floor behavior

The exact archive initializes Mesh::bsq_rms=0. Its src/hydro/new_blockdt.cpp:208 has `if (CGL_EOS && false)` around both an inactive field timestep constraint and the RMS accumulation (lines291–299); bsq_av initializes0 at65 and is assigned to block_bsq_rms at361. Thus global bsq_rms remains zero in these runs. GetBSqrFloor() returns zero regardless of configured bsqr_floor. Both EOS and face cutoffs are zero throughout, as recorded percell in CSVs. The nominal facefloor2e-20 and magnetized2e-28 controls are consequently identical. This behavior is preserved; no floor fix was applied. The initial physical mean B² of the sharp contact is .5 but is not reflected in archived bsq_rms.

The source additionally has an EOS-versus-face floor-unit inconsistency if RMS were nonzero: EOS compares |B| against bsqr_floor*bsq_rms, faces against its square root. This experiment does not exercise that branch because actual RMS stays zero. Read-only inspection of the local modified Athena++ tree confirms its src/hydro/new_blockdt.cpp:137–149 restores active RMS accumulation, assigning block_bsq_rms at293; that modified source was not run.

## Initial conditions and outputs

Sharp: rho=1, v=±10, Bx=Bz=0, By=1e-12 on incoming side and1 on the outgoing side, isotropic p=1.5−By²/2. The v<0 state is the reflected fixture, matching current AthenaK. PLM RK2 CFL.4,128 cells,one block,outflow, optional kinetic limiters and nu_coll off, no heat flux. There is no Squire archived FOFC implementation. Initialization calls unmodified PrimitiveToConserved with actual B; weak A/rho=82.89306334778564.

Smooth conditional follow-up: v=10, By=1e-12+(1−1e-12)*(1+tanh((x−.5)/.04))/2, p=1.5−By²/2; nx128/256/512, fixed t=.012. All states above actual zero floor. Raw percell everycycle data and first-five-stage dumps retain rho,momentum,E,A,pressures,B and actual floor diagnostics.

No files in the archived reference directory or production worktree were modified. Source/binary SHA256 manifest and experiment.patch are retained.
