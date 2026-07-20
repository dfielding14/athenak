<!-- BEGIN build-memory-table -->
# Frontier Control-Plane Navigation

## Scope and safety boundary

This package is the fail-closed control plane for authorized Frontier PIC
campaigns. Production mutations must run through an installed, checksummed,
read-only snapshot via `run_control_plane.py`. A mutable source checkout may
install or capture explicitly source-only evidence where the runner permits it;
it must not create real manifests, ledger entries, reservations, submissions,
or reconciliation records.

Read `README.md`, the parent `../AGENTS.md`, and the active readiness policy
before changing or operating this package.

## Task lookup

| Task | Start in | Also inspect or validate |
| --- | --- | --- |
| Change installed-inventory verification or command dispatch | `run_control_plane.py` | `install_control_plane.py`, prepared inventory, and source-checkout rejection tests |
| Change shared paths, site/account policy, or fail-closed checks | `control_plane_common.py` | Every caller plus readiness registries and schemas |
| Change installation or storage preflight | `install_control_plane.py`, `capture_storage_preflight_evidence.py` | Installed checksums, authorized roots, and storage schema |
| Change ledger format or locking | `ledger.py`, `initialize_frontier_ledger.py` | Mirroring, append-only semantics, recovery, and ledger tests |
| Change reservation or submission | `create_pre_submit_manifest.py`, `validate_and_reserve_frontier_job.py`, `launch_trampoline.py` | Shell launch wrappers, attestation, ledger transitions, and exact policy bindings |
| Change post-run reconciliation | `reconcile_frontier_job.py`, `reconcile_*_registered_execution.py` | External artifact inventory, checksums, and claim-specific readiness records |
| Change build/freeze or policy promotion | `write_orion_build_profile.py`, `create_clean_candidate_freeze.py`, `revalidate_clean_candidate.py`, `promote_active_policy.py` | Schemas, clean-HEAD binding, and installed control-plane version |

## Validation

Run source-level tests from the repository root on a Linux/Frontier-compatible
host without invoking production mutations:

- `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=. python3 tst/publication/frontier_control_plane/test_control_plane.py`
- `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=. python3 tst/publication/frontier_control_plane/test_ledger.py`
- `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=. python3 tst/publication/frontier_control_plane/test_q011_pressure_review_packet_verifier.py`

The control-plane and ledger suites exercise Linux file sealing, canonical-path,
and locking behavior; macOS is not a valid substitute for those platform tests.

## Local constraints

- Preserve exact authorized roots, account/site rules, lock ordering, mirrored
  ledger semantics, and reason strings unless the governing policy changes.
- Fail closed on dirty or mismatched source, mutable installed files, checksum
  drift, stale reservations, missing evidence, or ambiguous state transitions.
- Do not demonstrate a mutating command from the source tree. Use fixtures and
  temporary roots in tests; production operation requires the installed runner.
- Keep campaign-specific reconcilers distinct unless their schemas and claim
  contracts are demonstrably identical.
<!-- END build-memory-table -->
