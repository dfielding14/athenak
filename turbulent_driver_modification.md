# Spec: Concurrent Turbulence Driving (Impulsive + Sustained) in Athena‑K

**Owner:** @<your-handle>
**Audience:** AI code agent (and future maintainers)
**Goal:** Allow *simultaneous* operation of two turbulence driving channels with **different properties**:
1. **Impulsive / decaying** driving at the beginning of the run (a one‑off “kick”; currently **type 0**).
2. **Sustained** driving throughout the simulation (continuous in time).

> ⚠️ **Do not guess. Ask first.** Before writing code, confirm all items in **§0. Required confirmations** with the Owner. Be surgical and minimally invasive. Preserve existing behavior when the new feature is not used. Document the change clearly here and in `~/Work/Research/AthenaK/athenak-df/docs/source/modules/srcterms.md`.

---

## 0) Required confirmations (ask the Owner first)

1. **Which ID denotes “continuous”?**
   The current code documents: type **0** = impulsive & decaying, type **1** = finite‑duration then stop, type **2** = continuous.
   The request says “do type **0 and 1** simultaneously” *and* mentions “continuous driving.”
   - **Confirm** whether the second driver should be **type 2 (continuous)** or **type 1 (finite duration)**.
   - If it is type 1, **confirm duration** and stop time semantics.

2. **Config format currently in use** (e.g., Athena++-style `input.par` blocks vs. YAML/TOML).
   - **Confirm** how the turbulence driver is currently parameterized and where (block/section names, keys).

3. **Spectral & forcing details** to be independent per driver:
   - Solenoidal/compressive mix (`f_sol` or equivalent).
   - Spectrum parameters (`kmin`, `kmax`, `kpeak` or equivalent).
   - Target amplitude / power injection / `rms_accel` (whichever the code uses).
   - Correlation time `t_corr` (for OU/continuous).
   - Time schedule (`t_start`, optional `t_stop`).
   - Optional exponential decay constant for type 0 (if supported/desired beyond “hydro/MHD decay”).

4. **RNG & determinism requirements** (restarts, MPI/GPU invariance).
   - **Confirm** that each driver needs its own seed and that reproducibility across MPI ranks/GPU is required.

5. **Diagnostics to expose**
   - Per‑driver injection power, amplitude, seeds, and spectral band in logs?
   - Per‑driver contributions in output fields/diagnostics?

6. **Restart I/O**
   - **Confirm** whether to persist per‑driver internal states (e.g., OU process state) in restart dumps.

---

## 1) Capability to add

### 1.1 Summary
Add **multi‑driver support** for turbulence source terms so the code can hold **N ≥ 1** independent driver instances and **sum** their contributions every step. Preserve backwards compatibility:
- If the user configures a single driver as before, behavior is **identical**.
- If multiple drivers are configured, the total forcing is the **sum of each driver’s forcing** (applied in deterministic order).

### 1.2 Minimal design
Introduce a **thin composite wrapper** that reuses existing driver implementations:

```

+-------------------------+

| CompositeTurbDriver         |                                                  |
| --------------------------- | ------------------------------------------------ |
| std::vector<DriverPtr>      | ----> existing Type0 driver (impulse)            |
|                             | ----> existing Type1 or Type2 driver (sustained) |
| methods mirror current      | (no behavior change to existing classes)         |
| +-------------------------+ |                                                  |

````

- The composite obeys the *current* turbulence driver public interface (constructor signature, `Init(...)`, `ApplyForcing(...)` or equivalent), so **call sites remain unchanged.**
- Summation order is **stable** (index order) for deterministic results.
- Each child driver has its **own seed, spectrum, amplitude, schedule**, etc.

---

## 2) User‑facing behavior

### 2.1 What the user can do after this change
- Configure **two concurrent drivers with different properties**:
  - **Driver A:** *impulsive* at `t_start_A` (one‑shot), then passively decays in the fluid (and/or optional exponential decay if supported).
  - **Driver B:** *sustained* in time (continuous OU type or finite‑duration type per Owner confirmation).
- Allocate fully independent **spectral windows** and **forcing modes** to each driver.
- Use **separate RNG seeds** per driver for reproducibility.
- See **per‑driver diagnostics** (init summary in logs; optional runtime power injection).

### 2.2 Backwards compatibility
- Existing inputs continue to work with **no changes**.
- If the config contains only a single driver definition, the composite wrapper collapses to that single driver.

---

## 3) Configuration schema (proposed)

> **Important:** Use the project’s existing configuration pattern. If blocks are named differently, map these keys accordingly and keep legacy keys untouched.

Two *equivalent* minimal, backward‑compatible approaches are acceptable—choose the one that matches the current parser:

**A) Indexed blocks (preferred if multi‑section parsing is already supported):**
```ini
# Example (adjust block names/keys to current repo conventions)
[turbulence.drivers.0]
type        = 0        # impulsive
t_start     = 0.0
# Optional—only if supported beyond passive decay:
decay_model = "none"   # or "exp"
tau_decay   = 0.0
kmin        = 2
kpeak       = 3
kmax        = 4
f_sol       = 1.0
amplitude   = 0.5
seed        = 12345

[turbulence.drivers.1]
# CONFIRM with Owner: type = 2 (continuous) or type = 1 (finite duration)
type        = 2
t_start     = 0.0
# For type 1 only:
# t_stop   = <float>
t_corr      = 0.5
kmin        = 2
kpeak       = 3
kmax        = 4
f_sol       = 1.0
amplitude   = 0.05
seed        = 67890
````

**B) Single block with repeated, suffixed keys (use only if A is infeasible):**

```ini
[turbulence]
n_drivers             = 2

driver0.type          = 0
driver0.t_start       = 0.0
driver0.kmin          = 2
driver0.kpeak         = 3
driver0.kmax          = 4
driver0.f_sol         = 1.0
driver0.amplitude     = 0.5
driver0.seed          = 12345

driver1.type          = 2      # or 1, per Owner
driver1.t_start       = 0.0
driver1.t_corr        = 0.5    # for continuous OU
driver1.kmin          = 2
driver1.kpeak         = 3
driver1.kmax          = 4
driver1.f_sol         = 1.0
driver1.amplitude     = 0.05
driver1.seed          = 67890
# driver1.t_stop     = <float> # if type 1
```

**Notes**

* **Seeds:** If a seed is omitted, derive a deterministic default as `global_seed + index * large_prime`.
* **Amplitude / Power:** Use existing code’s controls (e.g., `rms_accel`, target energy injection); keep names unchanged.
* **Schedules:** Respect `t_start`; for **type 1** also respect `t_stop`; for **type 0**, apply once at/after `t_start` then stop injecting.

---

## 4) Implementation plan (surgical)

1. **Locate current turbulence driver**

   * Identify the *existing* driver class (e.g., `TurbDriver`, `TurbulenceSrcTerm`, etc.) and its constructor & apply methods.
   * Identify the parsing site for turbulence driver parameters.

2. **Add composite wrapper** (new ~100–150 LOC total)

   * New class `CompositeTurbDriver` in the same module as the current driver (or adjacent).
   * Holds `std::vector<std::unique_ptr<BaseTurbDriver>> drivers_;` (use the current base type).
   * Implements the **same** public API as the current single driver.
   * `ApplyForcing(...)` simply loops `for (auto& d : drivers_) d->ApplyForcing(...);`

     * Ensure **GPU safety**: if existing drivers already offload kernels, the composite loop must not introduce extra device<->host transfers. Prefer fusing where trivial; otherwise keep the minimal loop and rely on existing kernels.
   * The composite replaces the single driver at **exactly one** instantiation point (factory / constructor), leaving the rest untouched.

3. **Parsing multiple drivers**

   * Extend the existing parsing function to:

     * Detect multiple driver blocks or suffixed keys (Approach A or B in §3).
     * Build one `BaseTurbDriver` instance per block and emplace into `CompositeTurbDriver`.
     * If only **one** driver is configured, still construct the composite with a single element (or conditionally skip the composite but expose the same type).
   * Keep **all legacy keys** working. No renames.

4. **Per‑driver scheduling**

   * Reuse existing logic:

     * **type 0:** create once on/after `t_start`; no further updates. (If exponential decay exists in the codebase and Owner wants it, add `decay_model` and `tau_decay` gated behind config; default is “none”.)
     * **type 1:** active only for `t ∈ [t_start, t_stop)`.
     * **type 2:** continuous for `t ≥ t_start` with OU correlation `t_corr`.

5. **RNG separation & determinism**

   * Each driver uses an independent RNG stream.
   * If the repository has a *global* RNG, derive per‑driver sequences deterministically (document the scheme).
   * Verify restarts reproduce bitwise-identical results for fixed seeds, MPI layout, and GPU/CPU (if that is already true today).

6. **Diagnostics & logging**

   * On init, print **one line per driver**: index, type, k‑band, f_sol, amplitude or injection target, `t_start`, `t_stop` (if any), `t_corr` (if any), and seed.
   * (Optional) If there is an existing per‑step “power injected” diagnostic, add **per‑driver** contributions and the sum. Make it **opt‑in** via a new config flag to avoid overhead.

7. **Restart I/O**

   * If the current driver persists internal states (e.g., OU fields), persist **each driver’s state** under a namespaced key (e.g., `/turbulence/driver/0`, `/1`, …).
   * Bump a **minor** version in the restart metadata if needed.

8. **Unit & integration tests**

   * **Construction test:** Multiple drivers parse correctly; logs match expectations.
   * **Sum correctness (determinism):** For small 3D grids with fixed seeds, verify the composite + summed single runs are identical to machine precision where applicable.
   * **Schedule edges:** `t == t_start`, `t == t_stop`.
   * **Restart round‑trip:** Write → read → continue produces identical results vs. uninterrupted run.

---

## 5) Code sketch (illustrative, adapt to repo)

> Use existing types/names; this is **pseudocode** to anchor the approach.

```cpp
// CompositeTurbDriver.hpp
class CompositeTurbDriver : public BaseTurbDriver {
 public:
  explicit CompositeTurbDriver(std::vector<std::unique_ptr<BaseTurbDriver>> drivers)
  : drivers_(std::move(drivers)) {}

  void Initialize(const GridMeta& grid) override {
    for (auto& d : drivers_) d->Initialize(grid);
  }

  // Called each step to add body force / source terms
  void ApplyForcing(State& S, Real time, Real dt) override {
    for (auto& d : drivers_) d->ApplyForcing(S, time, dt);
  }

  void WriteRestart(RestartWriter& w) const override {
    for (size_t i = 0; i < drivers_.size(); ++i) {
      w.PushGroup(fmt::format("turbulence/driver/{}", i));
      drivers_[i]->WriteRestart(w);
      w.PopGroup();
    }
  }

  void ReadRestart(RestartReader& r) override {
    for (size_t i = 0; i < drivers_.size(); ++i) {
      r.PushGroup(fmt::format("turbulence/driver/{}", i));
      drivers_[i]->ReadRestart(r);
      r.PopGroup();
    }
  }

 private:
  std::vector<std::unique_ptr<BaseTurbDriver>> drivers_;
};
```

**Factory wiring (single touchpoint):**

```cpp
std::unique_ptr<BaseTurbDriver> MakeTurbDriver(const Params& P) {
  auto blocks = ParseDriverBlocks(P);  // returns N configs (N>=1)
  std::vector<std::unique_ptr<BaseTurbDriver>> ds;
  ds.reserve(blocks.size());
  for (size_t i = 0; i < blocks.size(); ++i) {
    ds.emplace_back(MakeSingleDriver(blocks[i], /*index=*/i));
  }
  if (ds.size() == 1) return std::move(ds[0]); // optional micro‑opt
  return std::make_unique<CompositeTurbDriver>(std::move(ds));
}
```

---

## 6) Performance & safety

* **Performance:** The composite adds only a small loop over *N* drivers. Each driver’s existing kernels remain unchanged. If profiling shows launch overhead for many drivers, consider trivial fusion; **not required** for N=2.
* **Thread/MPI/GPU safety:** No new global state. All per‑driver state lives inside its instance. Summation order is fixed.
* **Numerical stability:** Since the total forcing increases, document that users should retune amplitudes to keep the desired Mach/Alfvénic regime.

---

## 7) Documentation updates

Update `~/Work/Research/AthenaK/athenak-df/docs/source/modules/srcterms.md`:

* **New subsection:** *“Composite Turbulence Driving”*

  * Explain that multiple drivers can be active simultaneously and their forces are summed.
  * Document configuration (Approach A or B in §3), including a **worked example** showing a type 0 impulse plus a sustained driver.
  * Clarify the semantics of types (0, 1, 2) and scheduling with `t_start`, `t_stop`, `t_corr`.
  * Note reproducibility considerations (per‑driver seeds) and restart behavior.
  * Add a short **FAQ**: “Can I use more than two drivers?” (Yes), “What if spectra overlap?” (They sum; tune amplitudes), “How do I keep old runs identical?” (Use a single driver exactly as before).

**Doc snippet to include verbatim (adapt names to actual parser):**

````markdown
### Composite Turbulence Driving

Athena‑K supports multiple turbulence driving channels **active at the same time**. Each channel can be
configured independently (type, spectrum, solenoidal fraction, amplitude/power, correlation time, schedule, seed).
At each step the body force is the **sum** of all active channels.

**Example: impulsive kick at start + continuous forcing**
```ini
[turbulence.drivers.0]
type        = 0
t_start     = 0.0
kmin        = 2; kpeak = 3; kmax = 4
f_sol       = 1.0
amplitude   = 0.5
seed        = 12345

[turbulence.drivers.1]
type        = 2            # continuous OU
t_start     = 0.0
t_corr      = 0.5
kmin        = 2; kpeak = 3; kmax = 4
f_sol       = 1.0
amplitude   = 0.05
seed        = 67890
````

> Backwards compatibility: Existing single‑driver inputs behave identically to previous versions.

```

---

## 8) Testing checklist

- [ ] Parser accepts **two** driver definitions and rejects malformed ones with clear errors.
- [ ] Logs print the **two** drivers with independent properties.
- [ ] For a small 3D box and fixed seeds, composite run equals the sum of separate runs to within round‑off when fields are added offline (sanity check).
- [ ] Type 0 fires **once** at `t_start`, then no further updates (unless `decay_model` explicitly enabled).
- [ ] Type 1 respects `t_stop`. Type 2 is continuous after `t_start`.
- [ ] Restart from a checkpoint reproduces the trajectory (if the single‑driver version already did).
- [ ] Optional per‑driver diagnostics emit when enabled.

---

## 9) Deliverables

1. **Code changes** (minimal surface area):
   - New composite class in the turbulence source‑term module.
   - Parser extension for multiple driver blocks/keys.
   - Factory wiring change at a single call site.
   - Optional: per‑driver diagnostics gated by a config flag.

2. **Docs**:
   - Update `~/Work/Research/AthenaK/athenak-df/docs/source/modules/srcterms.md` per §7.

3. **Tests**:
   - Add/extend unit tests for parser and composite behavior.
   - A small regression input that exercises the dual‑driver case.

---

## 10) PR hygiene

- One focused PR with clear title: “srcterms: add composite turbulence driving (impulsive + sustained)”.
- Explain **why** and **what changed**. Include the config example.
- Link to this spec. Summarize acceptance tests and results.
- Keep diffs small; avoid renames that churn blame history.

---

## 11) Open questions for the Owner (must be answered before coding)

- Confirm **type ID** for the sustained driver (**2** continuous vs **1** finite‑duration).
- Provide desired **default values** for the second driver’s amplitude and `t_corr` (or duration if type 1).
- Confirm **configuration style** (blocks vs suffixed keys) to match the repo and avoid surprises.
- Specify whether **exponential decay** for type 0 is desired beyond passive decay, and if so, provide `tau_decay`.
- Confirm the exact list of **diagnostics** you want to see and their output cadence.

---

## 12) Non‑goals (explicitly out of scope)

- Changing physics or math of existing single drivers.
- Altering legacy key names or breaking old input decks.
- Deep kernel fusion or performance retuning (unless profiling shows a clear regression).

---

*Please ask the Owner about any ambiguity. Do not guess. Keep the changes surgical, with a clear explanation in code comments and in the docs update.*
```

---

If you want, I can also generate a ready‑to‑paste doc patch for `srcterms.md` with exactly the subsection from §7 already formatted for your Sphinx build.
