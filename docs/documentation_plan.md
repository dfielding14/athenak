# Athena-K Documentation Project — Developer Brief

This document is your starting point for creating **modern, navigable, and aesthetically polished documentation** for the Athena-K AMR-MHD codebase.

Your primary goals:

1. **A clickable, high-level system flowchart** of Athena-K’s runtime and major components.
2. **Per-file reference pages** with short, accurate descriptions for every file in the tree.
3. **Narrative module guides** that explain the purpose, structure, and interactions of key subsystems.
4. **A migration guide for Athena++ users**, focusing on Athena-K’s use of Kokkos and `par_for`-style wrappers, along with patterns, gotchas, and idioms specific to this code.
5. **A polished, modern documentation site** that is easy to navigate, visually appealing, and sustainable to maintain.

---

## 0. Output Requirements

- **Docs site**: Fast, searchable, responsive, with dark mode, built from a reproducible toolchain.
- **Flowcharts**: At least one full runtime diagram plus focused sub-diagrams for key workflows (AMR, I/O, integrators, physics).
- **Per-file docs**: Every source and header file described in 1–3 sentences.
- **Module pages**: Concise, readable descriptions of major components (mesh, meshblock, integrators, par_for usage, AMR, I/O, boundary conditions, problem generators).
- **Athena++ → Athena-K migration guide**: Practical and specific, not generic Kokkos docs.
- **CI integration**: Automated build, link check, and doc coverage check.

---

## 1. Tooling

We will use:

- **Doxygen** to parse the C++ and generate XML and call graphs.
- **Sphinx** (with `breathe` and `myst-parser`) to consume Doxygen XML and render narrative + API docs together.
- **Mermaid** for curated, human-readable flowcharts.
- **Graphviz** for auto-generated call and include graphs.
- **Furo theme** for modern look and responsive design.
- **Optional CSS tweaks** to improve aesthetics (spacing, colors, code block style).

Deliverables should live in a `docs/` folder with:
- `Doxyfile`
- `conf.py` for Sphinx
- `requirements.txt`
- `Makefile` (or equivalent build script)

---

## 2. Overall Approach

### Phase 1 — Codebase Inventory
- Scan the entire tree to generate a list of all `.hpp` and `.cpp` files.
- For each file, capture:
  - Location in the source tree.
  - Top-level comment summary (if present).
  - Key types/functions defined.
  - Include dependencies.
- Use this inventory to seed per-file stub pages.

### Phase 2 — High-Level Flow & Modules
- Identify the main runtime sequence:
  1. Startup and parameter parsing.
  2. Mesh and MeshBlock creation.
  3. Initialization of physics/problem generators.
  4. Main time integration loop.
  5. AMR prolong/restrict cycles.
  6. Boundary communication and synchronization.
  7. I/O and finalization.
- Draft a **Mermaid flowchart** showing this at a high level.
- Create sub-diagrams for:
  - The AMR workflow.
  - The integrator and task execution model.
  - The I/O pipeline (including device ↔ host transfers).

### Phase 3 — Per-File Documentation
For each file:
- Write a **1–3 sentence summary** explaining its role in the codebase.
- Link to related modules or symbols.
- Avoid low-level detail dumps; focus on what a developer needs to know before opening the file.
- Include:
  - **Defines**: key classes, structs, functions.
  - **Depends on**: major includes and why.
  - **Used by**: key modules/functions.

### Phase 4 — Module Pages
Create narrative pages for:
- Mesh and MeshBlock.
- AMR (refinement criteria, prolong/restrict).
- Integrators (Runge–Kutta, VL2, etc.).
- Physics modules (how to add/extend).
- Reconstruction & Riemann solvers.
- Boundary conditions & exchanges.
- I/O and restart mechanisms.
- Parallelism and tasking (including MPI + Kokkos).
- `par_for` usage in Athena-K:
  - Show how it wraps Kokkos.
  - Discuss execution policies, indexing patterns, and best practices for GPU portability.
  - Highlight differences from Athena++ loop structures.

### Phase 5 — Athena++ → Athena-K Migration Guide
- Start with a mapping table of Athena++ concepts → Athena-K equivalents.
- Focus on:
  - `par_for` in Athena-K vs raw loops in Athena++.
  - Memory management via Kokkos Views.
  - Execution policies and how they are hidden/abstracted in Athena-K.
  - Boundary between host and device code.
  - Common mistakes when porting code (e.g., capturing host-only data in device lambdas, missing fences before MPI).
- Include **side-by-side** code snippets:
  - Athena++ triple loops → Athena-K `par_for`.
  - Stencil operations with halos.
  - Reductions.
- Include **gotchas**:
  - Correct index ordering.
  - LayoutLeft/LayoutRight choices.
  - Temporary array handling.
  - Device ↔ host data movement.

Example snippet (adapt to Athena-K style):
```cpp
// Athena-K par_for example
par_for("UpdateFluxes", DevExeSpace(), 0, nk, 0, nj, 0, ni,
  KOKKOS_LAMBDA(const int k, const int j, const int i) {
    // loop body
});
````

### Phase 6 — Aesthetics & Navigation

* Use Furo theme with dark mode.
* Organize nav bar:

  * **Overview**
  * **Flowcharts**
  * **Modules**
  * **Per-File Reference**
  * **Migration Guide**
* Make all diagrams clickable (nodes link to relevant doc pages).
* Ensure syntax highlighting matches Athena-K’s coding style.
* Add consistent icons/callouts for:

  * Warnings (“gotchas”).
  * Tips.
  * Examples.

### Phase 7 — CI & Maintenance

* Add GitHub Actions to:

  * Build docs.
  * Check for broken links.
  * Enforce doc coverage for source files.
* Provide a short README in `docs/` with build instructions.

---

## 3. Style Guide

* Keep summaries short and factual.
  * Concise, factual, “engineering handbook” voice. No fluff. Avoid hype.
  * Use consistent terminology (Athena-K preferred names for modules/types).
  * Prefer short paragraphs and bullet lists.
* Avoid duplicating code; link to API pages instead.
* Use bullet points for lists of responsibilities or gotchas.
* Keep diagrams **readable**: no more than 8–10 nodes per diagram unless it’s a top-level overview.
* When in doubt, prefer clarity over completeness — the goal is navigation and understanding, not dumping every line.
* Doc anatomy
  * Each page starts with a 1–2 sentence summary.
  * Then context (where this fits in the runtime).
  * Then details (key types/functions, invariants, device/host aspects).
  * End with links to related modules/files and code references.

---

## 4. Flexibility for the Developer

These instructions set **direction and constraints**, but **you have freedom** to:

* Decide the best diagram granularity.
* Refine which files group naturally into modules.
* Adjust formatting and navigation if it improves clarity or aesthetics.
* Propose additional tooling or integrations if it will improve the result.

The key measure of success:
A new Athena-K contributor can start at the top-level diagram and, within 3–4 clicks, get to a clear explanation of **any** part of the code they need to modify.

---
## 5. Per-Module and Per-File Documentation Procedure

For **every module** and **every file** in the Athena-K codebase, follow a consistent, repeatable process.
Your first task is to design **robust, reusable templates** for:

1. **File pages** — concise summary, key types/functions, dependencies, and notes.
2. **Module pages** — high-level purpose, data flows, key APIs, extension points, gotchas.
3. **Side-by-side migration snippets** — Athena++ code vs. Athena-K equivalent.

Once the templates are approved, apply them systematically.

---

### 5.1 Step-by-Step: File Documentation

1. **Locate file** in the source tree and open it.
2. **Identify purpose**: summarize in 1–3 sentences.
3. **List key symbols**: public classes, structs, functions.
4. **Dependencies**: major includes and why they are needed.
5. **Relationships**: where this file is used (major modules/functions).
6. **Notes**: gotchas, common pitfalls, or special considerations.
7. Fill in the **File Page Template** with this information.

---

### 5.2 Step-by-Step: Module Documentation

1. **Define scope**: list all files that make up this module.
2. **Purpose**: one-paragraph overview of what this module does in Athena-K.
3. **Place in runtime**: when/where this module is active (link to flowchart node).
4. **Key data structures**: describe ownership, lifetimes, host/device usage.
5. **Key functions/kernels**: purpose, policies, and where they are invoked.
6. **Extension points**: how a developer could extend or replace functionality.
7. **Gotchas**: ordering, synchronization, memory layout, or portability concerns.
8. Populate the **Module Page Template** accordingly.

---

### 5.3 Step-by-Step: Side-by-Side Migration

1. Choose a representative **Athena++ pattern** relevant to this module/file.
2. Locate the **Athena-K equivalent** in current code.
3. Show the **Athena++ version** (condensed but accurate).
4. Show the **Athena-K version**, using its actual `par_for` style and indexing.
5. Annotate both sides with:
   - What changed structurally.
   - What changed in parallel execution policy.
   - Any data layout or memory management differences.
6. Use the **Migration Snippet Template** to keep presentation consistent.

---

### 5.4 Expectations

- **Templates must be well-structured** so that future contributors can easily add to or update docs.
- Apply templates to **100% of modules and files**.
- Keep formatting and terminology consistent across all pages.
- Link every documented file and module into the site navigation and flowcharts.
