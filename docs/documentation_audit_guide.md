# Documentation Audit Playbook

Use this playbook for every documentation update so future contributors follow the same workflow and nothing in the docs drifts away from the code.

---

## Core Principles
1. **Vet every statement against the code.** Never assume behaviour—verify it in the source tree, input decks, or build scripts.
2. **Never guess or speculate.** If you cannot confirm something from the repository, do not keep it in the docs.
3. **Use clear, direct language.** Prefer concise explanations over ambiguous wording.
4. **Ask for help when in doubt.** If scope or behaviour is unclear, flag it in the log, alert the maintainer, and move forward.

## Step 1 — Identify the Target
1. **Pinpoint the doc scope**  
   - Read the page or section under review.  
   - List every feature, option, command, parameter, or example it mentions.
2. **Map to the source tree**  
   - Search the repository (`rg`, `find`) to locate the implementation files and input decks the doc is referencing.  
   - If you cannot find a corresponding file/function/flag in the tree, flag it immediately—the doc may be stale.

## Step 2 — Cross-Check Against the Code
1. **Inspect the implementation**  
   - Read the relevant `src/` modules, problem generators, or scripts to confirm behaviour, defaults, and configuration options.  
   - Note any guardrails or fatal error conditions (e.g., `tile_nx` must divide `mesh.nx1` in the turbulence driver).
2. **Compare with the current documentation**  
   - Highlight mismatches (removed files, renamed options, additional parameters, changed defaults).  
   - Update the text, tables, examples, and code snippets so that terminology and values match exactly what the code uses today.  
   - Add missing explanations for supported features (e.g., tiled driving), and **remove** anything that no longer exists (obsolete modes, defunct parameters, speculative references).

## Step 3 — Prune Irrelevant Material
1. **Delete dead references**  
   - Remove bibliographic references, inputs, or workflow steps that cannot be traced back to code in the repository.  
   - Drop example commands that require build flags or generators that no longer exist.
2. **Keep examples executable**  
   - Ensure every command shown uses the supported out-of-source build (`cmake -S . -B build` / `cmake --build …`).  
   - Prefer runtime overrides (e.g., `mesh/nx1=512`) over editing input files in place.

## Step 4 — Update Progress Artifacts
1. **Checklist (`docs/documentation_audit_checklist.md`)**  
   - Mark the relevant checkbox as completed once the doc has been fully vetted.  
   - Leave unchecked or add notes for items still outstanding.
2. **Audit log (`docs/documentation_audit_log.md`)**  
   - Add a row summarising the change: date, file(s), issue, fix, and source line references (e.g., `src/srcterms/turb_driver.cpp:118`).  
   - This creates an audit trail for every documentation adjustment.

## Step 5 — Verify the Sphinx Build
1. Run `cd docs && make html` and ensure it exits cleanly (no warnings about missing lexers or broken references).
2. Spot-check the generated page (`docs/build/html/...`) to confirm the rendered output matches your intent.

---

### Quick Reference Template
1. Identify doc scope and map statements → code files.  
2. Confirm functionality and defaults in implementation; update wording/tables/snippets.  
3. Remove unsupported features/examples/references.  
4. Update checklist checkbox + append log entry with evidence.  
5. Rebuild docs (`make html`) and fix any warnings.

> **Escalation:** If a page is too inconsistent to repair quickly (missing code, unclear scope, or would require new assumptions), note the issue in the log, alert the maintainer, and move on rather than guessing.

Following this template keeps the documentation in lockstep with the codebase and saves the next auditor from rediscovering stale sections.
