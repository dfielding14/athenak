# AthenaK Documentation Audit Primer

Welcome to the AthenaK documentation audit project.  This primer summarizes the
expectations, workflow, and guardrails you must follow before touching any
documents.  Keep it beside you while you work and refer back to the detailed
instructions in the repository whenever you are unsure.

---

## 1. Read the Ground Rules First

Start with the existing guidance:

- `docs/documentation_audit_guide.md` – step‑by‑step verification workflow,
  required logging, and escalation instructions.
- `docs/documentation_audit_checklist.md` – running list of pages that still
  need to be vetted.
- `CODEBASE_OVERVIEW.md` – high-level architecture summary; use it as a map,
  but always confirm specifics in the source before editing docs.

Do **not** begin editing until you understand both files.

---

## 2. Absolute Expectations

- **Vet every statement against the code.** For each claim in the docs, identify
  the corresponding implementation (header, source file, input deck, build
  script, etc.) and confirm it matches current behavior.  Citations should point
  to the files you actually inspected.
- **Be 100 % certain before you commit a change.** If you cannot prove that a
  statement is correct from the repository contents, do not include it.  Remove
  stale information or leave TODOs rather than guessing.
- **Never guess or invent context.** If something is unclear, ask.  Escalate key
  questions instead of inferring how a subsystem “probably” works.
- **Respect the existing instructions.** They supersede any personal habits.

---

## 3. Typical Audit Workflow

1. **Pick an item** from `documentation_audit_checklist.md`.
2. **Read the relevant doc page** end‑to‑end to understand what it claims.
3. **Trace each statement** to the implementation:
   - Use `rg`, `sed`, and the existing citations to locate the code, inputs, or
     scripts.
   - Confirm defaults, parameter names, and behavior directly from the source.
4. **Update the documentation** if and only if you have confirmed the facts.
   - Provide updated references in the `[file:line]` format already used in
     neighboring pages.
5. **Record the work**:
   - Add an entry to `docs/documentation_audit_log.md`.
   - Check off the item in `docs/documentation_audit_checklist.md`.
6. **Rebuild the docs** (`cd docs && make html`) and resolve any warnings.

---

## 4. Communication & Escalation

- If you encounter missing context, unfamiliar subsystems, or ambiguous
  behavior, pause and request clarification.  Do **not** fill gaps with
  assumptions.
- Surface blockers early—especially when code behavior contradicts the
  documentation or when the instructions themselves appear stale.

---

## 5. Tooling & Conventions

- Work inside the provided repo tools (shell, `rg`, `sed`, etc.).  Avoid IDE‑only
  shortcuts so that every step can be reproduced.
- Stick to Markdown and existing style conventions (headings, tables, citation
  formats).  Follow the example set by recently audited files.
- Always run `make html` after edits to ensure Sphinx succeeds without warnings.

---

## 6. Final Checklist Before Submitting

- [ ] Every modified doc page has been cross‑checked against current code.
- [ ] Citations reference real files/lines that you personally inspected.
- [ ] `docs/documentation_audit_log.md` contains an entry summarizing your work.
- [ ] `docs/documentation_audit_checklist.md` reflects the updated status.
- [ ] Sphinx build passes cleanly.
- [ ] You are prepared to defend each change with specific code references.

Remember: the goal is to eliminate drift between the codebase and its
documentation.  Certainty matters more than speed.  When in doubt, stop and ask.***
