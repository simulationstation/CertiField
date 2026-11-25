# Decisions

Decisions taken under prompt_benchpack_and_audit.md, with the reason.

- 2026-11-24. Day 0, stop and archive: no worker was running (the previous session had ended); the unpushed work of the stopped tracks was committed and pushed to their own branches, nothing merged or deleted, and archive/STATE_*.md written. No Jakobson-Naud or fractal-uncertainty track exists in this repository or this session, so there was nothing to stop or archive for it.
- 2026-11-24. `outreach/` and `audit/drafts/` are git-ignored: the repository is public, and the prompt forbids outside actions; drafts stay on this machine for the owner.
- 2026-11-24. The pack is built in Python (reusing Reference/Curved3D and Reference/Elasticity) in new folders only; the paclet and its frozen Skyrme code are not touched.
- 2026-11-24. Owner's change of plan: no time caps or stop rules; run every item of the prompt to the end, then the optional ones in the order P1c, stretch targets (P1 to 1e-4, P2 Q2 to 1e-3), R and Julia in the audit, the NDSolve comparison, P3 below 0.1 %. Milestone reports are posted without waiting; only irreversible steps (deleting work, publishing, contacting anyone) would be asked first. Assignment: builder (P1, P2, P1c, stretch, P3 tightening), harness (codes, NDSolve), auditor (all targets, then R and Julia), verifier (started when the certificate format is published), coordinator (monitor, packaging, final report).
