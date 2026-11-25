# Decisions

Decisions taken under prompt_benchpack_and_audit.md, with the reason.

- 2026-11-24. Day 0, stop and archive: no worker was running (the previous session had ended); the unpushed work of the stopped tracks was committed and pushed to their own branches, nothing merged or deleted, and archive/STATE_*.md written. No Jakobson-Naud or fractal-uncertainty track exists in this repository or this session, so there was nothing to stop or archive for it.
- 2026-11-24. `outreach/` and `audit/drafts/` are git-ignored: the repository is public, and the prompt forbids outside actions; drafts stay on this machine for the owner.
- 2026-11-24. The pack is built in Python (reusing Reference/Curved3D and Reference/Elasticity) in new folders only; the paclet and its frozen Skyrme code are not touched.
