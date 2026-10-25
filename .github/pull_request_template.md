## What this changes

<!-- One or two sentences. What does this pull request do, and why? -->

## Type of change

- [ ] Bug fix
- [ ] New feature or mechanism
- [ ] Performance (including a native port)
- [ ] Documentation or tooling
- [ ] Other:

## Related issue

<!-- Link the issue, for example "Closes #12". Large changes should have an issue first. -->

## How it was checked

- [ ] `pytest` passes
- [ ] `APEX_USE_CPP=1 pytest` passes (or the native module could not be built: say so below)
- [ ] Tests added or updated, and a fix has a test that fails without it
- [ ] Parameters changed: `scripts/gen_params.py` was re-run

<!-- Paste the relevant output, or the numbers measured on the reference models in res/test_models/. -->

## Measurements (if the result of the slicing changes)

<!-- Which models, which metric, before and after. A mechanism that is not measured on all four reference models must stay opt-in. -->

## Licence and attribution

- [ ] My contribution is under GPL-3.0, or the source and its licence are named below
- [ ] New third-party libraries are added to the credits in the README
- [ ] Notices in `pipeline/bird/` are kept intact

## Notes for the reviewer

<!-- Anything untested, skipped, or uncertain. Say it here rather than leaving it to be found. -->
