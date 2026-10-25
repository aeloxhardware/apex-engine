# Contributing to Apex Algorithm

Thanks for your interest. Apex Algorithm is a pre-release slicing engine, so interfaces and defaults still move. Please read this page before opening a pull request.

## Ways to contribute

- **Report a bug or suggest an idea:** open an issue. For a bug, include the model (or a minimal one that shows it), the command you ran, the output, and your OS and Python version.
- **Propose a change:** fork the repository, make your change on a branch, and open a pull request against `main`. Nothing is merged without review.
- **Ask a question:** open an issue and say it is a question.

For anything large (a new step, a change of data contract, a new dependency), open an issue first so we can agree on the approach before you write the code.

## Setting up

```bash
python -m venv .venv
.venv/bin/pip install trimesh scipy tetgen rtree manifold3d networkx shapely pytest
```

On Windows use `.venv\Scripts\python.exe`. `pymeshfix` is optional. The native C++ core is optional too: `python scripts/build_cpp.py` builds it, and nothing requires a compiler to run the tests.

## Before you open a pull request

1. **Run the tests in both positions of the native switch:**
   ```bash
   pytest
   APEX_USE_CPP=1 pytest
   ```
   The second run needs the native module. If you cannot build it, say so in the pull request.
2. **Add or update tests.** Each step has its own file in `tests/`. Build test meshes on the fly with `trimesh.creation` instead of committing mesh files. A test for a fix should fail without the fix.
3. **Keep the code style of the file you edit.** Same naming, same comment density, same idioms.
4. **Keep each pull request focused.** One change per pull request is easier to review and to revert.

## Project conventions

- **Each step is a pure function.** `pipeline/stepN_*.py` exposes `run(input, params) -> output` with no hidden state. Shared data types live in `core/mesh_types.py`.
- **Parameters are dataclasses.** Anything exposed to the interface must serialise through `config/schema.py`. If you add a parameter, regenerate the generated files with `scripts/gen_params.py`, since a test fails when they are out of date.
- **The Python code is the reference for the native core.** A native port must return what Python returns, bit for bit or within a tolerance that you state and measure. Never compile with `/fp:fast`.
- **Measure before you conclude.** A mechanism that improves one indicator can leave the reference metric unchanged or worse. Report the numbers on the reference models in `res/test_models/`.
- **Do not turn a mechanism on by default** unless it has been measured on all four reference models. New mechanisms start opt-in.
- **Do not claim what you did not check.** If a test is skipped or a case is untested, say so.

## Licence and attribution

Apex Algorithm is released under the [GNU GPL 3.0](LICENSE). By contributing, you agree that your contribution is released under the same licence.

- Only submit code you wrote, or code whose licence is compatible with GPL-3.0. Name the source and its licence in the pull request.
- Keep the credit to [S4 Slicer](https://github.com/jyjblrd/S4_Slicer) by Joshua Bird (GPL-3.0) in the README and `CITATION.cff`. Code derived from it (branch `bird_old`) must keep its notice and attribution.
- Add every new third-party library to the credits table in the README, with its version and licence.

## Review

Maintainers review pull requests as time allows. Expect questions and requests for changes. A pull request may be declined if it does not fit the direction of the project, and this is not a judgement on its quality.
