# Contributing

Bug reports, questions and pull requests are welcome. Please open an issue first for larger changes so we can agree on the approach.

## Building and testing

See [Building](README.md#building) in the README. Before sending a change, build with warnings enabled (the default) and run the tests:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

```bash
cmake --build build
```

```bash
ctest --test-dir build --output-on-failure
```

CI runs the same steps on Linux, Windows and macOS, builds the GUI, and runs the Sod study of the analysis scripts.

## Guidelines

- **Keep it simple.** The code is meant to be read and learned from: prefer clear, direct code over clever abstractions.
- **Match the surrounding style**: braces on their own line, `if(` without a space, a short banner comment above each function that explains the method, not the syntax.
- **`src/core` stays free of Qt** and of third-party dependencies, so the solver, the CLI and the tests build anywhere with a C++17 compiler.
- **New numerics come with a test** in `tests/tests.cpp`: a consistency or conservation check, an exact solution, or an order-of-accuracy measurement.

## When results change

If a change affects the numerical results, regenerate the studies and commit the updated figures and tables together with the change:

```bash
python analysis/run_all.py
```

Delete `analysis/runs/` first, since cached runs are reused. Update the numbers quoted in the README if they changed.

## GUI screenshots

The screenshots in `docs/screenshots` are produced by a small program in [`tools/screenshots`](tools/screenshots/CMakeLists.txt) that drives the GUI. Regenerate them after visible GUI changes.

## Making a release (maintainers)

1. Update `CHANGELOG.md` and tag the commit (`git tag v1.x.y`).
2. Build the GUI in Release mode and install it into a staging folder. The install step runs Qt's deployment tool, which copies the Qt and compiler runtime libraries:

   ```bash
   cmake --install build --prefix stage/Euler2D_FVM
   ```

3. Add `LICENSE.txt`, a short `README.txt`, and a `third-party/` folder with the license texts of Qt (LGPL v3) and of the compiler runtime, then zip the folder and attach it to the GitHub release.
