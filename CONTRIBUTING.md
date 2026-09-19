# Contributing

Contributions are welcome, especially evaluated data improvements, independent test cases, additional decay paths, and analysis notebooks.

## Development setup

Install SciNumTools3 first, then configure this project with its prefix:

```bash
cmake -S . -B build -G Ninja -DCMAKE_PREFIX_PATH=/path/to/snt-prefix
cmake --build build
./build/nuclide-atlas --output build/inventory.csv
python3 python/plot_inventory.py build/inventory.csv build/inventory.png
```

The project requires CMake 3.22, C++20, an SNT installation that includes the host `snt` executable, and Python plus matplotlib only for plotting.

## Data changes

Data is code. Every new or changed numerical record should include:

- a primary evaluated source and access date in the pull request;
- units for every physical quantity;
- isotope identity, atomic mass, stable flag, half-life, daughter, branching, mode, and energy;
- enough significant figures to justify the intended use;
- a scenario or regression test that exercises the new chain.

Do not silently replace evaluated values with web snippets. For a record with multiple meaningful branches, open an issue or extend the schema and solver together; do not hide the branches by picking one arbitrarily.

## Pull requests

Keep PRs focused. Include the command you ran, state whether data are rounded, and call out any model limitation that affects interpretation. We value auditability over a deceptively broad database.
