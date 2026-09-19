# Architecture

```text
scenario.dip ──┐
              ├─ DIPL parser / validation ── resolved Environment ── C++ Bateman solver ────── CSV
isotope *.dip ─┤                                                                            │
constants.dip ─┘                                                                            └─ Python plot
                                      ▲
cmake/build.dip ── snt_dip_get() ────┘ (build-time policy)
```

`data/nuclear.dip` is a registry, not a monolithic data dump. Its `$source` records point at individual files and import their `isotope` tables as keyed items beneath `nuclear.isotopes[<id>]`. This makes diffs and provenance reviews local to one record.

At run time `src/main.cpp` reads quantities through a DIPL cursor, reconstructs a PUQ quantity from its numeric value and unit metadata, and explicitly converts it to the solver’s base unit. As a result, a data mistake such as assigning a mass to `half_life` is caught as a dimensional incompatibility before it can become a numerical result.

The solver makes decay constants from `ln(2)/half_life` and evaluates the serial-chain Bateman solution. Each snapshot is exact for the supplied one-daughter chain; the reporting grid therefore has no effect on stiffness or stability.

At configure time, SNT's `snt_dip_get()` evaluates `cmake/build.dip`. This is a real CMake/DIPL integration: build policy is not duplicated in a second config syntax.
