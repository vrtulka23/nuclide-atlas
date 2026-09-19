# Architecture

```text
scenario.dip ──┐
              ├─ DIPL parser / validation ── resolved Environment ── C++ Bateman solver ────── CSV
isotope *.dip ─┤                                                                            │
constants.dip ─┘                                                                            └─ Python plot
                                      ▲
dip/build.dip ── snt_dip_get() ──────┘ (build-time policy)

activation scenario.dip ── DIPL activation_simulation ── coupled rate matrix ── exp(A·t) ── CSV
```

`dip/data/nuclear.dip` is a registry, not a monolithic data dump. Its `$source` records point at individual files and import their `isotope` tables as keyed items beneath `nuclear.isotopes[<id>]`. This makes diffs and provenance reviews local to one record.

At run time `src/main.cpp` reads quantities through a DIPL cursor, reconstructs a PUQ quantity from its numeric value and unit metadata, and explicitly converts it to the solver’s base unit. As a result, a data mistake such as assigning a mass to `half_life` is caught as a dimensional incompatibility before it can become a numerical result.

The solver makes decay constants from `ln(2)/half_life` and evaluates the serial-chain Bateman solution. Each snapshot is exact for the supplied one-daughter chain; the reporting grid therefore has no effect on stiffness or stability.

At configure time, SNT's `snt_dip_get()` evaluates `dip/build.dip`. This is a real CMake/DIPL integration: build policy is not duplicated in a second config syntax.

Besides the C++ standard and feature switches, `dip/build.dip` selects the default scenario. CMake embeds that DIPL-selected path in the executable, while `--scenario` remains an explicit runtime override. Physical constants and CSV inclusion policy stay in the runtime DIPL data/scenario layer rather than being duplicated in CMake or C++.

When `build.enable_activation` is true, CMake emits `SNT_NUCLEAR_ENABLE_ACTIVATION` for the executable. The activation code is therefore absent, rather than merely hidden at runtime, when the DIPL switch is false. Activation scenarios are separately validated against `activation_simulation`: PUQ converts flux in `1/(cm2*s)`, microscopic cross sections in `cm2`, and both phase durations to SI before their product becomes a reaction rate in `s-1`.

The activation engine constructs a coupled transition matrix from each isotope’s validated decay daughter and the scenario’s capture target/product. It adds an absorbing loss state for omitted radioactive branches, then evaluates matrix exponentials by scaling-and-squaring with a converged Taylor series. Transition columns are normalized through each squaring and every propagated state is checked for positivity and atom balance. This replaces grid-dependent stepping for irradiation and cooldown while retaining the exact Bateman path for ordinary serial-decay scenarios; cancellation-level Bateman remnants are reported as zero.
