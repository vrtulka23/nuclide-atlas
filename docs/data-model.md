# Data model and provenance

Each isotope source has one `isotope` table:

```dip
isotope
  id str = "U238"
  atomic_mass float = 238.05078826 g/mol
  stable bool = false
  half_life float = 4.468e9 yr
  decay
    daughter str = "Th234"
    branching float = 1.0
    mode str = "alpha"
    energy float = 4.26975 MeV
```

All numerical fields have a unit. `data/nuclear.dip` imports a file per record and is the one explicit registry that determines which records are visible to a scenario.

## Isotope-record schema

The registry declares an `isotope_record` schema and instantiates it for every entry in `nuclear.isotopes`. It then injects the scalar values from each source record into that schema-backed entry. This is intentional: the schema validates the **resolved database** consumed by the C++ solver, while the small source files stay independently reviewable and retain their value-level provenance.

The shared contract checks isotope and daughter identifier formats, labels, positive molar mass, non-negative half-life and decay energy, a branching fraction in `[0, 1]`, and the supported serial-chain modes (`alpha`, `beta-`, `stable`). A malformed imported record therefore fails during DIPL parsing, before a calculation can begin.

## Simulation schema

`data/nuclear.dip` also declares the reusable `nuclear_simulation` schema. It defines the input contract once: a bounded title, sample-isotope format, mass in mass dimensions, duration and minimum time in time dimensions, a `minimum_time < duration` cross-field rule, a 2-to-1,000,000 reporting-point limit, supported time grids, and a CSV filename format. A scenario is therefore only an assignment layer:

```dip
simulation : nuclear_simulation
  sample.isotope = "U238"
  sample.mass = 1 g
  duration = 4.5 Gyr
  points = 1000
```

This is more than a template: DIPL applies the schema while parsing and rejects incompatible units or invalid choices before the C++ solver sees them.

## DIPL features used without custom infrastructure

The project gets the following capabilities from DIPL/PUQ rather than reimplementing them around a generic text format:

- `$source` and hierarchical imports compose one auditable database from individual isotope files.
- `$schema` applies one contract to every resolved isotope record and another to every scenario.
- `!format`, `!options`, and `!condition` validate strings, enumerations, units, ranges, and cross-field relationships.
- PUQ parses and converts physical quantities, so `4.5 Gyr` satisfies a schema declared in seconds and `100 ng` satisfies one declared in grams.
- `?authors`, `?doi`, `?url`, and related properties attach citations to a value without adding solver data fields.
- Expressions, references, arrays, custom units, tags, constants, and DIPL-to-CMake `snt_dip_get()` remain available as the model grows.

Provenance is attached to the `id` value through DIPL metadata properties—not model data nodes—using `?authors`, `?title`, `?journal`, `?year`, `?doi`, `?url`, and `?version`. That keeps citations available to DIPL-aware consumers while ensuring they never enter the solver's parameter namespace.

The initial seed is compiled from rounded commonly published evaluated values to keep DIPL files readable. It is a demonstrator dataset, not a claimed replacement for ENSDF/NUBASE. A production contribution should attach the precise evaluated source, evaluator/version, retrieval date, and any branch-selection rationale using these metadata properties.

The current solver follows one daughter branch per isotope record. The `branching` value preserves loss to omitted channels, but does not model their inventories. Full branching support is a natural next schema-and-solver extension.
