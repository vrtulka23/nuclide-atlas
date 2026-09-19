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

## Simulation schema

`data/nuclear.dip` also declares the reusable `nuclear_simulation` schema. It defines the input contract once: sample isotope format, mass in mass dimensions, duration and minimum time in time dimensions, positive-value constraints, reporting-point minimum, supported time grids, and output fields. A scenario is therefore only an assignment layer:

```dip
simulation : nuclear_simulation
  sample.isotope = "U238"
  sample.mass = 1 g
  duration = 4.5 Gyr
  points = 1000
```

This is more than a template: DIPL applies the schema while parsing and rejects incompatible units or invalid choices before the C++ solver sees them.

Provenance is attached to the `id` value through DIPL metadata properties—not model data nodes—using `?authors`, `?title`, `?journal`, `?year`, `?doi`, `?url`, and `?version`. That keeps citations available to DIPL-aware consumers while ensuring they never enter the solver's parameter namespace.

The initial seed is compiled from rounded commonly published evaluated values to keep DIPL files readable. It is a demonstrator dataset, not a claimed replacement for ENSDF/NUBASE. A production contribution should attach the precise evaluated source, evaluator/version, retrieval date, and any branch-selection rationale using these metadata properties.

The current solver follows one daughter branch per isotope record. The `branching` value preserves loss to omitted channels, but does not model their inventories. Full branching support is a natural next schema-and-solver extension.
