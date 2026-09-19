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

All numerical fields have a unit. `dip/data/nuclear.dip` imports a file per record and is the one explicit registry that determines which records are visible to a scenario.

## Isotope-record schema

The registry declares `nuclear.isotopes` as a keyed `map : isotope_record`. DIPL automatically applies the schema to every item, while each source is imported whole into its corresponding `nuclear.isotopes[<id>]` entry. This validates the **resolved database** consumed by the C++ solver without repeating field-by-field injections; the small source files remain independently reviewable and retain their value-level provenance.

The shared contract checks isotope and daughter identifier formats, labels, positive molar mass, non-negative half-life and decay energy, a branching fraction in `[0, 1]`, and the supported serial-chain modes (`alpha`, `beta-`, `stable`). A malformed imported record therefore fails during DIPL parsing, before a calculation can begin.

## Simulation schema

`dip/data/nuclear.dip` also declares the reusable `nuclear_simulation` schema. It defines the input contract once: a bounded title, sample-isotope format, mass in mass dimensions, duration and minimum time in time dimensions, a `minimum_time < duration` cross-field rule, a 2-to-1,000,000 reporting-point limit, supported time grids, CSV filename format, and flags for activity and Q-value-power columns. A scenario is therefore only an assignment layer:

```dip
simulation : nuclear_simulation
  sample.isotope = "U238"
  sample.mass = 1 g
  duration = 4.5 Gyr
  points = 1000
```

This is more than a template: DIPL applies the schema while parsing and rejects incompatible units or invalid choices before the C++ solver sees them.

`constants.dip` is equally part of the runtime contract: Avogadro’s constant sets the atom count, `seconds_per_year` sets CSV year conversion, and `joules_per_mev` converts decay Q values to joules. The C++ solver reads and unit-converts these values through PUQ; it does not maintain duplicate numerical literals.

## Activation-simulation schema

`activation_simulation` adds a separately typed request for a one-group irradiation followed by cooldown. It validates the sample isotope/mass, irradiation duration and reporting points, neutron flux in `1/(cm2*s)`, capture target/product identifiers, microscopic cross section in `cm2`, cooldown duration/points, and a CSV output name. The U-238 example is concise because units and bounds live centrally:

```dip
activation : activation_simulation
  sample.isotope = "U238"
  sample.mass = 1 g
  irradiation.neutron_flux = 1e14 1/(cm2*s)
  irradiation.capture.target = "U238"
  irradiation.capture.product = "U239"
  irradiation.capture.cross_section = 2.68e-24 cm2
```

This showcases the same SNT contract used by the decay model: a malformed unit or an invalid rate input is rejected by DIPL/PUQ before the C++ rate matrix is assembled. The cross section is intentionally scenario data, so its spectrum dependence is explicit and reviewable rather than hidden in solver code.

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

Each isotope record still supports one explicit radioactive daughter branch. The activation engine can add one scenario-defined neutron-capture feed edge to that decay network. The `branching` value preserves loss to omitted radioactive channels in an explicit absorbing balance state (`untracked_loss_atoms`), but does not model their inventories. Full decay branching and multiple reaction channels are natural next schema-and-solver extensions.
