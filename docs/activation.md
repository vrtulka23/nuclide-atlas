# One-group activation workflow

This workflow models a sample exposed to a constant, spectrum-averaged neutron flux and then allowed to cool. It is intended for transparent exploration and software prototyping, not reactor analysis or safety work.

## Run it

With the default build policy (`dip/build.dip` has `enable_activation = true`):

```bash
./setup.sh -arp -s dip/scenarios/u238_activation.dip
```

`-a` selects the activation CLI path, `-r` writes `build/inventory.csv`, and `-p` produces `build/inventory.png` using `.venv`. The equivalent direct command is:

```bash
./build/nuclide-atlas --activation \
  --scenario dip/scenarios/u238_activation.dip \
  --output build/activation.csv
```

To retain this dashboard alongside the ordinary decay dashboard, use distinct names:

```bash
./build/nuclide-atlas --activation \
  --scenario dip/scenarios/u238_activation.dip \
  --output build/u238-activation.csv
.venv/bin/python python/plot_inventory.py \
  build/u238-activation.csv build/u238-activation.png \
  --scenario dip/scenarios/u238_activation.dip --data-dir dip/data
```

If `build.enable_activation` is false, `--activation` is rejected and the matrix solver is not compiled. Edit `dip/build.dip` and rerun `./setup.sh -b` to change the build capability.

## Input contract

Activation scenarios are DIPL instances of `activation_simulation`, declared in `dip/data/nuclear.dip`:

```dip
activation : activation_simulation
  title = "U-238 activation to plutonium-239"
  sample.isotope = "U238"
  sample.mass = 1 g
  irradiation.duration = 30 day
  irradiation.points = 301
  irradiation.neutron_flux = 1e14 1/(cm2*s)
  irradiation.capture.target = "U238"
  irradiation.capture.product = "U239"
  irradiation.capture.cross_section = 2.68e-24 cm2
  cooldown.duration = 2 yr
  cooldown.points = 401
  output.csv = "activation.csv"
  output.include_activity = true
  output.include_q_power = true
```

DIPL validates identifiers, positive sample mass, point limits, nonnegative cooldown, and physical dimensions. PUQ converts flux, cross section, and durations before C++ receives them. The microscopic capture rate is `flux × cross_section`, in `s-1`.

## Model and numerical method

The solver assembles a coupled transition matrix. For every radioactive parent `j`, it adds a diagonal loss `-lambda_j` and a feed to the declared daughter `b_ji lambda_j`. During irradiation it also adds capture loss from `target` and equal feed to `product`:

```text
dN/dt = A N
A[target,target] -= flux × cross_section
A[product,target] += flux × cross_section
```

The same natural-decay network is used in cooldown, without the capture terms. Any undeclared fraction of a radioactive branch feeds an absorbing, unreported loss state. This makes the propagated system closed: the solver normalizes the stochastic transition columns during scaling-and-squaring and checks positivity and total atom balance after every propagation. Reporting points control output resolution only.

The bundled scenario joins the ordinary U-238 decay series to the capture path U-238 → U-239 → Np-239 → Pu-239. It uses an illustrative thermal capture cross section. The flux and cross section must be averaged over the same neutron spectrum; substituting either value independently is not physically meaningful.

## Output contract

Activation CSVs always include `time_s`, `time_yr`, `phase`, `total_atoms`, `capture_reactions_per_s`, `untracked_loss_atoms`, and each `<isotope>_atoms` column. `phase` is `irradiation` or `cooldown`; the capture reaction rate is zero during cooldown. `total_atoms + untracked_loss_atoms` is atom-balanced to the initial inventory (within floating-point tolerance); the latter makes omitted branches visible instead of silently changing the balance.

`output.include_activity` adds `total_activity_bq` and `<isotope>_activity_bq` columns. `output.include_q_power` adds `total_q_power_w`. The latter is a decay Q-value rate obtained from activity times the tabulated Q value and `nuclear.constants.joules_per_mev`; it is not deposited heat.

The bundled plot is a dashboard: a log-inventory heatmap shows which members dominate in each period; a companion log-log trajectory panel retains the individual nuclide curves; a dual-axis panel compares activity with Q-value power; capture production is shown only while flux is applied; and a final-composition bar chart ranks the remaining inventory. The dashed line marks the irradiation-to-cooldown transition.

## Limits

- One constant, spectrum-averaged flux and one capture reaction per scenario.
- No fission, depletion through competing reactions, self-shielding, resonance treatment, multigroup transport, spatial effects, or changing flux.
- One explicit radioactive daughter per isotope record; omitted decay branches are lost from the reported inventory.
- No uncertainty propagation, radiation transport, dose calculation, or thermal model.

Use evaluated, spectrum-appropriate reaction data for any serious investigation. The bundled records are compact demonstrator data with [NUBASE2020](https://www-nds.iaea.org/amdc/ame2020/NUBASE2020.pdf) provenance; NNDC’s [NuDat documentation](https://www.nndc.bnl.gov/nudat3/guide/) describes its evaluated reaction and decay data coverage.

## One-group deuteron reactions

`--deuteron-activation` accepts a separate `deuteron_activation_simulation` contract. It propagates one target-residual reaction edge while a constant deuteron flux is applied, then follows ordinary radioactive decay during cooldown:

```dip
deuteron_activation : deuteron_activation_simulation
  title = "Deuterium target: H-2(d,p)H-3 demonstrator"
  sample.isotope = "H2"
  sample.mass = 1 ug
  irradiation.duration = 1 day
  irradiation.points = 101
  irradiation.deuteron_flux = 1e12 1/(cm2*s)
  irradiation.reaction.emitted_particle = "p"
  irradiation.reaction.target = "H2"
  irradiation.reaction.product = "H3"
  irradiation.reaction.cross_section = 1e-27 cm2
  cooldown.duration = 30 yr
  cooldown.points = 301
  output.csv = "h2-deuteron-activation.csv"
  output.include_activity = true
  output.include_q_power = true
```

The supported labels are `d,p`, `d,n`, and `d,alpha`; each denotes the emitted particle in `target(d,particle)product` notation. The target-to-product edge is a residual-inventory model, not a transport simulation of emitted particles. Its rate is `deuteron_flux × cross_section`, and the output column is `deuteron_reactions_per_s`.

The model is appropriate for transparent, effective-rate demonstrations when the supplied cross section has already been averaged over the actual incident-energy distribution and target depth. The dashboard adds a cumulative-yield panel by trapezoidally integrating `deuteron_reactions_per_s`, alongside the instantaneous reaction-rate panel. It does not calculate the energy dependence of a deuteron reaction, beam slowing or straggling, charged-particle transport, angular distributions, heating, or several competing channels. The included H-2(d,p)H-3 cross section is illustrative only, not evaluated data.
