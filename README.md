DeadWood_DWDDecay
================
2026-05-15

## Overview

`DeadWood_DWDDecay` manages the downed woody debris (DWD) pool. Each
5-year timestep it receives fallen snags from `DeadWood_snagDecay`,
stochastically maps their snag decay class (DC1–DC5) to a starting DWD
decay class (DC1–DC4), then advances all DWD pieces through DC1–DC4+
using a time- and diameter-dependent logistic transition model. Pieces
exiting DC4 are removed as fully decomposed.

## Inputs

| Object | Class | Description |
|----|----|----|
| `fallenSnags` | `data.table` | Snags that fell in the current 5-year timestep, produced by `DeadWood_snagDecay`. Required columns: `pixelID` (integer), `species` (character), `DC` (integer, snag DC 1–5), `initBiomass` (numeric, Mg ha⁻¹), `diameter_cm` (numeric, cm). |

**Minimum diameter:** Any piece with `diameter_cm < 7.5 cm` will cause
an error.

## Outputs

| Object | Class | Description |
|----|----|----|
| `DWDTable` | `data.table` | Current DWD inventory, updated every 5 years. Fully decomposed pieces (exited DC4) are removed. |

`DWDTable` schema:

| Column | Type | Description |
|----|----|----|
| `pixelID` | integer | Raster pixel of origin |
| `species` | character | Tree species |
| `DC` | integer | Current DWD decay class (1–4) |
| `ageInDC` | integer | Years spent in the current DWD decay class |
| `ageSinceEntry` | integer | Total years since the piece entered the DWD pool |
| `initBiomass` | numeric | Biomass at time of death (Mg ha⁻¹); never modified after entry |
| `diameter_cm` | numeric | Piece diameter (cm); never modified after entry |

**Note on the DWD decay class scale:** The DWD pool uses a 4-class
scale. DWD DC4 combines what would be snag DC4 and DC5 — both represent
highly decomposed material on the forest floor.

**`ageInDC` vs `ageSinceEntry`:** These two columns track time
differently. `ageInDC` counts years spent in the *current* decay class
and resets to zero every time a piece advances to a new class — it
measures how long a piece has been at its present level of decay.
`ageSinceEntry` counts total years since the piece first entered the DWD
pool and never resets — it is the age variable (*A*) used in the
logistic transition probability calculation.

## Parameters

| Parameter            | Type         | Default source                 |
|----------------------|--------------|--------------------------------|
| `DWD_logisticParams` | `data.table` | Vanderwel et al. 2006          |
| `snagToDWD_DCmat`    | `matrix`     | Vanderwel et al. 2006, Table 4 |

### `DWD_logisticParams`

A 4-row `data.table` with columns `a0`, `a1`, `b0`, `b1` controlling the
transition probability for each DC. One row per DC transition:

| Row | Transition                 |
|-----|----------------------------|
| 1   | DC1 → DC2                  |
| 2   | DC2 → DC3                  |
| 3   | DC3 → DC4                  |
| 4   | DC4 → decomposed (removed) |

Parameters `a` and `b` are piece-specific and scale with diameter (*D*,
cm):

$$a = a_0 + a_1 \ln(D), \qquad b = b_0 + b_1 \ln(D)$$

These feed into a logistic survival function used to compute conditional
transition probabilities — see the Events section for the full
derivation.

**Note:** The transition calculation assumes a fixed 5-year interval
matching the module’s timestep. This interval is hardcoded in
`computeDWDTransProb()` (`R/dwd-transition.R`). If the module timestep
were changed, this function would also need to be updated.

Default fitted coefficients (Vanderwel et al. 2006):

| Transition |     a0 |     a1 |     b0 |     b1 |
|------------|-------:|-------:|-------:|-------:|
| DC1 → DC2  | −7.295 |  1.888 |  1.325 | −0.239 |
| DC2 → DC3  | −9.663 |  0.000 |  0.934 |  0.000 |
| DC3 → DC4  |  3.882 | −2.986 | −0.076 |  0.093 |
| DC4 → out  | −5.584 |  0.000 |  0.149 |  0.000 |

When `a1 = 0` and `b1 = 0` the transition rate is diameter-independent.

### `snagToDWD_DCmat`

A 5×4 probability matrix that stochastically maps the snag decay class
at the moment of fall (rows 1–5) to the starting DWD decay class
(columns 1–4). DWD is measured on a 4-class scale; snag DC4 and DC5 both
map probabilistically into DWD DC4.

Default values (Vanderwel et al. 2006, Table 4); rows sum to 1:

|              | DWD DC1 | DWD DC2 | DWD DC3 | DWD DC4 |
|--------------|--------:|--------:|--------:|--------:|
| **Snag DC1** |   0.000 |   0.162 |   0.815 |   0.023 |
| **Snag DC2** |   0.000 |   0.078 |   0.871 |   0.052 |
| **Snag DC3** |   0.000 |   0.036 |   0.854 |   0.111 |
| **Snag DC4** |   0.000 |   0.016 |   0.762 |   0.222 |
| **Snag DC5** |   0.000 |   0.007 |   0.599 |   0.395 |

Most fallen snags enter DWD at DC3 regardless of their snag DC; the
probability of entering a higher DWD DC increases as snag DC increases,
reflecting more advanced pre-fall decomposition.

## Events

The module fires three event types.

### `init` (once, at simulation start)

- Validates that `DWD_logisticParams` is a 4-row `data.table` and
  `snagToDWD_DCmat` is a 5×4 matrix.
- Creates an empty `DWDTable` with the correct schema.
- Schedules `receive` (priority 2) and `transition` (priority 3) at
  `start(sim) + 5`.

### `receive` (every 5 years, priority 2)

Runs before `transition` within the same timestep. For each piece in
`sim$fallenSnags`, samples a DWD entry DC from the row of
`snagToDWD_DCmat` corresponding to the piece’s snag DC:

$$\text{DWD-DC}_\text{entry} \sim \text{Categorical}\!\left(\text{snagToDWD\_DCmat}[\text{snag-DC},\, \cdot\,]\right)$$

Sets `ageInDC = 0` and `ageSinceEntry = 0` for all incoming pieces, then
appends them to `DWDTable`. Does nothing if `fallenSnags` is empty or
`NULL`.

### `transition` (every 5 years, priority 3)

Advances all pieces currently in `DWDTable` through the logistic decay
model.

**Step 1 — Compute piece-specific logistic parameters**

For each piece, `a` and `b` are derived from diameter using the
parameter row for the piece’s current DC:

$$a = a_0 + a_1 \cdot \ln(D), \qquad b = b_0 + b_1 \cdot \ln(D)$$

where $D$ is `diameter_cm` and $(a_0, a_1, b_0, b_1)$ come from
`DWD_logisticParams` indexed by the piece’s current DC.

**Step 2 — Logistic progression function**

The probability that a piece has *not yet* transitioned to the next DC
by age $t$ (years since DWD entry) is:

$$Pr = \frac{1}{1 + e^{a + b \cdot t}}$$

**Step 3 — Conditional transition probability**

For a piece at age $A$ (`ageSinceEntry`), the probability of
transitioning during the next 5-year interval $[A,\, A+5]$, given it has
not yet transitioned, is:

$$p_\text{transition} = \frac{Pr(A) - Pr(A+5)}{Pr(A)}$$

This is clamped to $[0, 1]$. If $Pr(A) \approx 0$ the piece is set to
transition with certainty.

**Step 4 — Stochastic transition and removal**

Each piece independently draws:

$$\text{Transition} \sim \text{Bernoulli}(p_\text{transition})$$

- Pieces that transition: `DC` increments by 1; `ageInDC` resets to 0.
- Pieces that do not transition: `DC` unchanged; `ageInDC` increments by
  5.
- All pieces: `ageSinceEntry` increments by 5 regardless.
- Pieces with `DC > 4` are removed from `DWDTable` (fully decomposed).

## Event scheduling and module interactions

`DeadWood_Mortality` fires every year; the remaining modules fire every
5 years. Priority determines order only among events scheduled at the
same simulation time.

| Timestep | Priority | Module | Event | Purpose |
|----|----|----|----|----|
| every 1 yr | 0 | `DeadWood_Mortality` | `transition` | Apply baseline mortality, append to `cohortData` |
| every 5 yr | 1 | `DeadWood_snagDecay` | `transition` | Advance snag DC, produce `fallenSnags` |
| every 5 yr | 2 | `DeadWood_DWDDecay` | `receive` | Accept `fallenSnags` into DWD pool ← **this module** |
| every 5 yr | 3 | `DeadWood_DWDDecay` | `transition` | Advance DWD DC via logistic model ← **this module** |
| every 5 yr | 4 | `DeadWood_Biomass` | `transition` | Compute biomass rasters from updated pools |
| every 5 yr | 5 | `DeadWood_Biomass` | `plot` | Capture biomass snapshots |

The `receive` event fires before `transition` within the same timestep,
ensuring newly fallen snags are incorporated before any DC advancement
occurs.

## Example

``` r
library(SpaDES.project)
set.seed(42)

# Fallen snags spanning all snag DCs and a range of piece diameters.
# The logistic transition model is diameter-dependent: fine wood (< 20 cm)
# advances through decay classes faster than coarse wood (> 30 cm).
fallenSnags <- data.table::data.table(
  pixelID     = 1:50,
  species     = rep(c("Pinus strobus", "Pinus resinosa"), 25),
  DC          = rep(1:5, 10),
  ageInDC     = 0L,
  initBiomass = 15.0,
  diameter_cm = rep(c(10.0, 15.0, 25.0, 35.0, 45.0), 10)
)

out <- SpaDES.project::setupProject(
  paths       = list(modulePath = "modules"),
  modules     = "BosunForestEcology/DeadWood_DWDDecay@main",
  times       = list(start = 0, end = 50),
  fallenSnags = fallenSnags
)

mySim <- SpaDES.core::simInitAndSpades2(out)

library(ggplot2)

dcSummary <- mySim$DWDTable[, .(biomass = sum(initBiomass)), by = DC][order(DC)]

ggplot(dcSummary, aes(x = DC, y = biomass)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 1:4) +
  labs(x = "Decay class", y = "Total DWD biomass (Mg/ha)",
       title = "DWD biomass by decay class at year 50") +
  theme_bw()
```

## References

Vanderwel, M.C., Malcolm, J.R., Smith, S.M., and Islam, N. (2006). An
integrated model for snag and downed woody debris decay class
transition. *Forest Ecology and Management*, 234(1–3), 48–59.
<https://doi.org/10.1016/j.foreco.2006.06.020>

## Package dependencies

- [`SpaDES.core`](https://github.com/PredictiveEcology/SpaDES.core) (\>=
  3.0.0)
- [`data.table`](https://CRAN.R-project.org/package=data.table)
