# MissLearn

Teaching application for the four-step framework for identifying and handling
missing data in longitudinal cohorts.

Part of the **Longitudinal Methods Platform**, Institute of Psychiatry,
Psychology and Neuroscience, King's College London.

All data in this application is simulated. No patient data is used, and nothing
is stored.

## What it does

Simulates a longitudinal cohort with a known group effect and a known
missing-data mechanism, then works through four steps:

1. **Describe** — how much is missing, at which visits, in what pattern.
2. **Diagnose** — model the probability of dropping out on observed values.
3. **Handle** — compare complete-case analysis and multiple imputation against
   a truth only a simulation can supply.
4. **Sensitivity** — shift imputed values under a departure from MAR and find
   the point at which the conclusion changes.

Imputation is by fully conditional specification with predictive mean matching
(`mice`), pooled with Rubin's rules.

## Run locally

```r
install.packages(c("shiny", "mice"))
shiny::runApp()
```

## Deploy to Posit Connect Cloud

Connect Cloud needs a public GitHub repository containing `app.R` and a
`manifest.json`. Generate the manifest from this directory:

```r
install.packages("rsconnect")
rsconnect::writeManifest()
```

Commit `manifest.json` alongside `app.R`, push, then publish from Connect Cloud
and select this repository. Regenerate the manifest whenever a package is added,
or the build will fail on a dependency that exists only on your machine.

## Licence

MIT. See `LICENSE`.
