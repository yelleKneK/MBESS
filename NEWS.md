# MBESS 5.1.0 (2026)

Maintenance release. Corrects the seed handling in the cluster
randomized design AIPE sample size planners.

* The `ss.aipe.crd.es.*` planners (`ss.aipe.crd.es.nclus.fixedwidth`,
  `ss.aipe.crd.es.nindiv.fixedwidth`, `ss.aipe.crd.es.nclus.fixedbudget`,
  `ss.aipe.crd.es.nindiv.fixedbudget`, `ss.aipe.crd.es.both.fixedbudget`,
  and `ss.aipe.crd.es.both.fixedwidth`) now default `seed` to `NULL`
  instead of the fixed value `123321`. With the previous default, every
  call silently reused the same seed, which made the a priori Monte Carlo
  results look perfectly reproducible and misled users about the sampling
  variability of the planned sample size. With `seed = NULL` the
  simulation now draws from the current state of the random number
  generator (honest variability), and the planner restores the caller's
  RNG state on exit so that calling it does not perturb the user's random
  number stream. Supplying an explicit integer `seed` still produces
  reproducible results exactly as before.

# MBESS 5.0.1 (2026)

Maintenance release. Corrects the sampling-variance scaling in the
polynomial change model sample size planner `ss.aipe.pcm()`.

* `ss.aipe.pcm()`: the single-subject sampling variance of the change
  coefficient now includes the factor `frequency^(2m)` (where `m` is the
  order of the polynomial trend), matching Raudenbush and Liu (2001) and
  the companion power function `ss.power.pcm()`. Previously this factor
  was omitted, so for designs measured more than once per unit of time
  (`frequency` greater than 1) the planned sample size was too small and
  the resulting confidence interval was wider than requested. Results
  for `frequency = 1` are unchanged, including the published Kelley and
  Rausch (2011) examples, which still return 278 and 316.
* `ss.aipe.pcm()`: corrected the constant for a cubic trend from
  `1/1000800` to the correct `1/100800` (that is, `(3!)^2 / (6! 7!)`).
  The earlier value, a transcription error with an extra digit, inflated
  the implied error variance for `trend = "cubic"` by a factor of about
  9.93 and so overstated the planned sample size. Linear and quadratic
  trends were unaffected.
* `ss.aipe.pcm()`: in the assurance branch, the t critical value is now
  recomputed as the sample size iterates rather than reusing the value
  from the expected-width search. The numerical effect is small, but the
  previous code used a critical value carrying the wrong degrees of
  freedom.

# MBESS 5.0.0 (2026)

Maintenance release. MBESS is mature, stable, and widely used, and it
remains fully supported. Maintenance mode does not mean the package is
finished or abandoned: MBESS will continue to be maintained so that it
keeps working with current versions of R and its dependencies, and it
will remain on CRAN. New methodological development now takes place in
the successor package, DMAR (Design, Measurement, and Analysis in R),
which is available on CRAN. MBESS and DMAR are designed to coexist;
existing code that relies on MBESS will keep working, so there is no
need to migrate.

* Documentation and references brought up to date. Citations that had
  previously appeared as "submitted", "in press", or "under review" now
  cite the published articles (for example, Kelley, 2008; Kelley & Lai,
  2011; Lai & Kelley, 2012; Chattopadhyay & Kelley, 2016, 2017;
  Lachowicz, Preacher, & Kelley, 2018; Li, McLouth, & Delaney, 2020).
* Numerous spelling and typographical corrections throughout the help
  pages.
* `ci.reliability()`: fixed a misnamed variable that caused the
  lower/upper logit confidence limit swap to be skipped in an edge case.
* Internal calls into Suggested packages were brought up to their
  current, fully namespaced forms. `ss.aipe.rmsea.sensitivity()` now
  calls `sem::sem()` directly, so it no longer fails when the `sem`
  package is installed but not attached. `ci.reliability()` now uses
  `lavaan::lavInspect()` in place of the older `lavaan::inspect()`,
  with identical numerical results.
* Added an `inst/CITATION` file, this `NEWS` file, and an expanded
  `README`.
* Contributor acknowledgments now ship with the package (see the
  `CONTRIBUTORS` file in the installed package) rather than depending on
  an external web page.

# MBESS 4.0.0 (2016)

Major release.

# MBESS 3.0.0 (2010)

Major release. Versions 3.0.0 through 3.3.3 included Keke (Brian) Lai as
coauthor and maintainer. Lai, Ken Kelley's doctoral advisee at Indiana
University and the University of Notre Dame, has since stepped away from
package maintenance and development, but his contributions to individual
functions remain and he is credited as their coauthor.

# MBESS 2.0.0 (2008)

Major release.

# MBESS 1.0.0 (2007)

Major release.

# MBESS 0.0.1 (2006-05-06)

First public release.
