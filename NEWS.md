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
