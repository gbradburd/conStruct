
## Summary

This is a minor update from 1.0.6 to 1.0.7. It fixes a bug in
`calc.covariance` (it used the sample rather than population covariance,
which could produce non-positive-definite covariance matrices in rare
small-sample cases) and allows user-specified pdf sizes in
`compare.two.runs`. See NEWS.md for details.

## Test environments
* local macOS (aarch64-apple-darwin20), R 4.5.2
* win-builder

## R CMD check results
0 errors | 1 warning | 1 note (local check)

* WARNING: "unknown warning group '-Wfixed-enum-extension', ignored"
  from R's own R_ext/Boolean.h header, emitted by the local machine's
  (pre-release) Apple clang 21 toolchain. Not related to conStruct code
  and not expected to reproduce on CRAN's build machines.
* NOTE: "unable to verify current time" -- a transient failure to reach
  the time-check service from the local machine, unrelated to the package.

The GNU-make and installed-library-size NOTEs seen in earlier CRAN
submissions (GNU make is required because the package depends on rstan/
rstantools and is noted in SystemRequirements; the installed size is large
because the C++ Stan models are compiled at installation) did not appear
in this check.

## Downstream dependencies

* There are currently no downstream dependencies for this package.