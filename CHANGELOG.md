Changes to Strata are listed in this file.

## [1.0.2] - Unreleased - 2022-07-?

- Added a check when computing the curl of the MGF to make sure `MGF_settings::compute_curl` was set to `true`, so that all curl-related initializations are done (in response to issue #1 opened by @Riarrieta).
- Fixed a mistake in README.md for installing Strata with a pre-existing OpenBLAS build (in response to issue #2 opened by @pedrohnv).

## [1.0.1] - 2022-03-13

- Added another example with digitized data from Yuan, Sarkar, Salazar-Palma, TMTT, 2006; modified `testMGF.cpp` accordingly.
- Added `testSpectralMGF.cpp` to demonstrate how Strata can be used to compute the MGF in spectral (spatial frequency) domain.
- Modified the default Sommerfeld integration path deformation height above the real axis based on experimentation.
- Modified README to include a section on related projects by other authors.

## [1.0.0] - 2021-11-08

- Initial release of Strata
