Strata is a C++ library for computing the dyadic multilayer Green's function (MGF), intended for use in computational electromagnetics codes for solving Maxwell's equations in integral form.
Strata can be used as a standalone library to compute the MGF as a function of frequency and/or spatial separation.
It can also be incorporated easily into existing integral equation solvers.

Source code, documentation and examples will be released in December 2021.

## Prerequisites

* [CMake](https://cmake.org/)

## Installation and Usage

Detailed instructions coming in December 2021.

## Features

The functionality of Strata can be divided into the computation of the MGF in (a) the spectral domain, and (b) the spatial domain.

### a) Computation of the spectral-domain MGF

Strata can compute, in spectral domain:

* Each unique component of the dyadic $`\mathbf{G}^{(\mathrm{A})}`$, and of $`{G}^{(\phi)}`$, as defined in formulation C of [[1]](#mgf01)
* Each unique component of the dyadics $`\mathbf{G}^{(\mathrm{EM})}`$ and $`\mathbf{G}^{(\mathrm{HJ})}`$, as defined in [[2]](#mgf02)

In addition, one can optionally extract the quasistatic terms as defined in [[3]](#qse01).

The code is written with modularity and readability as a priority, and can be extended to incorporate other formulations of the spectral-domain MGF.

### b) Computation of the spatial-domain MGF

Strata allows computing the spatial domain MGF via the following methods:

* Direct numerical integration of the Sommerfeld integrals, as described in [[4]](#mgf03)
* The discrete complex image method (DCIM) [[5]](#dcim01)
* Precomputation of interpolation tables constructed using either one of the above methods

In addition, one can optionally extract (in spectral domain) and add (in spatial domain) the quasistatic terms as defined in [[3]](#qse01).

To easily incorporate into existing codes, one can also extract out the singular behaviour of the MGF, so that the singularities can be treated separately.

The code is written with modularity and readability as a priority, to allow easily extending its functionality.

## Citing Strata

We request that you acknowledge the authors of Strata by citing the following:

```latex
@INPROCEEDINGS{strata,
	author={S. {Sharma} and P. {Triverio}},
	booktitle={2021 {IEEE} International Symposium on Antennas and Propagation and {USNC-URSI} Radio Science Meeting},
	title={Strata: An Open-Source {C++} Library for Computing {Green's} Functions for Layered Media},
	year={2021},
	month={Dec.},
	address = {Singapore}}
```

## References

<a name="mgf01"></a>[1] K. A. Michalski and D. Zheng, "Electromagnetic scattering and radiation by surfaces of arbitrary shape in layered media. I. theory," *IEEE Trans. Antennas Propag.*, vol. 38, no. 3, pp. 335–344, Mar. 1990.

<a name="mgf02"></a>[2] K. A. Michalski and J. R. Mosig, "Multilayered media Green's functions in integral equation formulations," *IEEE Trans. Antennas Propag.*, vol. 45, no. 3, pp. 508–519, Mar. 1997.

<a name="qse01"></a>[3] E. Simsek, Q. H. Liu, and B. Wei, "Singularity Subtraction for Evaluation of Green's Functions for Multilayer Media," *IEEE Trans. Microw. Theory Tech.*, vol. 54, pp. 216–225, Jan 2006.

<a name="mgf03"></a>[4] K. A. Michalski and J. R. Mosig, "Efficient computation of Sommerfeld integral tails - methods and algorithms," *J. Electromagn. Waves Appl.*, vol. 30, no. 3, pp. 281–317, 2016.

<a name="dcim01"></a>[5] M. I. Aksun, "A Robust Approach for the derivation of closed-form Green's functions," *IEEE Trans. Microw. Theory Tech.*, vol. 44, pp. 651–658, May 1996.

## Contributors

Strata was developed as part of a PhD project at the [Modelics Lab](http://modelics.org/), University of Toronto, Canada.

* Shashwat Sharma (Architect & Developer)
* Piero Triverio (Principal Investigator)

