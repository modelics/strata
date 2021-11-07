.. strata documentation main file.
.. Author: Shashwat Sharma
.. Created on: Nov. 06, 2021

Introduction
============

Strata is a C++ library for computing the multilayer Green's function (MGF) for computational electromagnetics codes based on the boundary element method. The source code can be found on `GitHub <https://github.com/modelics/strata>`_. 

If you are new to Strata, we recommend that you start from the :ref:`basicusage` section. A series of examples is provided in ``strata/test/`` for a hands-on introduction to various aspects of Strata, and we strongly recommend using these examples as a starting point in your application.

Installation instructions can be found in the :ref:`build` section.

Table of Contents
-----------------

.. toctree::
   :maxdepth: 2

   build
   basic_usage
   settings_list
   tech_example
   people

.. _cite:

Citing Strata
-------------

We request that you acknowledge the authors of Strata by citing the following:

.. code-block:: latex

	@INPROCEEDINGS{strata,
		author={S. {Sharma} and P. {Triverio}},
		booktitle={2021 {IEEE} International Symposium on Antennas and Propagation and {USNC-URSI} Radio Science Meeting},
		title={Strata: An Open-Source {C++} Library for Computing {Green's} Functions for Layered Media},
		year={2021},
		month={Dec.},
		address = {Singapore}}

Acknowledgment
--------------

* Strata ships with files from the libAmosBessel library by Homer Reid, which provides C++ wrappers for computing Bessel functions with the original fortran routines by T. E. Amos.
* Also included are selected files from the `Boost <https://www.boost.org/>`_  C++ library for numerical integration.


.. Indices and tables
.. ==================

.. * `GitHub page <https://github.com/modelics/strata>`_
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
