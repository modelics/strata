.. Author: Shashwat Sharma
.. Created on: Nov 06, 2021

.. _basicusage:

Basic Usage
===========

There are three basis phases: defining the layer stackup, initializing the MGF engine, and computing the MGF.

.. _stackupdef:

Stackup Definition
------------------

The geometric extents and material parameters of each layer are stored in and managed by objects of the class ``LayerManager``.
This object must be provided with stackup information, which can be done via a text file in the ``.yaml`` file format.
The formatting of the layer definition file is discussed here: :ref:`egtech`.
Once the layer definition file is created, the ``LayerManager`` object can be initialized as follows:

.. code-block:: C++

    LayerManager lm;
    lm.ProcessTechFile(name_of_tech_file);

The argument ``name_of_tech_file`` is of type ``std::string`` and contains the path to the layer definition file, including the ``.yaml`` extension.

.. _initialization:

Initialization
--------------

Next, the desired settings are provided to the MGF engine and various precomputations are performed internally.
The settings are supplied via the ``MGF_settings`` data structure. A list of basic settings is provided here: :ref:`mgfsettings`.
The operating frequency ``f``, the ``LayerManager`` object, and the ``MGF_settings`` object are provided as arguments to the ``MGF`` class via the ``MGF::Initialize`` function:

.. code-block:: C++

    double f = 30.0e9;

    MGF_settings s;
    s.method = MGF_INTEGRATE;

    MGF mgf;
    mgf.Initialize(f, lm, s);

The object ``mgf`` is now responsible for the actual computation of the MGF.
Various MGF computation methods are available; direct Sommerfeld integration is considered here.
See the examples in ``testMGF.cpp``, ``testInterp.cpp``, and ``testDCIM.cpp`` for details on other methods.
The initialization step must be performed at every frequency, so it is recommended to keep the ``mgf`` object local to the frequency loop, and call ``MGF::Initialize`` at each frequency point.

.. _computation:

Computation
-----------

The MGF can now be computed for given source and test coordinates, (``x_obs``, ``y_obs``, ``z_obs``) and (``x_src``, ``y_src``, ``z_src``), respectively:

.. code-block:: C++

    double x_diff = x_obs - x_src;
    double y_diff = y_obs - y_src;

    std::array<std::complex<double>, 9> G_dyadic;
    std::complex<double> G_phi;
    mgf.ComputeMGF(x_diff, y_diff, z_obs, z_src, G_dyadic, G_phi);

Strata uses the MGF formulation-C proposed by Michalski and Zheng [1]_, which involves a dyadic with 9 components (stored in row-major order in ``G_dyadic``), and one scalar component (stored in ``G_phi``).
Some entries of ``G_dyadic`` are always 0 in formulation-C, but all 9 components are returned for compatibility with other MGF formulations which one may want to implement in the future.
Note that the MGF computed above **does include** the cos and sin pre-factors [1]_.

References
----------

.. [1] K\. A\. Michalski and D\. Zheng, "Electromagnetic scattering and radiation by surfaces of arbitrary shape in layered media\. I\. Theory," in *IEEE Trans. Antennas Propag.*, vol\. 38, no\. 3, pp\. 335-344, March 1990\. 

