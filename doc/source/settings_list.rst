.. Author: Shashwat Sharma
.. Created on: Nov 07, 2021

.. _mgfsettings:

List of MGF Settings
====================

.. _basicsettings:

Basic Settings
--------------

method : {MGF_INTEGRATE, MGF_INTERPOLATE, MGF_DCIM, MGF_QUASISTATIC}
    The MGF computation technique [1]_, [2]_, [3]_.
    Default: MGF_INTERPOLATE.

extract_quasistatic : {true, false}
    Choose whether to extract quasistatic contributions in the spectral domain, and then add them back analytically in the spatial domain [3]_.
    Default: false.

extract_homogeneous : {true, false}
    Choose whether to extract the homogeneous material Green's function in the spectral domain, and then add it back in the spatial domain.
    Default: false.

extract_singularities : {true, false}
    Only relevant when either ``extract_quasistatic`` or ``extract_homogeneous`` is ``true``. Choose whether to extract weakly and strongly singular terms, and return the multiplicative coefficients necessary to later add back the singular terms. See the example ``testSingularity.cpp``.
    Default: *false*.

verbose : {true, false}
    Allow messages to be written to ``std::cout``.
    Default: *true*.

.. \item \lstinline{DCIM_method}:  One of \{\lstinline{DCIM_ONE_LEVEL}, \lstinline{DCIM_TWO_LEVEL}~\cite{DCIM01}, \lstinline{DCIM_THREE_LEVEL}~\cite{DCIM03}\}. Choose the flavour of DCIM to be used. Default: \lstinline{DCIM_TWO_LEVEL}.

.. \end{itemize}

.. \subsection{Advanced Settings}

.. \begin{itemize}

.. \item \lstinline{switching_point}: \lstinline{double} greater than 0. Only relevant for the \lstinline{method} \lstinline{MGF_INTEGRATE}. Set the extent along the positive real axis to which the Sommerfeld integration path is to be deformed, to avoid poles. This corresponds to the variable $a$ defined in~\cite{SI_PE}. Set to a negative number to use the default setting. Default: \lstinline{-1.0}.

.. \item \lstinline{tol_svd}: \lstinline{double} greater than 0. Only relevant for the \lstinline{method} \lstinline{MGF_DCIM}. Set the relative tolerance below which singular values are to be ignored in the application of the GPOF method~\cite{gpof} in the DCIM.\@ Default: \lstinline{1.0e-4}.

.. \item \lstinline{tol_eig}: \lstinline{double} greater than 0. Only relevant for the \lstinline{method} \lstinline{MGF_DCIM}. Set the relative tolerance below which eigenvalues are to be ignored in the application of the GPOF method~\cite{gpof} in the DCIM.\@ Default: \lstinline{1.0e-16}.

.. .. components : \lstinline{std::vector<bool>} of size 5. Choose which of the 5 unique formulation-C components~\cite{MGF02} should be computed. Default: \lstinline{1, 1, 1, 1, 1}.


References
----------

.. [1] K\. A\. Michalski and D\. Zheng, "Electromagnetic scattering and radiation by surfaces of arbitrary shape in layered media\. I\. Theory," in *IEEE Trans. Antennas Propag.*, vol\. 38, no\. 3, pp\. 335-344, March 1990\. 

.. [2] M\. I\. Aksun, "A robust approach for the derivation of closed-form Green's functions," in *IEEE Trans. Microw. Theory Tech.*, vol\. 44, no\. 5, pp\. 651-658, May 1996\.

.. [3] E\. Simsek, Q\. H\. Liu and B\. Wei, "Singularity subtraction for evaluation of Green's functions for multilayer media," in *IEEE Trans. Microw. Theory Tech.*, vol\. 54, no\. 1, pp. 216-225, Jan. 2006\.


