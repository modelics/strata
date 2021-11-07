.. Author: Shashwat Sharma
.. Created on: Aug 24, 2020

.. _egtech:

Example Technology File
=======================

Below is an annotated example of a layer definition file.

.. code-block:: yaml

	unit: um # applies to all spatial coordinates, and can be 'm', 'mm', 'um', or 'nm'

	dielectric_layers:

	    # The list of dielectric layers should include the dielectric fill corresponding to metal / via / ground layers.

	    # Layer properties can either be specified in this format...
	    L1:
	        zmin: 0
	        h: 50
	        epsr: 2.1
	        mur: 1
	        sigma: 0 # S/m

	    # ...or in this format
	    L2:
	        {zmin: -50, h: 50, epsr: 12.5, mur: 1, sigma: 0}

	# The following entries specify the half-spaces above ("top") and below ("bottom") the stack-up.
	# If an infinite PEC ground plane is desired, set the conductivity "sigma" to a negative number, e.g. -1. Usually, this would apply to the bottom halfspace.

	top_halfspace:
	    {epsr: 1, mur: 1, sigma: 0}

	bottom_halfspace:
	    {epsr: 1, mur: 1, sigma: -1}

