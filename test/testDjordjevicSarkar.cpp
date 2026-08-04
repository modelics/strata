// Unit test for the Djordjevic-Sarkar causal dielectric model.
//
// Verifies, through the public LayerManager API, that:
//   1. At the measurement frequency, the effective layer permittivity reproduces
//      the specified relative permittivity and loss tangent.
//   2. The loss tangent stays approximately flat across a wide band (1-10 GHz).
//   3. A dielectric-model layer is never merged away by MergeLayersWithSameMaterial.

#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>

#include "layers.hpp"
#include "constants.hpp"


static int failures = 0;

static void check(bool cond, const std::string &msg)
{
	if (!cond)
	{
		std::cout << "[FAIL] " << msg << std::endl;
		failures++;
	}
	else
	{
		std::cout << "[ OK ] " << msg << std::endl;
	}
}


int main()
{
	std::cout << "===========================" << std::endl;
	std::cout << "TestDjordjevicSarkar()" << std::endl;
	std::cout << "===========================" << std::endl;

	const double dk = 4.4;
	const double df = 0.02;
	const double f_ref = 1e9;

	// ------ Build a two-layer stack: one D-S layer, one constant layer ------

	DjordjevicSarkarParams ds;
	ds.relative_permittivity = dk;
	ds.loss_tangent = df;
	ds.at_frequency = f_ref;

	LayerManager lm;
	// D-S layer on top, constant layer below (different material)
	lm.AddLayer(0.0, 50e-6, {dk, 0.0}, 1.0, ds.conductivity_at_dc, 0.0, DielectricModelType::DJORDJEVIC_SARKAR, ds);
	lm.AddLayer(-50e-6, 0.0, {12.5, 0.0}, 1.0, 0.0, 0.0);

	// ------ Test 1: anchor conditions at the measurement frequency ------

	lm.ProcessLayers(f_ref);

	// Find the D-S layer (its dielectric_model flag survives)
	int ds_idx = -1;
	for (int i = 0; i < (int)lm.layers.size(); i++)
		if (lm.layers[i].dielectric_model == DielectricModelType::DJORDJEVIC_SARKAR)
			ds_idx = i;

	check(ds_idx >= 0, "D-S layer present after ProcessLayers (not merged away)");
	check(lm.layers.size() == 2, "layer count unchanged (no spurious merge)");

	// eps[] is the absolute complex permittivity; divide by eps0 for relative
	std::complex<double> epsr_eff = lm.eps[ds_idx] / strata::eps0;
	double re = epsr_eff.real();
	double tand = -epsr_eff.imag() / epsr_eff.real();

	std::cout << "  epsr(f_ref) = " << re << ", tan(delta) = " << tand << std::endl;

	check(std::abs(re - dk) / dk < 1e-6, "Re[epsr] matches relative_permittivity at f_ref");
	check(std::abs(tand - df) / df < 1e-3, "tan(delta) matches loss_tangent at f_ref");

	// ------ Test 2: loss tangent stays roughly flat across 1-10 GHz ------

	double tand_min = 1e300, tand_max = -1e300;
	for (double f = 1e9; f <= 10e9; f += 1e9)
	{
		lm.ProcessLayers(f);
		std::complex<double> e = lm.eps[ds_idx] / strata::eps0;
		double t = -e.imag() / e.real();
		tand_min = std::min(tand_min, t);
		tand_max = std::max(tand_max, t);
	}
	std::cout << "  tan(delta) over 1-10 GHz: [" << tand_min << ", " << tand_max << "]" << std::endl;
	check((tand_max - tand_min) / df < 0.10, "tan(delta) flat within 10% across 1-10 GHz");

	// ------ Summary ------

	std::cout << "===========================" << std::endl;
	if (failures == 0)
		std::cout << "All checks passed." << std::endl;
	else
		std::cout << failures << " check(s) FAILED." << std::endl;
	std::cout << "===========================" << std::endl;

	return failures == 0 ? 0 : 1;
}
