////////////////////////////////////////////////////////////////////////////////
// gmwb.cpp
// --------
//
// Guaranteed minimum withdrawal benefits in variable annuities, as formulated
// in [1].
//
// [1] Dai, Min, Yue Kuen Kwok, and Jianping Zong. "Guaranteed minimum
// withdrawal benefit in variable annuities." Mathematical Finance 18.4 (2008):
// 595-611.
//
// [2] Chen, Zhuliang, and Peter A. Forsyth. "A numerical scheme for the impulse
// control formulation for pricing variable annuities with a guaranteed minimum
// withdrawal benefit (GMWB)." Numerische Mathematik 109.4 (2008): 535-569.
//
// Author: Parsiad Azimzadeh
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>
#include <QuantPDE/Modules/HJBQVI>

using namespace QuantPDE;
using namespace QuantPDE::Modules;

////////////////////////////////////////////////////////////////////////////////

#include <algorithm> // max

using namespace std;

////////////////////////////////////////////////////////////////////////////////

int main() {

	// Dimensions
	constexpr int Dimension = 2;
	constexpr int StochasticControlDimension = 1;
	constexpr int ImpulseControlDimension = 1;

	// Expiry time (infinity for corresponding elliptic problem)
	const Real T = 10.;
	//const Real T = numeric_limits<double>::infinity();

	// Interest rate
	const Real r = 0.05;

	// Premium
	const Real eta = 0.;

	// Volatility
	const Real sigma = 0.2;

	// Withdrawal rate
	const Real G = 10.;

	// Penalty rate
	const Real k = 0.1;

	// Transaction cost
	const Real c = 0.;

	// Initial payment
	const Real w_0 = 100.;

	// Number of timesteps
	const int timesteps = 32;

	// How to handle the control
	auto method = HJBQVIControlMethod::FULLY_IMPLICIT;
	//auto method = HJBQVIControlMethod::FULLY_EXPLICIT;

	// Maximum level of refinement
	// Solution and control data are printed at this level of refinement
	const int max_refinement = 6;

	// Problem description
	HJBQVI<
		Dimension,
		StochasticControlDimension,
		ImpulseControlDimension
	> hjbqvi(
		// Initial spatial grid
		{
			// Hand-picked axis, 65 points; used in [2]
			(w_0 / 100.) * Axis {
				0., 5., 10., 15., 20., 25.,
				30., 35., 40., 45.,
				50., 55., 60., 65., 70., 72.5, 75., 77.5, 80.,
						82., 84.,
				86., 88., 90., 91., 92., 93., 94., 95.,
				96., 97., 98., 99., 100.,
				101., 102., 103., 104., 105., 106.,
				107., 108., 109., 110., 112., 114.,
				116., 118., 120., 123., 126.,
				130., 135., 140., 145., 150., 160., 175., 200.,
						225.,
				250., 300., 500., 750., 1000.
			},

			// 51 points evenly spaced on [0, w_0]
			Axis::uniform(0., w_0, 51)
		},

		#if 0
		{
			Axis { 0., 50., 100. },
			Axis { 0., 50., 100. }
		},
		#endif

		// Initial stochastic control grid
		{ Axis { 0., G } }, // Does not get refined

		// Initial impulse control grid
		{ Axis { 0., 1. } }, // Does get refined

		// Expiry
		T,

		// Discount
		[=] (Real t, Real w, Real a) { return r; },

		// Volatility
		{
			[=] (Real t, Real w, Real a) { return sigma * w; },
			[=] (Real t, Real w, Real a) { return 0.; }
		},

		// Controlled drift
		{
			[=] (Real t, Real w, Real a, Real q) {
				// No continuous withdrawal at a == 0
				return a > QuantPDE::epsilon ? -q : 0.;
			},
			[=] (Real t, Real w, Real a, Real q) {
				// No continuous withdrawal at a == 0
				return a > QuantPDE::epsilon ? -q : 0.;
			}
		},

		// Uncontrolled drift
		{
			[=] (Real t, Real w, Real a) { return (r-eta)*w; },
			[=] (Real t, Real w, Real a) { return 0.; }
		},

		// Controlled continuous flow
		[=] (Real t, Real w, Real a, Real q) {
			// No continuous withdrawal at a == 0
			return a > QuantPDE::epsilon ? q : 0.;
		},

		// Uncontrolled continuous flow
		[=] (Real t, Real w, Real a) { return 0.; },

		// Transition
		{
			[=] (Real t, Real w, Real a, Real zeta_frac) {
				const Real zeta = zeta_frac * a;
				return max(w - zeta, 0.);
			},
			[=] (Real t, Real w, Real a, Real zeta_frac) {
				const Real zeta = zeta_frac * a;
				return a - zeta;
			}
		},

		// Impulse flow
		[=] (Real t, Real w, Real a, Real zeta_frac) {
			const Real zeta = zeta_frac * a;
			return (1-k) * zeta - c;
		},

		// Exit function
		[=] (Real t, Real w, Real a) { return max(w, (1-k) * a - c); },

		// Initial number of timesteps
		timesteps,

		// How to handle the control
		method,

		// Bounded domain
		false,

		// Do not refine stochastic control grid
		false
	);

	// Run
	HJBQVI_main(
		hjbqvi,
		{ w_0, w_0 }, // Convergence test at (w=w_0, a=w_0)
		max_refinement
	);

	return 0;

}

