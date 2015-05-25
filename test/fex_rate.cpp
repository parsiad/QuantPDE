////////////////////////////////////////////////////////////////////////////////
// fex_rate.cpp
// ------------
//
// Optimal combined control of the exchange rate, as formulated in [1].
//
// [1] Mundaca, Gabriela, and Bernt Øksendal. "Optimal stochastic intervention
// control with application to the exchange rate." Journal of Mathematical
// Economics 29.2 (1998): 225-243.
//
// Author: Parsiad Azimzadeh
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>
#include <QuantPDE/Modules/HJBQVI>

using namespace QuantPDE;
using namespace QuantPDE::Modules;

////////////////////////////////////////////////////////////////////////////////

#include <cmath>  // fabs
#include <limits> // numeric_limits

using namespace std;

////////////////////////////////////////////////////////////////////////////////

int main() {

	// Dimensions
	constexpr int Dimension = 1;
	constexpr int StochasticControlDimension = 1;
	constexpr int ImpulseControlDimension = 1;

	// Boundaries
	const Real x_min = -6.;
	const Real x_max = +6.;

	// Number of points to use in discretizations
	const int points = 32;

	// Parameter that controls the clustering of nodes to the optimal parity
	const Real intensity = 10.;

	// Optimal parity
	const Real m = 0.;

	// Interest rate differential bounds
	const Real q_min = 0.;
	const Real q_max = 0.05;

	// Expiry time (infinity for corresponding elliptic problem)
	//const Real T = 10.;
	const Real T = numeric_limits<double>::infinity();

	// Discount
	const Real rho = 0.02;

	// Volatility
	const Real sigma = 0.3;

	// Other parameters
	const Real a = 0.25;
	const Real b = 3.;
	const Real lambda = 1.;
	const Real c = 0.;

	// Number of timesteps
	const int timesteps = 32;

	// How to handle the control
	auto method = HJBQVIControlMethod::FULLY_IMPLICIT;
	//auto method = HJBQVIControlMethod::FULLY_EXPLICIT;

	// Maximum level of refinement
	// Solution and control data are printed at this level of refinement
	const int max_refinement = 6;

	// x axis
	const Axis x = Axis::cluster(
		x_min,
		m,
		x_max,
		points,
		intensity
	);

	// Problem description
	HJBQVI<
		Dimension,
		StochasticControlDimension,
		ImpulseControlDimension
	> hjbqvi(
		// Initial spatial grid
		{ x },

		// Initial stochastic control grid
		{ Axis::uniform( q_min, q_max, points ) },

		// Initial impulse control grid
		{ x },

		// Expiry
		T,

		// Discount
		[=] (Real t, Real x) {
			return rho;
		},

		// Volatility
		{
			[=] (Real t, Real x) {
				return sigma;
			}
		},

		// Controlled drift
		{
			[=] (Real t, Real x, Real q) {
				return -a * q;
			}
		},

		// Uncontrolled drift
		{
			[=] (Real t, Real x) {
				return 0.;
			}
		},

		// Controlled continuous flow
		[=] (Real t, Real x, Real q) {
			return -q * q * b;
		},

		// Uncontrolled continuous flow
		[=] (Real t, Real x) {
			const Real tmp = max(x - m, 0.);
			return -tmp * tmp;
		},

		// Transition
		{
			[=] (Real t, Real x, Real x_new) {
				return x_new;
			}
		},

		// Impulse flow
		[=] (Real t, Real x, Real x_new) {
			const Real zeta = x_new - x;
			return - lambda * fabs(zeta) - c;
		},

		// Exit function
		[=] (Real t, Real x) {
			return 0.;
		},

		// Initial number of timesteps
		timesteps,

		// How to handle the control
		method,

		// Bounded domain?
		false,

		// Refine stochastic control?
		true,

		// Refine impulse control?
		true,

		// Are the coefficients time independent?
		true
	);

	// Run
	HJBQVI_main(
		hjbqvi,
		{ m }, // Convergence test at optimal parity
		max_refinement
	);

	return 0;

}

