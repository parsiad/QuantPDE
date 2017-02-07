////////////////////////////////////////////////////////////////////////////////
// exchange_rate.cpp
// -----------------
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
	const int stochastic_control_points = 9;
	const int impulse_control_points = 16;
	const int points = 32;

	// Parameter that controls the clustering of nodes to the optimal parity
	const Real intensity = 10.;

	// Optimal parity
	const Real m = 0.;

	// Interest rate differential bounds
	const Real q_min = 0.;
	const Real q_max = 0.07;

	// Expiry time
	const Real T = 10.;
	//const Real T = numeric_limits<Real>::infinity();

	// Discount
	const Real rho = 0.02;

	// Volatility
	const Real sigma = 0.3;

	// Other parameters
	const Real a = 0.25;
	const Real b = 3.;
	const Real lambda = 1.;
	const Real c = 0.1;

	// Number of timesteps
	const int timesteps = 16;

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

	// q axis
	const Axis w = Axis::uniform( q_min, q_max, stochastic_control_points );

	// Impulse axis
	const Axis x_impulse = Axis::cluster(
		x_min,
		m,
		x_max,
		impulse_control_points,
		intensity
	);

	// Problem description
	HJBQVI<
		Dimension,
		StochasticControlDimension,
		ImpulseControlDimension
	> hjbqvi(
		// Initial number of timesteps
		timesteps,

		// Initial spatial grid
		{ x },

		// Initial stochastic control grid
		{ w },

		// Initial impulse control grid
		{ x_impulse },

		// Expiry time
		T,

		// Discount factor
		[=] (Real t, Real x) { return rho; },

		// Volatility
		{ [=] (Real t, Real x) { return sigma; } },

		// Drift
		{ [=] (Real t, Real x, Real q) { return -a * q; } },

		// Continuous cash/utility/etc. flow
		[=] (Real t, Real x, Real q) {
			const Real tmp = max(x - m, 0.);
			return - tmp * tmp - q * q * b;
		},

		// Impulse transition
		{ [=] (Real t, Real x, Real x_new) { return x_new; } },

		// Impulse cash/utility/etc. reward
		[=] (Real t, Real x, Real x_new) {
			const Real zeta = x_new - x;

			// Prune bad controls for direct control scheme
			if(zeta >= 0.) {
				return -numeric_limits<Real>::infinity();
			}
			// End Prune

			return - lambda * fabs(zeta) - c;
		},

		// Cash/utility/etc. reward at expiry
		[=] (Real t, Real x) { return 0.; }
	);

	// Scheme to use
	hjbqvi.usePenalizedScheme();
	//hjbqvi.useDirectControlScheme();
	//hjbqvi.useSemiLagrangianScheme();

	// Tell module that coefficients are time independent for optimization
	hjbqvi.coefficientsAreTimeIndependent();

	// Run
	HJBQVI_main(
		hjbqvi,
		{ m }, // Convergence test at optimal parity
		max_refinement
	);

	return 0;

}

