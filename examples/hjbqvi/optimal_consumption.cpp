////////////////////////////////////////////////////////////////////////////////
// optimal_consumption.cpp
// -----------------------
//
// The optimal consumption problem subject to transaction costs, as formulated
// in [1].
//
// [1] Chancelier, Jean-Philippe, Bernt Øksendal, and Agnès Sulem. "Combined
// stochastic control and optimal stopping, and application to numerical
// approximation of combined stochastic and impulse control." Proceedings of the
// Steklov Institute of Mathematics. VA Steklov. 237.0 (2002): 149-172.
//
// Author: Parsiad Azimzadeh
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>
#include <QuantPDE/Modules/HJBQVI>

using namespace QuantPDE;
using namespace QuantPDE::Modules;

////////////////////////////////////////////////////////////////////////////////

#include <cmath> // pow

using namespace std;

////////////////////////////////////////////////////////////////////////////////

int main() {

	// Dimensions
	constexpr int Dimension = 2;
	constexpr int StochasticControlDimension = 1;
	constexpr int ImpulseControlDimension = 1;

	// Expiry time (infinity for corresponding elliptic problem)
	const Real T = 40.;
	//const Real T = numeric_limits<Real>::infinity();

	// Discount factor
	const Real rho = 0.1;

	// Interest rate
	const Real r = 0.07;

	// Drift
	const Real mu = 0.11;

	// Volatility
	const Real sigma = 0.3;

	// 1 - relative risk aversion
	const Real gamma = 0.3;

	// Scaled transaction cost
	const Real lambda = 0.1;

	// Fixed transaction cost
	const Real c = 0.05;

	// Maximum withdrawal rate
	const Real q_max = 100.;

	// Initial values
	const Real w_0 = 45.2;
	const Real b_0 = 45.2;

	// Number of timesteps
	const int timesteps = 32;

	// Number of control points
	const int control_points = 16;

	// How to handle the control
	auto method = HJBQVIControlMethod::PENALTY_METHOD;
	//auto method = HJBQVIControlMethod::DIRECT_CONTROL;
	//auto method = HJBQVIControlMethod::EXPLICIT_CONTROL;

	// Maximum level of refinement
	// Solution and control data are printed at this level of refinement
	const int max_refinement = 6;

	// Used for both dimensions
	const Real boundary = 200.;
	const int spatial_points = 20;
	const Axis axis = Axis::uniform(0., boundary, spatial_points);

	// Gets the value of an impulse control
	auto zeta_bounded = [=] (Real t, Real w, Real b, Real zeta_frac) {
		const Real lo = max( -w, -(boundary-b+c)/(1-lambda) );
		const Real hi = min( (b-c)/(1+lambda), boundary-w );
		if(lo > hi) {
			// No impulse
			return 0.;
		}
		return zeta_frac * (hi - lo) + lo;
	};

	// Problem description
	HJBQVI<
		Dimension,
		StochasticControlDimension,
		ImpulseControlDimension
	> hjbqvi(
		// Initial number of timesteps
		timesteps,

		// Initial spatial grid
		{ axis, axis },

		// Initial stochastic control grid
		{ Axis::uniform(0., q_max, control_points) },

		// Initial impulse control grid
		{ Axis::uniform(0., 1., control_points) },

		// Expiry time
		T,

		// Discount factor
		[=] (Real t, Real w, Real b) { return rho; },

		// Volatility
		{
			[=] (Real t, Real w, Real b) { return sigma * w; },
			[=] (Real t, Real w, Real b) { return 0.; }
		},

		// Drift
		{
			[=] (Real t, Real w, Real b, Real q) { return mu * w; },
			[=] (Real t, Real w, Real b, Real q) {
				// No consumption at boundaries
				return
					r * b
					+ ((0 < b && b < boundary) ? -q : 0.)
				;
			}
		},

		// Continuous cash/utility/etc. flow
		[=] (Real t, Real w, Real b, Real q) {
			// No consumption at boundaries
			return (0 < b && b < boundary) ? pow(q, gamma)/gamma:0.;
		},

		// Impulse transition
		{
			[=] (Real t, Real w, Real b, Real zeta_frac) {
				const Real zeta = zeta_bounded(t,w,b,zeta_frac);
				return w + zeta;
			},
			[=] (Real t, Real w, Real b, Real zeta_frac) {
				const Real zeta = zeta_bounded(t,w,b,zeta_frac);
				return b - zeta - lambda * fabs(zeta) - c;
			}
		},

		// Impulse cash/utility/etc. reward
		[=] (Real t, Real w, Real b, Real zeta_frac) {

			// Prune bad controls for direct control scheme
			const Real zeta = zeta_bounded(t,w,b,zeta_frac);
			const Real w_plus = w + zeta;
			const Real b_plus = b - zeta - lambda * fabs(zeta) - c;
			if(
				fabs(w - w_plus) + fabs(b - b_plus)
						< QuantPDE::epsilon
				|| zeta == 0.
			) {
				return -numeric_limits<Real>::infinity();
			}
			// End prune

			return 0.;
		},

		// Cash/utility/etc. reward at expiry
		[=] (Real t, Real w, Real b) {
			const Real liquidate = max(b + (1-lambda) * w - c, 0.);
			return pow(liquidate, gamma) / gamma;
		}
	);

	// Scheme to use
	hjbqvi.usePenalizedScheme();
	//hjbqvi.useDirectControlScheme();
	//hjbqvi.useSemiLagrangianScheme();

	// Tell module that coefficients are time independent for optimization
	hjbqvi.coefficientsAreTimeIndependent();

	// Tell module to ignore controls that move the state off the grid
	hjbqvi.ignoreExtrapolatoryControls();

	// Run
	HJBQVI_main(
		hjbqvi,
		{ w_0, b_0 }, // Convergence test at (w=w_0, b=b_0)
		max_refinement
	);

	return 0;

}

