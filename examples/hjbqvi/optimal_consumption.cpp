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

// Uncomment for stochastic volatility
//#define SV 1

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

	#ifdef SV
	constexpr int Dimension = 3;
	#else
	constexpr int Dimension = 2;
	#endif

	// Dimensions
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

	// Scaled transaction cost (0 <= lambda < 1)
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

	// Maximum level of refinement
	// Solution and control data are printed at this level of refinement
	const int max_refinement = 6;

	// Used for both account dimensions
	const Real boundary = 200.;
	const int spatial_points = 20;
	const Axis axis = Axis::uniform(0., boundary, spatial_points);
	const Axis volatility_axis = Axis::uniform(0., 1., spatial_points);

	// Volatility parameters
	const Real volvol = 0.4;   // Volatility of volatility
	const Real kappa  = 5.;    // Mean reversion speed
	const Real theta  = sigma; // Mean
	const Real v_0    = theta;

	// Gets the value of an impulse control
	auto z_bounded = [=] (Real t, Real w, Real b, Real z_frac) {
		const Real tmp1 = b-c;
		const Real tmp2 = tmp1-boundary;
		const Real lo = max(
			-w,
			(tmp2 > 0) ? (tmp2/(1 - lambda)) : (tmp2/(1 + lambda))
		);
		const Real hi = min(
			boundary-w,
			(tmp1 > 0) ? (tmp1/(1 + lambda)) : (tmp1/(1 - lambda))
		);
		const Real z = z_frac * (hi - lo) + lo;
		if(fabs(z) < QuantPDE::epsilon) { return 0.; }
		return z;
	};

	#ifdef SV
	#define ARGS Real w, Real b, Real v
	#else
	#define ARGS Real w, Real b
	#endif

	// Problem description
	HJBQVI<
		Dimension,
		StochasticControlDimension,
		ImpulseControlDimension
	> hjbqvi(
		// Initial number of timesteps
		timesteps,

		// Initial spatial grid
		#ifdef SV
		{ axis, axis, volatility_axis },
		#else
		{ axis, axis },
		#endif

		// Initial stochastic control grid
		{ Axis::uniform(0., q_max, control_points) },

		// Initial impulse control grid
		{ Axis::uniform(0., 1., control_points) },

		// Expiry time
		T,

		// Discount factor
		[=] (Real t, ARGS) { return rho; },

		// Volatility
		{
			#ifdef SV
			[=] (Real t, ARGS) { return sqrt(v) * w; },
			[=] (Real t, ARGS) { return 0.; },
			[=] (Real t, ARGS) { return volvol * sqrt(v); }
			#else
			[=] (Real t, ARGS) { return sigma * w; },
			[=] (Real t, ARGS) { return 0.; }
			#endif
		},

		// Drift
		{
			[=] (Real t, ARGS, Real q) { return mu * w; },
			[=] (Real t, ARGS, Real q) {
				// No consumption at boundaries
				return
					r * b
					+ ((0 < b && b < boundary) ? -q : 0.)
				;
			}
			#ifdef SV
			, [=] (Real t, ARGS, Real q) {
				return kappa * (theta - v);
			}
			#endif
		},

		// Continuous cash/utility/etc. flow
		[=] (Real t, ARGS, Real q) {
			// No consumption at boundaries
			return (0 < b && b < boundary) ? pow(q, gamma)/gamma:0.;
		},

		// Impulse transition
		{
			[=] (Real t, ARGS, Real z_frac) {
				const Real z = z_bounded(t,w,b,z_frac);
				return w + z;
			},
			[=] (Real t, ARGS, Real z_frac) {
				const Real z = z_bounded(t,w,b,z_frac);
				return b - z - lambda * fabs(z) - c;
			}
			#ifdef SV
			, [=] (Real t, ARGS, Real z_frac) { return v; }
			#endif
		},

		// Impulse cash/utility/etc. reward
		[=] (Real t, ARGS, Real z_frac) {
			// Impulse
			const Real z = z_bounded(t,w,b,z_frac);
			return z == 0. ? -numeric_limits<Real>::infinity() : 0.;
		},

		// Cash/utility/etc. reward at expiry
		[=] (Real t, ARGS) {
			const Real liquidate = max(b + (1-lambda) * w - c, 0.);
			return pow(liquidate, gamma) / gamma;
		}
	);

	// Scheme to use
	hjbqvi.usePenalizedScheme();
	//hjbqvi.useDirectControlScheme();
	//hjbqvi.useSemiLagrangianScheme();

	// Linear solver to use
	hjbqvi.useBiCGSTABSolver();
	//hjbqvi.useSparseLUSolver();

	// Tell module that coefficients are time independent for optimization
	hjbqvi.coefficientsAreTimeIndependent();

	// Tell module to ignore controls that move the state off the grid
	hjbqvi.ignoreExtrapolatoryControls();

	// Run
	HJBQVI_main(
		hjbqvi,

		// Convergence test at (w=w_0, b=b_0, v=v_0)
		#ifdef SV
		{ w_0, b_0, v_0 },
		#else
		{ w_0, b_0 },
		#endif

		max_refinement
	);

	return 0;

}

