////////////////////////////////////////////////////////////////////////////////
// stock_loan.cpp
// --------------
//
// Author: Parsiad Azimzadeh
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>
#include <QuantPDE/Modules/Lambdas>
#include <QuantPDE/Modules/Operators>
#include <QuantPDE/Modules/Utilities>

using namespace QuantPDE;
using namespace QuantPDE::Modules;

////////////////////////////////////////////////////////////////////////////////

#include <algorithm> // max, min
#include <cassert>   // assert
#include <cmath>     // exp
#include <limits>    // numeric_limits
#include <numeric>   // accumulate

using namespace std;

////////////////////////////////////////////////////////////////////////////////

Real
	T, r, vol, S_0,
	q_hat, xi, rho, P, beta, eta, arrival, theta,
	up_probability, up_mean_r, down_mean_r
;
int N;
bool explicit_impulses, margin_calls_enabled;
RectilinearGrid1 *grid;

////////////////////////////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////////////////////////////

constexpr Real inf = numeric_limits<Real>::infinity();

/**
 * Interpolation function for similarity reduction.
 * @param v The value of the stock loan for fixed q_hat.
 * @param s The value of the collateral.
 * @param q The outstanding loan.
 * @param q_hat A fixed, positive number.
 * @return The interpolated value of v at (s, q).
 */
inline Real simv(const Interpolant1 &v, Real s, Real q, Real q_hat) {
	assert(q > 0.);
	const Real alpha = q_hat / q;
	return v(alpha * s) / alpha;
}

/**
 * Computes the loan growth factor for some time.
 * @param t A time.
 * @return The growth factor.
 */
inline Real growth(Real t) { return exp((r+xi) * t); }

////////////////////////////////////////////////////////////////////////////////
// Impulse functions
////////////////////////////////////////////////////////////////////////////////

Real repay(Real s) {
	return max(s - growth(T) * q_hat, 0.);
}

Real prepay(Real t, Real s) {
	return (t >= P) ? max(s - rho * growth(t) * q_hat, 0. ) : (-inf);
}

Real liquidation(Real t, Real s) {
	const Real g = growth(t);
	return (g * q_hat >= beta * s) ? max(s - g * q_hat, 0.) : inf;
}

Real margin_call(const Interpolant1 &v, Real t, Real s, Real z) {
	const Real g = growth(t);
	return (g * q_hat >= eta * s) ?
		(simv(v, s, z * s / g, q_hat) - (g * q_hat - z * s)) :
		inf
	;
}

////////////////////////////////////////////////////////////////////////////////

ResultsTuple1 run(int k) {

	// 2^k
	int factor = 1;
	for(int i = 0; i < k; ++i) {
		factor *= 4;
	}
	const int Nck = N * factor;

	////////////////////////////////////////////////////////////////////////
	// Spatial grid
	////////////////////////////////////////////////////////////////////////

	auto refined_grid = grid->refined(k); // Refine grid R times

	////////////////////////////////////////////////////////////////////////
	// Iteration tree
	////////////////////////////////////////////////////////////////////////

	const Real dt = T / Nck;
	ReverseConstantStepper stepper(0., T, dt);
	ToleranceIteration tolerance;

	// Explicit events
	for(int n = 1; n < Nck; ++n) { // Skip n = Nck as it has no effect
		const Real t = n * dt;

		if(explicit_impulses) {
			// Client
			stepper.add(
				t,
				[=] (const Interpolant1 &v, Real s) {
					const Real v0 = v(s);
					const Real v1 = prepay(t, s);
					return max(v0, v1);
				},
				refined_grid
			);
		}

		// Lender
		stepper.add(
			t,
			[=] (const Interpolant1 &v, Real s) {
				const Real v0 = v(s);
				const Real v1 = liquidation(t, s);
				const Real v2 = margin_calls_enabled ?
					margin_call(v, t, s, theta) :
					inf
				;
				return min(v0, min(v1, v2));
			},
			refined_grid
		);
	}

	////////////////////////////////////////////////////////////////////////
	// Linear system tree
	////////////////////////////////////////////////////////////////////////

	IterationNode *root;

	BlackScholesJumpDiffusion1 bs(
		refined_grid,

		r,   // Interest
		vol, // Volatility
		0.,  // Dividend rate

		arrival, // Mean arrival time
		doubleExponential(up_probability, up_mean_r, down_mean_r)
	);
	bs.setIteration(stepper);

	typedef ReverseBDFOne Discretization;
	Discretization discretization(refined_grid, bs);
	discretization.setIteration(stepper);

	PenaltyMethodDifference1 penalty(
		refined_grid,
		discretization,
		prepay
	);
	penalty.setIteration(tolerance);

	if(explicit_impulses) {
		root = &discretization;
	} else {
		stepper.setInnerIteration(tolerance);
		root = &penalty;
	}

	////////////////////////////////////////////////////////////////////////
	// Running
	////////////////////////////////////////////////////////////////////////

	SparseLUSolver solver;

	auto solution = stepper.solve(
		refined_grid, // Domain
		repay,        // Initial condition
		*root,        // Root of linear system tree
		solver        // Linear system solver
	);

	////////////////////////////////////////////////////////////////////////

	// Timesteps
	unsigned timesteps = stepper.iterations()[0];

	// Average number of policy iterations
	Real policy_its = nan("");
	auto its = tolerance.iterations();
	policy_its = accumulate(its.begin(), its.end(), 0.) / its.size();

	return ResultsTuple1(
		{(Real) refined_grid.size(), (Real) timesteps, policy_its },
		solution, S_0
	);
}

int main(int argc, char **argv) {
	// Parse configuration file
	Configuration configuration = getConfiguration(argc, argv);

	// Get options
	int kn, k0;
	Real S_max, S_min, dS;
	kn = getInt(configuration, "maximum_refinement", 7);
	k0 = getInt(configuration, "minimum_refinement", 0);
	T = getReal(configuration, "time_to_expiry", 10.);
	r = getReal(configuration, "interest_rate", .02);
	vol = getReal(configuration, "volatility", .3);
	S_0 = getReal(configuration, "asset_price", 100.);
	q_hat = getReal(configuration, "loan_value", 80.);
	xi = getReal(configuration, "spread", 0.);
	rho = getReal(configuration, "penalty_scaling", 1.);
	P = getReal(configuration, "lockout_time", 0.);
	beta = getReal(configuration, "liquidation_trigger", 80. / 89.);
	eta = getReal(configuration, "margin_call_trigger", 80. / 94.);
	arrival = getReal(configuration, "jump_arrival_rate", 0.2);
	up_probability = getReal(configuration, "jump_up_probability", 0.09);
	up_mean_r = getReal(configuration, "jump_up_mean_reciprocal", 2.3);
	down_mean_r = getReal(configuration, "jump_down_mean_reciprocal", 1.8);
	theta = getReal(configuration, "margin_call_loan_to_value_ratio", .8);
	S_min = getReal(configuration, "print_asset_price_minimum", 0.);
	S_max = getReal(configuration, "print_asset_price_maximum", 300.);
	dS = getReal(configuration, "print_asset_price_step_size", 10.);
	N = getInt(configuration, "initial_number_of_timesteps", 2);
	explicit_impulses = getInt(configuration, "handle_all_impulses_explicitly", false);
	margin_calls_enabled = getInt(configuration, "margin_calls_enabled", false);
	RectilinearGrid1 default_grid( S_0 * Axis::special );
	RectilinearGrid1 tmp = getGrid(configuration, "initial_grid", default_grid);
	grid = &tmp;

	// Print configuration file
	cerr << configuration << endl << endl;

	// Run and print results
	ResultsBuffer1 buffer(
		run,
		{ "Nodes", "Steps", "Mean Policy Iterations" },
		kn, k0
	);
	buffer.addPrintGrid( RectilinearGrid1(Axis::range(S_min, dS, S_max)) );
	buffer.stream();

	return 0;
}

