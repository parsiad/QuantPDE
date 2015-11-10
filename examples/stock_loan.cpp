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

using namespace std;

////////////////////////////////////////////////////////////////////////////////

Real
	T, r, vol, S_0,
	q_hat, xi, rho, P, beta, eta, arrival, jump_mean, jump_std, theta
;
int N;
RectilinearGrid1 *grid;

/**
 * Interpolation function for similarity reduction.
 * @param v The value of the stock loan for fixed q_hat.
 * @param s The value of the collateral.
 * @param q The outstanding loan.
 * @param q_hat A fixed, positive number.
 */
inline Real similarity_value(
	const Interpolant1 &v,
	Real s, Real q,
	Real q_hat
) {
	assert(q > 0.);
	const Real alpha = q_hat / q;
	return v(alpha * s) / alpha;
}

ResultsTuple1 run(int k) {

	// 2^k
	int factor = 1;
	for(int i = 0; i < k; ++i) {
		factor *= 2;
	}

	////////////////////////////////////////////////////////////////////////
	// Spatial grid
	////////////////////////////////////////////////////////////////////////

	auto refinedGrid = grid->refined(k); // Refine grid R times

	////////////////////////////////////////////////////////////////////////
	// Payoffs and impulses
	////////////////////////////////////////////////////////////////////////

	auto repay = [=] (Real s) {
		const Real grow = exp((r+xi) * T);
		return max(s - grow * q_hat, 0.);
	};

	auto prepay = [=] (Real t, Real s) {
		const Real grow = exp((r+xi) * t);
		if(t >= P) { return max(s - rho * grow * q_hat, 0.); }
		else { return -numeric_limits<Real>::infinity(); }
	};

	auto liquidation = [=] (Real t, Real s) {
		const Real grow = exp((r+xi) * t);
		if(grow * q_hat >= beta * s) {
			return max(s - grow * q_hat, 0.);
		} else { return numeric_limits<Real>::infinity(); }
	};

	/*
	auto margin_call = [=] (const Interpolant1 &v, Real t, Real s, Real z) {
		const Real grow = exp((r+xi) * t);
		if(grow * q_hat >= eta * s) {
			return similarity_value(v, s, z * s / grow, q_hat)
					- (grow * q_hat - z * s);
		} else { return numeric_limits<Real>::infinity(); }
	};
	*/

	////////////////////////////////////////////////////////////////////////
	// Iteration tree
	////////////////////////////////////////////////////////////////////////

	const Real dt = T / (N * factor);
	ReverseConstantStepper stepper(0., T, dt);
	ToleranceIteration tolerance;
	stepper.setInnerIteration(tolerance);

	// Lender's events (handled explicitly)
	for(int n = 1; n < N; ++n) {
		const Real t = n * dt;
		stepper.add(
			n * dt,
			[=] (const Interpolant1 &v, Real s) {
				const Real v0 = v(s);
				const Real v1 = liquidation(t, s);
				//const Real v2 = margin_call(v, t, s, theta);
				//return min(v0, min(v1, v2));
				return min(v0, v1);
			},
			refinedGrid
		);
	}

	////////////////////////////////////////////////////////////////////////
	// Linear system tree
	////////////////////////////////////////////////////////////////////////

	BlackScholesJumpDiffusion1 bs(
		refinedGrid,

		r,   // Interest
		vol, // Volatility
		0.,  // Dividend rate

		arrival, // Mean arrival time
		lognormal(jump_mean, jump_std) // Log-normal probability density
	);
	bs.setIteration(stepper);

	typedef ReverseBDFOne Discretization;
	Discretization discretization(refinedGrid, bs);
	discretization.setIteration(stepper);

	PenaltyMethodDifference1 penalty(refinedGrid, discretization, prepay);
	penalty.setIteration(tolerance);

	////////////////////////////////////////////////////////////////////////
	// Running
	////////////////////////////////////////////////////////////////////////

	BiCGSTABSolver solver;

	auto solution = stepper.solve(
		refinedGrid, // Domain
		repay,       // Initial condition
		penalty,     // Root of linear system tree
		solver       // Linear system solver
	);

	////////////////////////////////////////////////////////////////////////

	// Timesteps
	unsigned timesteps = stepper.iterations()[0];

	// Average number of policy iterations
	Real policyIts = nan("");
	auto its = tolerance.iterations();
	policyIts = accumulate(its.begin(), its.end(), 0.) / its.size();

	return ResultsTuple1(
		{(Real) refinedGrid.size(), (Real) timesteps, policyIts },
		solution, S_0
	);
}

int main(int argc, char **argv) {
	// Parse configuration file
	Configuration configuration = getConfiguration(argc, argv);

	// Get options
	int kn, k0;
	kn = getInt(configuration, "maximum_refinement", 5);
	k0 = getInt(configuration, "minimum_refinement", 0);
	T = getReal(configuration, "time_to_expiry", 1.);
	r = getReal(configuration, "interest_rate", .04);
	vol = getReal(configuration, "volatility", .2);
	S_0 = getReal(configuration, "asset_price", 100.);
	q_hat = getReal(configuration, "loan_value", 80.);
	xi = getReal(configuration, "spread", 0.);
	P = getReal(configuration, "lockout_time", 0.);
	beta = getReal(configuration, "liquidation_trigger", 80. / 89.);
	eta = getReal(configuration, "margin_call_trigger", 80. / 94.);
	arrival = getReal(configuration, "jump_arrival_rate", .05);
	jump_mean = getReal(configuration, "jump_amplitude_mean", -.8);
	jump_std = getReal(configuration, "jump_amplitude_deviation", .42);
	theta = getReal(configuration, "margin_call_loan_to_value_ratio", .8);
	N = getInt(configuration, "initial_number_of_timesteps", 12);
	RectilinearGrid1 defGrid( (S_0*Axis::special) + (q_hat*Axis::special) );
	RectilinearGrid1 tmp = getGrid(configuration, "initial_grid", defGrid);
	grid = &tmp;

	// Print configuration file
	cerr << configuration << endl << endl;

	// Run and print results
	ResultsBuffer1 buffer(
		run,
		{ "Nodes", "Steps", "Mean Policy Iterations" },
		kn, k0
	);
	buffer.addPrintGrid( RectilinearGrid1(Axis::range(0., 10., 300.)) );
	buffer.stream();

	return 0;
}

