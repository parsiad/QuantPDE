////////////////////////////////////////////////////////////////////////////////
// unequal_borrowing_lending_rates.cpp
// -----------------------------------
//
// Prices a long/short position straddle option under the Black-Scholes model
// assuming unequal borrowing/lending rates [1].
//
// The pricing problem is given by
//
// V_t = \sup_{r \in \{ r_l, r_b \}} (-\sigma^2 S^2 V_{SS} / 2 + r( V - S V_S ))
// V(0, S) = max(S - K, K - S)
//
// for the short position. The \sup becomes an \inf for the long position.
//
// [1] Forsyth, Peter A., and George Labahn. "Numerical methods for controlled
// Hamilton-Jacobi-Bellman PDEs in finance." Journal of Computational Finance
// 11.2 (2007): 1.
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

#include <iostream> // cerr
#include <memory>   // unique_ptr
#include <numeric>  // accumulate

using namespace std;

////////////////////////////////////////////////////////////////////////////////

Configuration configuration;
Real T, r_l, r_b, vol, divs, S_0, K;
int N;
bool long_position;
RectilinearGrid1 *grid;

std::vector<Real> run(int k) {

	// 2^k
	Real factor = 1;
	for(int i = 0; i < k; ++i) { factor *= 2; }

	////////////////////////////////////////////////////////////////////////
	// Spatial grid
	////////////////////////////////////////////////////////////////////////

	// Refine grid R times
	auto refinedGrid = grid->refined(k);

	////////////////////////////////////////////////////////////////////////
	// Control grid
	////////////////////////////////////////////////////////////////////////

	// Control can be two interest rates: lending or borrowing
	RectilinearGrid1 controls( Axis { r_l, r_b } );

	////////////////////////////////////////////////////////////////////////
	// Payoff
	// ------
	//
	// Payoffs are lambda functions. The following is equivalent to
	//
	// auto payoff = [K] (Real S) {
	// 	return S < K ? K - S : S - K;
	// };
	////////////////////////////////////////////////////////////////////////

	auto payoff = straddlePayoff(K);

	////////////////////////////////////////////////////////////////////////
	// Iteration tree
	// --------------
	//
	// Sets up the loop structure:
	// for(int n = 0; n < N; ++n) {
	// 	for(int k = 0; ; ++k) {
	// 		// Solve a linear system
	// 		if(error < tolerance) break;
	// 	}
	// }
	////////////////////////////////////////////////////////////////////////

	ReverseConstantStepper stepper(
		0.,             // Initial time
		T,              // Expiry time
		T / N / factor  // Timestep size
	);
	ToleranceIteration tolerance;
	stepper.setInnerIteration(tolerance);

	////////////////////////////////////////////////////////////////////////
	// Linear system tree
	// ------------------
	//
	// Makes the linear system to solve at each iteration
	////////////////////////////////////////////////////////////////////////

	BlackScholes1 bs(
		refinedGrid,

		// Interest rate (passed as a control)
		Control1(refinedGrid),

		vol, // Volatility
		divs  // Dividend rate
	);

	// Policy iteration (a node in the linear system tree)
	//
	// The notation [N][_M] at the end of class names is used when the
	// problem is N-dimensional (in space) with an M-dimensional control.
	// In the following, we instantiate a policy iteration class meant for
	// one-dimensional (in space) problems with a one-dimensional control.
	//
	// Pick min/max policy iteration depending on whether we are considering
	// the long or short position problem.
	unique_ptr<IterationNode> policy(long_position
		? (IterationNode*)
			// sup[ V_tau - LV ]
			new MaxPolicyIteration1_1(refinedGrid, controls, bs)
		: (IterationNode*)
			// inf[ V_tau - LV ]
			new MinPolicyIteration1_1(refinedGrid, controls, bs)
	);
	policy->setIteration(tolerance); // Associate with k-iteration

	// BDF2 (timestepping)
	ReverseBDFTwo bdf2(refinedGrid, *policy);
	bdf2.setIteration(stepper); // Associate with n-iteration

	////////////////////////////////////////////////////////////////////////
	// Running
	// -------
	//
	// Everything prior to this was setup. Now we run the method.
	////////////////////////////////////////////////////////////////////////

	// Linear system solver
	BiCGSTABSolver solver;

	// Calculate the solution at time zero
	auto solution = stepper.solve(
		refinedGrid,   // Domain
		payoff, // Initial condition
		bdf2,   // Root of linear system tree
		solver  // Linear system solver
	);

	////////////////////////////////////////////////////////////////////////

	unsigned outer;
	Real inner = nan("");
	Real value;

	// Outer iterations
	outer = stepper.iterations()[0];

	// Average number of inner iterations
	auto its = tolerance.iterations();
	inner = accumulate(its.begin(), its.end(), 0.) / its.size();

	value = solution(S_0);

	return
		{
			(Real) refinedGrid.size(),
			(Real) outer,
			inner,
			value
		}
	;
}

int main(int argc, char **argv) {
	// Parse configuration file
	configuration = getConfiguration(argc, argv);

	// Get options
	int kn, k0;
	kn = getInt(configuration, "maximum_refinement", 5);
	k0 = getInt(configuration, "minimum_refinement", 0);
	T = getReal(configuration, "time_to_expiry", 1.);
	r_b = getReal(configuration, "interest_rate_long", .05);
	r_l = getReal(configuration, "interest_rate_short", .03);
	vol = getReal(configuration, "volatility", .3);
	divs = getReal(configuration, "dividend_rate", 0.);
	S_0 = getReal(configuration, "asset_price", 100.);
	K = getReal(configuration, "strike_price", 100.);
	N = getInt(configuration, "initial_number_of_timesteps", 12);
	long_position = getBool(configuration, "long_position", false);
	RectilinearGrid1 defGrid( (S_0 * Axis::special) + (K * Axis::special) );
	RectilinearGrid1 tmp = getGrid(configuration, "initial_grid", defGrid);
	grid = &tmp;

	// Print configuration file
	cerr << configuration << endl << endl;

	// Run and print results
	results(
		run,
		{
			"Nodes",
			"Steps",
			"Mean Inner Iterations",
			"Value"
		},
		kn, k0
	);

	return 0;
}

