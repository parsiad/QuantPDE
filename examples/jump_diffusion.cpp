////////////////////////////////////////////////////////////////////////////////
// jump_diffusion.cpp
// ------------------
//
// Computes the price of a European call with jump-diffusion driven by a Poisson
// process. The jump amplitude is assumed to be lognormally distributed.
//
// Author: Parsiad Azimzadeh
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>
#include <QuantPDE/Modules/Lambdas>
#include <QuantPDE/Modules/Operators>
#include <QuantPDE/Modules/Utilities>

using namespace QuantPDE;
using namespace QuantPDE::Modules;

///////////////////////////////////////////////////////////////////////////////

#include <iostream> // cerr
#include <memory>   // unique_ptr

using namespace std;

///////////////////////////////////////////////////////////////////////////////

enum JumpType { merton, kou };

JumpType jump_type;
Real
	T, r, vol, divs, S_0, K,
	jump_arrival_rate,
	jump_mean, jump_std,
	jump_up_probability,
	jump_up_mean_reciprocal, jump_down_mean_reciprocal
;
int N;
bool call, digital, variable;
RectilinearGrid1 *grid;

ResultsTuple1 run(int k) {

	// 2^k
	int factor = 1;
	for(int i = 0; i < k; ++i) {
		factor *= 2;
	}

	////////////////////////////////////////////////////////////////////////
	// Spatial grid
	////////////////////////////////////////////////////////////////////////

	// Refine grid R times
	auto refined_grid = grid->refined(k);

	////////////////////////////////////////////////////////////////////////
	// Payoff
	// ------
	//
	// Payoffs are lambda functions. The following is equivalent to
	//
	// auto payoff = [K] (Real S) {
	// 	return S < K ? K - S : 0.;
	// };
	////////////////////////////////////////////////////////////////////////

	auto payoff = callPayoff(K);

	////////////////////////////////////////////////////////////////////////
	// Iteration tree
	// --------------
	//
	// Sets up the loop structure:
	// for(int n = 0; n < N; ++n) {
	// 	// Solve a linear system
	// }
	////////////////////////////////////////////////////////////////////////

	ReverseConstantStepper stepper(
		0.,              // Initial time
		T,               // Expiry time
		T / (N * factor) // Timestep size
	);

	////////////////////////////////////////////////////////////////////////
	// Linear system tree
	// ------------------
	//
	// Makes the linear system to solve at each iteration
	////////////////////////////////////////////////////////////////////////

	auto jump_density = (jump_type == JumpType::merton) ?
		lognormal(
			jump_mean,
			jump_std
		) :
		doubleExponential(
			jump_up_probability,
			jump_up_mean_reciprocal,
			jump_down_mean_reciprocal
		)
	;

	BlackScholesJumpDiffusion1 bs(
		refined_grid,
		r,    // Interest
		vol,  // Volatility
		divs, // Dividend rate
		jump_arrival_rate, jump_density
	);
	bs.setIteration(stepper);

	typedef ReverseBDFOne Discretization;
	Discretization discretization(refined_grid, bs);
	discretization.setIteration(stepper);

	////////////////////////////////////////////////////////////////////////
	// Running
	// -------
	//
	// Everything prior to this was setup. Now we run the method.
	////////////////////////////////////////////////////////////////////////

	// Linear system solver
	SparseLUSolver solver;

	auto solution = stepper.solve(
		refined_grid,    // Domain
		payoff,         // Initial condition
		discretization, // Root of linear system tree
		solver          // Linear system solver
	);

	////////////////////////////////////////////////////////////////////////

	// Timesteps
	unsigned timesteps = stepper.iterations()[0];

	return ResultsTuple1(
		{ (Real) refined_grid.size(), (Real) timesteps, },
		solution, S_0
	);
}

int main(int argc, char **argv) {
	// Parse configuration file
	Configuration configuration = getConfiguration(argc, argv);

	// Get options
	int kn, k0;
	Real S_max, S_min, dS;
	kn = getInt(configuration, "maximum_refinement", 8);
	k0 = getInt(configuration, "minimum_refinement", 3);
	T = getReal(configuration, "time_to_expiry", .25);
	r = getReal(configuration, "interest_rate", .05);
	vol = getReal(configuration, "volatility", .15);
	divs = getReal(configuration, "dividend_rate", 0.);
	S_0 = getReal(configuration, "asset_price", 100.);
	K = getReal(configuration, "strike_price", 100.);
	S_min = getReal(configuration, "print_asset_price_minimum", 0.);
	S_max = getReal(configuration, "print_asset_price_maximum", S_0 * 2.);
	dS = getReal(configuration, "print_asset_price_step_size", S_0 / 10.);
	jump_arrival_rate = getReal(configuration, "jump_arrival_rate", .1);
	N = getInt(configuration, "initial_number_of_timesteps", 12);
	RectilinearGrid1 default_grid( (S_0 * Axis::special) + (K * Axis::special) + (S_0 * Axis { 1000. }) + (K * Axis { 1000. }) );
	RectilinearGrid1 tmp = getGrid(configuration, "initial_grid", default_grid);
	grid = &tmp;
	string jump_type_str = getString(configuration, "jump_type", "lognormal");
	if(jump_type_str == "lognormal") {
		jump_type = JumpType::merton;
		jump_mean = getReal(configuration, "jump_amplitude_mean", -.1);
		jump_std = getReal(configuration, "jump_amplitude_deviation", .45);
	} else if(jump_type_str == "double_exponential") {
		jump_type = JumpType::kou;
		jump_up_probability = getReal(configuration, "jump_up_probability", .3445);
		jump_up_mean_reciprocal = getReal(configuration, "jump_up_mean_reciprocal", 3.0465);
		jump_down_mean_reciprocal = getReal(configuration, "jump_down_mean_reciprocal", 3.0775);
		// Exact solution with default parameters is 3.973479 at S = 100
	} else {
		cerr << "error: jump_model can either be \"lognormal\" or \"double_exponential\"" << endl;
		return 2;
	}

	// Print configuration file
	cerr << configuration << endl << endl;

	// Run and print results
	ResultsBuffer1 buffer(
		run,
		{ "Nodes", "Steps" },
		kn, k0
	);
	buffer.addPrintGrid( RectilinearGrid1(Axis::range(S_min, dS, S_max)) );
	buffer.stream();

	return 0;
}

