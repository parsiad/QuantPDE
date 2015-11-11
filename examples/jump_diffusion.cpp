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

Real T, r, vol, divs, S_0, K, arrival, jump_mean, jump_std;
int N;
bool call, digital, american, variable;
RectilinearGrid1 *grid;

ResultsTuple1 run(int k) {

	// 2^k or 2^(2k)
	int factor = 1;
	for(int i = 0; i < k; ++i) {
		if(american) { factor *= 4; }
		else         { factor *= 2; }
	}

	////////////////////////////////////////////////////////////////////////
	// Spatial grid
	////////////////////////////////////////////////////////////////////////

	// Refine grid R times
	auto refinedGrid = grid->refined(k);

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

	BlackScholesJumpDiffusion1 bs(
		refinedGrid,

		r,    // Interest
		vol,  // Volatility
		divs, // Dividend rate

		arrival, // Mean arrival time
		lognormal(jump_mean, jump_std) // Log-normal probability density
	);
	bs.setIteration(stepper);

	typedef ReverseBDFOne Discretization;
	Discretization discretization(refinedGrid, bs);
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
		refinedGrid,    // Domain
		payoff,         // Initial condition
		discretization, // Root of linear system tree
		solver          // Linear system solver
	);

	////////////////////////////////////////////////////////////////////////

	// Timesteps
	unsigned timesteps = stepper.iterations()[0];

	return ResultsTuple1(
		{ (Real) refinedGrid.size(), (Real) timesteps, },
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
	T = getReal(configuration, "time_to_expiry", 1.);
	r = getReal(configuration, "interest_rate", .05);
	vol = getReal(configuration, "volatility", .3);
	divs = getReal(configuration, "dividend_rate", 0.);
	S_0 = getReal(configuration, "asset_price", 100.);
	K = getReal(configuration, "strike_price", 100.);
	S_min = getReal(configuration, "print_asset_price_minimum", 0.);
	S_max = getReal(configuration, "print_asset_price_maximum", S_0 * 2.);
	dS = getReal(configuration, "print_asset_price_step_size", S_0 / 10.);
	arrival = getReal(configuration, "jump_arrival_rate", .05);
	jump_mean = getReal(configuration, "jump_amplitude_mean", -.8);
	jump_std = getReal(configuration, "jump_amplitude_deviation", .42);
	N = getInt(configuration, "initial_number_of_timesteps", 12);
	RectilinearGrid1 defGrid( (S_0 * Axis::special) + (K * Axis::special) );
	RectilinearGrid1 tmp = getGrid(configuration, "initial_grid", defGrid);
	grid = &tmp;

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

