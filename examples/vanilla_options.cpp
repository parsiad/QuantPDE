////////////////////////////////////////////////////////////////////////////////
// convergence_table.cpp
// ---------------------
//
// Outputs the rate of convergence for computing the price of a
// European/American digital/nondigital call/put.
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
#include <numeric>  // accumulate

using namespace std;

///////////////////////////////////////////////////////////////////////////////

Real T, r, vol, divs, S_0, K;
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

	Function1 payoff = digital ?
		  ( call ? digitalCallPayoff(K) : digitalPutPayoff(K) )
		: (call ? callPayoff(K) : putPayoff(K))
	;

	auto refinedGrid = grid->refined(k);

	// Black-Scholes operator (L in V_t = LV)
	BlackScholes1 bs(
		refinedGrid,
		r, vol, divs
	);

	// Timestepping method
	unique_ptr<Iteration> stepper(variable
		? (Iteration *) new ReverseVariableStepper(
			0.,               // startTime
			T,                // endTime
			T / (N * factor), // dt
			1. / factor       // target
		)
		: (Iteration *) new ReverseConstantStepper(
			0.,              // startTime
			T,               // endTime
			T / (N * factor) // dt
		)
	);

	// Time discretization method
	typedef ReverseRannacher Discretization;
	Discretization discretization(refinedGrid, bs);
	discretization.setIteration(*stepper);

	// American-specific components; penalty method or not?
	IterationNode *root;
	unique_ptr<ToleranceIteration> tolerance;
	unique_ptr<PenaltyMethodDifference1> penalty;
	if(american) {
		// American case
		penalty = unique_ptr<PenaltyMethodDifference1>(
			new PenaltyMethodDifference1(
				refinedGrid,
				discretization,
				payoff
			)
		);

		tolerance = unique_ptr<ToleranceIteration>(
			new ToleranceIteration()
		);

		penalty->setIteration(*tolerance);
		stepper->setInnerIteration(*tolerance);

		root = penalty.get();
	} else {
		// European case
		root = &discretization;
	}

	// Linear system solver
	SparseLUSolver solver;

	// Compute solution
	auto solution = stepper->solve(
		refinedGrid,
		payoff,
		*root,
		solver
	);

	// Timesteps
	unsigned timesteps = stepper->iterations()[0];

	// Average number of policy iterations
	Real policyIts = nan("");
	if(american) {
		auto its = tolerance->iterations();
		policyIts = accumulate(its.begin(), its.end(), 0.) / its.size();
	}

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
	Real S_max, S_min, dS;
	kn = getInt(configuration, "maximum_refinement", 5);
	k0 = getInt(configuration, "minimum_refinement", 0);
	T = getReal(configuration, "time_to_expiry", 1.);
	r = getReal(configuration, "interest_rate", .04);
	vol = getReal(configuration, "volatility", .2);
	divs = getReal(configuration, "dividend_rate", 0.);
	S_0 = getReal(configuration, "asset_price", 100.);
	K = getReal(configuration, "strike_price", 100.);
	S_min = getReal(configuration, "print_asset_price_minimum", 0.);
	S_max = getReal(configuration, "print_asset_price_maximum", S_0 * 2.);
	dS = getReal(configuration, "print_asset_price_step_size", S_0 / 10.);
	N = getInt(configuration, "initial_number_of_timesteps", 12);
	call = getBool(configuration, "is_call", true);
	digital = getBool(configuration, "is_digital", false);
	american = getBool(configuration, "is_american", false);
	variable = getBool(configuration, "use_variable_timestepping", false);
	RectilinearGrid1 defGrid( (S_0 * Axis::special) + (K * Axis::special) );
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
	buffer.addPrintGrid( RectilinearGrid1(Axis::range(S_min, dS, S_max)) );
	buffer.stream();

	return 0;
}

