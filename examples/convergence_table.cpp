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
#include <QuantPDE/Modules/Configuration>
#include <QuantPDE/Modules/Lambdas>
#include <QuantPDE/Modules/Operators>

using namespace QuantPDE;
using namespace QuantPDE::Modules;

///////////////////////////////////////////////////////////////////////////////

#include <memory>   // unique_ptr
#include <numeric>  // accumulate

using namespace std;

///////////////////////////////////////////////////////////////////////////////

Configuration configuration;
Real T, r, vol, divs, S_0, K;
int N, M;
bool call, digital, american, variable;
RectilinearGrid1 *grid;

std::vector<Real> run(int k) {

	// 2^k or 2^(2k)
	Real factor = 1;
	for(int i = 0; i < k; ++i) {
		if(american) { factor *= 4; }
		else         { factor *= 2; }
	}

	Function1 payoff = digital ?
		  ( call ? digitalCallPayoff(K) : digitalPutPayoff(K) )
		: (call ? callPayoff(K) : putPayoff(K))
	;

	auto refined_grid = grid->refined(k);

	unsigned outer;
	Real inner = nan("");
	Real value;

	// Black-Scholes operator (L in V_t = LV)
	BlackScholes1 bs(
		refined_grid,
		r, vol, divs
	);

	// Timestepping method
	unique_ptr<Iteration> stepper(variable
		? (Iteration *) new ReverseVariableStepper(
			0.,              // startTime
			T,               // endTime
			T / N / factor,  // dt
			1. / factor      // target
		)
		: (Iteration *) new ReverseConstantStepper(
			0.,            // startTime
			T,             // endTime
			T / N / factor // dt
		)
	);

	// Time discretization method
	ReverseRannacher discretization(refined_grid, bs);
	discretization.setIteration(*stepper);

	// American-specific components; penalty method or not?
	IterationNode *root;
	unique_ptr<ToleranceIteration> tolerance;
	unique_ptr<PenaltyMethodDifference1> penalty;
	if(american) {
		// American case
		penalty = unique_ptr<PenaltyMethodDifference1>(
			new PenaltyMethodDifference1(
				refined_grid,
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
	BiCGSTABSolver solver;

	// Compute solution
	auto solution = stepper->solve(
		refined_grid,
		payoff,
		*root,
		solver
	);

	// Outer iterations
	outer = stepper->iterations()[0];

	// Average number of inner iterations
	if(american) {
		auto its = tolerance->iterations();
		inner = accumulate(its.begin(), its.end(), 0.)
				/ its.size();
	}

	// Solution at S = stock (default is 100.)
	// Linear interpolation is used to get the value not on the grid
	value = solution(S_0);

	return
		{
			(Real) refined_grid.size(),
			(Real) outer,
			(Real) inner,
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
	r = getReal(configuration, "interest_rate", .04);
	vol = getReal(configuration, "volatility", .2);
	divs = getReal(configuration, "dividend_rate", 0.);
	S_0 = getReal(configuration, "asset_price", 100.);
	K = getReal(configuration, "strike_price", 100.);
	N = getInt(configuration, "initial_number_of_timesteps", 12);
	call = getBool(configuration, "is_call", true);
	digital = getBool(configuration, "is_digital", false);
	american = getBool(configuration, "is_american", false);
	variable = getBool(configuration, "use_variable_timestepping", false);
	RectilinearGrid1 defGrid( (S_0 * Axis::special) + (K * Axis::special));
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

