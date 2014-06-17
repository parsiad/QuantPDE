#include <QuantPDE/Core>
#include <QuantPDE/Modules/Payoffs>

// TODO: Change these includes; shouldn't include src directory explicitly
#include <QuantPDE/src/Modules/BlackScholes.hpp>

using namespace QuantPDE;
using namespace QuantPDE::Modules;

///////////////////////////////////////////////////////////////////////////////

#include <iostream>  // cout, cerr

using namespace std;

///////////////////////////////////////////////////////////////////////////////

int main() {

	const Real
		K = 100.,
		T = 1.,
		r_l = .03,
		r_b = .05,
		vol = .3
	;

	const unsigned
		N = 100
	;

	// Initial discretization
	RectilinearGrid1 grid(
		Axis {
			0., 10., 20., 30., 40., 50., 60., 70.,
			75., 80.,
			84., 88., 92.,
			94., 96., 98., 100., 102., 104., 106., 108., 110.,
			114., 118.,
			123.,
			130., 140., 150.,
			175.,
			225.,
			300.,
			750.,
			2000.,
			10000.
		}
	);

	// Control can be two interest rates: lending or borrowing
	RectilinearGrid1 controls( Axis { r_l, r_b } );

	// Payoff is a straddle
	auto payoff = QUANT_PDE_MODULES_PAYOFFS_STRADDLE_FIXED_STRIKE( K );

	// Iteration tree
	ReverseConstantStepper stepper(
		0.,
		T,  // Expiry time
		N   // Number of steps
	);
	ToleranceIteration tolerance;
	stepper.setInnerIteration(tolerance);

	// Linear system tree
	BlackScholes bs(
		grid,

		// Interest rate (passed as a control)
		QUANT_PDE_CONTROL1(grid),

		vol, // Volatility
		0.   // Dividend rate
	);
	MaxPolicyIteration1_1 policy(grid, controls, bs);
	policy.setIteration(tolerance);
	ReverseLinearBDFTwo bdf2(grid, policy);
	bdf2.setIteration(stepper);

	// Initial solution
	Vector initial = (PointwiseMap1(grid))(payoff);

	// Linear system solver
	BiCGSTABSolver solver;

	// Get solution
	Vector solutionVector = stepper.iterateUntilDone(
		initial,
		bdf2,
		solver
	);

	// Print solution
	cout << grid.accessor(solutionVector);

	return 0;

}

