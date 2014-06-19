////////////////////////////////////////////////////////////////////////////////
// unequal_borrowing_lending_rates.cpp
// -----------------------------------
//
// Prices a long-position straddle option under the Black-Scholes model assuming
// unequal borrowing/lending rates. The pricing problem is given by:
//
// V_t = \sup_{r \in \{ r_l, r_b \}} (-\sigma^2 S^2 V_{SS} / 2 + r( V - S V_S ))
// V(0, S) = max(S - K, K - S)
//
// Author: Parsiad Azimzadeh
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>
#include <QuantPDE/Modules/Payoffs>

// TODO: Change these includes; shouldn't include src directory explicitly
#include <QuantPDE/src/Modules/BlackScholes.hpp>

using namespace QuantPDE;
using namespace QuantPDE::Modules;

////////////////////////////////////////////////////////////////////////////////

#include <iostream>  // cout, cerr

using namespace std;

////////////////////////////////////////////////////////////////////////////////

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

	////////////////////////////////////////////////////////////////////////
	// Iteration tree
	//
	// Sets up the loop structure:
	// for(int n = 0; n < N; n++) {
	// 	for(int k = 0; ; k++) {
	// 		// do stuff
	// 		if(error < tolerance) break;
	// 	}
	// }
	////////////////////////////////////////////////////////////////////////

	ReverseConstantStepper stepper(
		0., // Initial time
		T,  // Expiry time
		N   // Number of steps
	);
	ToleranceIteration tolerance;
	stepper.setInnerIteration(tolerance);

	////////////////////////////////////////////////////////////////////////
	// Linear system tree
	//
	// Makes the linear system to solve at each iteration
	////////////////////////////////////////////////////////////////////////

	BlackScholes bs(
		grid,

		// Interest rate (passed as a control)
		Control1::make(grid),

		vol, // Volatility
		0.   // Dividend rate
	);

	// Policy iteration
	// The notation [N][_M] at the end of class names is used when the
	// problem is N-dimensional (in space) with an M-dimensional control.
	// In the following, we instantiate a policy iteration class meant for
	// one-dimensional (in space) problems with a one-dimensional control.
	MaxPolicyIteration1_1 policy(grid, controls, bs);
	policy.setIteration(tolerance); // Associate with k-iteration

	// BDF2 (timestepping)
	ReverseLinearBDFTwo bdf2(grid, policy);
	bdf2.setIteration(stepper); // Associate with n-iteration

	////////////////////////////////////////////////////////////////////////

	// Initial solution (maps the payoff function to a vector)
	Vector initial = (PointwiseMap1(grid))(payoff);

	// Linear system solver
	BiCGSTABSolver solver;

	// Get the solution vector (not a function)
	Vector solution = stepper.iterateUntilDone(
		initial, // Initial iterand
		bdf2,    // Root of linear system tree
		solver   // Linear system solver
	);

	////////////////////////////////////////////////////////////////////////

	// Print solution
	cout << grid.accessor(solution);

	return 0;

}

