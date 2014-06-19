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
#include <QuantPDE/Modules/Payoffs>

// TODO: Change these includes; shouldn't include src directory explicitly
#include <QuantPDE/src/Modules/BlackScholes.hpp>

using namespace QuantPDE;
using namespace QuantPDE::Modules;

////////////////////////////////////////////////////////////////////////////////

#include <iostream>  // cout, cerr
#include <unistd.h>  // getopt

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/**
 * Prints help to stderr.
 */
void help() {
	cerr <<
"unequal_borrowing_lending_rates [OPTIONS]" << endl << endl <<
"Prices a long/short position straddle under the Black-Scholes model assuming" << endl <<
"unequal borrowing/lending rates." << endl <<
endl <<
"-b REAL" << endl <<
endl <<
"    sets the borrowing interest rate (default is 0.05)" << endl <<
endl <<
"-d REAL" << endl <<
endl <<
"    sets the dividend rate (default is 0.)" << endl <<
endl <<
"-K REAL" << endl <<
endl <<
"    sets the strike price (default is 100.)" << endl <<
endl <<
"-l REAL" << endl <<
endl <<
"    sets lending interest rate (default is 0.03)" << endl <<
endl <<
"-L" << endl <<
endl <<
"    long position (default is short)" << endl <<
endl <<
"-N POSITIVE_INTEGER" << endl <<
endl <<
"    sets the number of steps to take in time (default is 100)" << endl <<
endl <<
"-R NONNEGATIVE_INTEGER" << endl <<
endl <<
"    controls the coarseness of the grid, with 0 being coarsest (default is 0)" << endl <<
endl <<
"-T POSITIVE_REAL" << endl <<
endl <<
"    sets the expiry time (default is 1.)" << endl <<
endl <<
"-v REAL" << endl <<
endl <<
"    sets the volatility" << endl << endl;
}

int main(int argc, char **argv) {

	Real K     = 100.;
	Real T     = 1.;
	Real r_l   = .03;
	Real r_b   = .05;
	Real vol   = .3;
	Real div   = 0.;
	int N      = 100;
	int R      = 0;
	bool L     = false;

	// Setting options with getopt
	{ char c;
	while((c = getopt(argc, argv, "b:d:hK:l:LN:R:T:v:")) != -1) {
		switch(c) {
			case 'b':
				r_b = atof(optarg);
				break;
			case 'd':
				div = atof(optarg);
				break;
			case 'h':
				help();
				return 0;
			case 'K':
				K = atof(optarg);
				break;
			case 'l':
				r_l = atof(optarg);
				break;
			case 'L':
				L = true;
				break;
			case 'N':
				N = atoi(optarg);
				if(N <= 0) {
					cerr <<
"error: the number of steps must be positive" << endl;
					return 1;
				}
				break;
			case 'R':
				R = atoi(optarg);
				if(R < 0) {
					cerr <<
"error: the maximum level of refinement must be nonnegative" << endl;
					return 1;
				}
				break;
			case 'T':
				if((T = atof(optarg)) <= 0.) {
					cerr <<
"error: expiry time must be positive" << endl;
					return 1;
				}
				break;
			case 'v':
				vol = atof(optarg);
				break;
			case ':':
			case '?':
				cerr << endl;
				help();
				return 1;
		}
	} }

	////////////////////////////////////////////////////////////////////////

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

	// Refine grid
	for(int r = 0; r < R; r++) {
		grid.refine( RectilinearGrid1::NewTickBetweenEachPair() );
	}

	// Control can be two interest rates: lending or borrowing
	RectilinearGrid1 controls( Axis { r_l, r_b } );

	// Payoff is a straddle
	auto payoff = QUANT_PDE_MODULES_PAYOFFS_STRADDLE_FIXED_STRIKE( K );
	// Whenever possible, QuantPDE uses lambda functions to specify
	// functions.
	// For example, payoffs (in general: initial conditions), are always
	// specified as lambda functions.
	// The above is just a macro that expands to
	//
	// auto payoff = [K] (Real S) {
	// 	return S < K ? K - S : S - K;
	// };

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
		div  // Dividend rate
	);

	// Policy iteration
	//
	// The notation [N][_M] at the end of class names is used when the
	// problem is N-dimensional (in space) with an M-dimensional control.
	// In the following, we instantiate a policy iteration class meant for
	// one-dimensional (in space) problems with a one-dimensional control.
	//
	// Pick min/max policy iteration depending on whether we are considering
	// the long or short position problem.
	LinearSystemIteration *policy = L
		? (LinearSystemIteration*)
				new MinPolicyIteration1_1(grid, controls, bs)
		: (LinearSystemIteration*)
				new MaxPolicyIteration1_1(grid, controls, bs)
	;
	policy->setIteration(tolerance); // Associate with k-iteration

	// BDF2 (timestepping)
	ReverseLinearBDFTwo bdf2(grid, *policy);
	bdf2.setIteration(stepper); // Associate with n-iteration

	////////////////////////////////////////////////////////////////////////

	// Initial solution (maps the payoff function to a vector)
	Vector initial = (PointwiseMap1(grid))(payoff);

	// Linear system solver
	BiCGSTABSolver solver;

	// Get the solution vector
	Vector solution = stepper.iterateUntilDone(
		initial, // Initial iterand
		bdf2,    // Root of linear system tree
		solver   // Linear system solver
	);

	////////////////////////////////////////////////////////////////////////

	// Print solution
	cout << grid.accessor(solution);

	// Cleanup
	delete policy;

	return 0;

}

