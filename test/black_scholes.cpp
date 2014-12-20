////////////////////////////////////////////////////////////////////////////////
// black_scholes.cpp
// -----------------
//
// Computes a numerical solution to the Black-Scholes equation with constant
// coefficients.
//
// Author: Parsiad Azimzadeh
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>
#include <QuantPDE/Modules/Lambdas>
#include <QuantPDE/Modules/Operators>

using namespace QuantPDE;
using namespace QuantPDE::Modules;

///////////////////////////////////////////////////////////////////////////////

#include <iostream> // cout, cerr
#include <unistd.h> // getopt

using namespace std;

///////////////////////////////////////////////////////////////////////////////

/**
 * Prints help to stderr.
 */
void help() {
	cerr <<
"unequal_borrowing_lending_rates [OPTIONS]" << endl << endl <<
"Prices a long/short position put under the Black-Scholes model assuming" << endl <<
"unequal borrowing/lending rates." << endl <<
endl <<
"-d REAL" << endl <<
endl <<
"    Sets the dividend rate (default is 0.)." << endl <<
endl <<
"-K REAL" << endl <<
endl <<
"    Sets the strike price (default is 100.)." << endl <<
endl <<
"-r REAL" << endl <<
endl <<
"    Sets the interest rate (default is 0.04)." << endl <<
endl <<
"-R NONNEGATIVE_INTEGER" << endl <<
endl <<
"    Controls the coarseness of the grid, with 0 being coarsest (default is 0)." << endl <<
endl <<
"-t POSITIVE_REAL" << endl <<
endl <<
"    Sets the timestep size (default is 0.01)." << endl <<
endl <<
"-T POSITIVE_REAL" << endl <<
endl <<
"    Sets the expiry time (default is 1.)." << endl <<
endl <<
"-v REAL" << endl <<
endl <<
"    Sets the volatility (default is 0.2)." << endl << endl;
}

int main(int argc, char **argv) {

	Real K  = 100.;
	Real T  = 1.;
	Real r  = .04;
	Real v  = .2;
	Real q  = 0.;
	Real dt = .01;

	int R = 0;

	// Setting options with getopt
	{ char c;
	while((c = getopt(argc, argv, "d:hK:r:R:t:T:v:")) != -1) {
		switch(c) {
			case 'd':
				q = atof(optarg);
				break;
			case 'h':
				help();
				return 0;
			case 'K':
				K = atof(optarg);
				break;
			case 'r':
				r = atof(optarg);
				break;
			case 'R':
				R = atoi(optarg);
				if(R < 0) {
					cerr <<
"error: the maximum level of refinement must be nonnegative" << endl;
					return 1;
				}
				break;
			case 't':
				dt = atof(optarg);
				if(dt <= 0.) {
					cerr <<
"error: the timestep size must be positive" << endl;
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
				v = atof(optarg);
				break;
			case ':':
			case '?':
				cerr << endl;
				help();
				return 1;
		}
	} }

	////////////////////////////////////////////////////////////////////////
	// Spatial grid
	////////////////////////////////////////////////////////////////////////

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
	for(int i = 0; i < R; ++i) {
		grid = grid.refined();
	}

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

	auto payoff = putPayoff(K);

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
		0., // Initial time
		T,  // Expiry time
		dt  // Timestep size
	);

	////////////////////////////////////////////////////////////////////////
	// Linear system tree
	// ------------------
	//
	// Makes the linear system to solve at each iteration
	////////////////////////////////////////////////////////////////////////

	BlackScholes1 bs(
		grid,

		r, // Interest
		v, // Volatility
		q  // Dividend rate
	);

	ReverseBDFTwo1 bdf2(grid, bs);
	bdf2.setIteration(stepper);

	////////////////////////////////////////////////////////////////////////
	// Running
	// -------
	//
	// Everything prior to this was setup. Now we run the method.
	////////////////////////////////////////////////////////////////////////

	BiCGSTABSolver solver;

	auto V = stepper.solve(
		grid,   // Domain
		payoff, // Initial condition
		bdf2,   // Root of linear system tree
		solver  // Linear system solver
	);

	////////////////////////////////////////////////////////////////////////
	// Print solution
	////////////////////////////////////////////////////////////////////////

	RectilinearGrid1 printGrid( Axis::range(0., 10., 200.) );
	cout << accessor( printGrid, V );

	return 0;

}
