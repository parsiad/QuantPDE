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
"jump_diffusion [OPTIONS]" << endl << endl <<
"Computes the price of a European put with jump-diffusion driven by a Poisson " << endl <<
"process. The jump amplitude is assumed to be lognormally distributed." << endl <<
"unequal borrowing/lending rates." << endl <<
endl <<
"-d REAL" << endl <<
endl <<
"    Sets the dividend rate (default is 0.)." << endl <<
endl <<
"-g REAL" << endl <<
endl <<
"    Sets the jump amplitude standard deviation (default is 0.42)." << endl <<
endl <<
"-K REAL" << endl <<
endl <<
"    Sets the strike price (default is 100.)." << endl <<
endl <<
"-l NONNEGATIVE_REAL" << endl <<
endl <<
"    Sets the mean arrival time for the jump process (default is 0.05)." << endl <<
endl <<
"-m REAL" << endl <<
endl <<
"    Sets the mean jump amplitude (default is -0.8)." << endl <<
endl <<
"-N POSITIVE_REAL" << endl <<
endl <<
"    Sets the number of timesteps (default is 100)." << endl <<
endl <<
"-r REAL" << endl <<
endl <<
"    Sets the interest rate (default is 0.05)." << endl <<
endl <<
"-R NONNEGATIVE_INTEGER" << endl <<
endl <<
"    Controls the coarseness of the grid, with 0 being coarsest (default is 0)." << endl <<
endl <<
"-T POSITIVE_REAL" << endl <<
endl <<
"    Sets the expiry time (default is 1.)." << endl <<
endl <<
"-v REAL" << endl <<
endl <<
"    Sets the volatility (default is 0.3)." << endl << endl;
}

int main(int argc, char **argv) {

	Real K  = 100.; // Strike
	Real T  = 1.;   // Expiry
	Real r  = .05;  // Interest
	Real v  = .3;   // Volatility
	Real q  = 0.;   // Dividend rate

	Real l = 0.05;  // Mean jump arrival time
	Real m = -.8;   // Mean jump amplitude
	Real g = .42;   // Jump amplitude standard deviation

	int N = 100;    // Number of timesteps
	int R = 0;      // Level of refinement

	// Setting options with getopt
	{ char c;
	while((c = getopt(argc, argv, "d:g:hK:m:N:r:R:T:v:")) != -1) {
		switch(c) {
			case 'd':
				q = atof(optarg);
				break;
			case 'g':
				g = atof(optarg);
				break;
			case 'h':
				help();
				return 0;
			case 'K':
				K = atof(optarg);
				break;
			case 'l':
				l = atof(optarg);
				if(l < 0) {
					cerr <<
"error: the mean arrival time must be nonnegative" << endl;
					return 1;
				}
				break;
			case 'm':
				m = atof(optarg);
				break;
			case 'N':
				N = atoi(optarg);
				if(N <= 0) {
					cerr <<
"error: the number of steps must be positive" << endl;
					return 1;
				}
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
		grid.refine( RectilinearGrid1::NewTickBetweenEachPair() );
	}

	////////////////////////////////////////////////////////////////////////
	// Payoff
	// ------
	//
	// Payoffs are lambda functions. The following is just a macro that
	// expands to
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
		0., // Initial time
		T,  // Expiry time
		T/N // Timestep size
	);

	////////////////////////////////////////////////////////////////////////
	// Linear system tree
	// ------------------
	//
	// Makes the linear system to solve at each iteration
	////////////////////////////////////////////////////////////////////////

	BlackScholesJumpDiffusion bs(
		grid,

		r, // Interest
		v, // Volatility
		q, // Dividend rate

		l, // Mean arrival time (once every ten years)
		lognormal(m, g) // Log-normal probability density
	);

	bs.setIteration(stepper);

	ReverseLinearBDFTwo bdf2(grid, bs);
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
