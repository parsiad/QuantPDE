////////////////////////////////////////////////////////////////////////////////
// bermudan_put.cpp
// ----------------
//
// Computes the price of a Bermudan put with exercise occuring yearly.
//
// Author: Parsiad Azimzadeh
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>
#include <QuantPDE/Modules/Payoffs>

// TODO: Change these includes; shouldn't include src directory explicitly
#include <QuantPDE/src/Modules/BlackScholes.hpp>

using namespace QuantPDE;
using namespace QuantPDE::Modules;

///////////////////////////////////////////////////////////////////////////////

#include <iostream>  // cout, cerr
//#include <iomanip>   // setw
//#include <unistd.h>  // getopt

using namespace std;

///////////////////////////////////////////////////////////////////////////////

int main() {

	const Real K = 100.;
	const Real T = 1.;
	const Real r = 0.04;
	const Real v = 0.2;
	const Real q = 0.;

	const unsigned M = 10; // Must be > 0
	const unsigned N = 25;

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
	grid.refine( RectilinearGrid1::NewTickBetweenEachPair() );

	auto payoff = QUANT_PDE_MODULES_PAYOFFS_PUT_FIXED_STRIKE( K );

	ReverseConstantStepper::Factory factory(N);
	ReverseEventIteration1 stepper(
		0., // Initial time
		T,  // Expiry time
		factory
	);

	////////////////////////////////////////////////////////////////////////
	// Exercise events
	////////////////////////////////////////////////////////////////////////

	for(unsigned m = 0; m < M; m++) {
		stepper.add(
			// Time at which the event takes place
			T / M * m,

			// Take the maximum of the continuation and exercise
			// values
			[K] (const Interpolant1 &V, Real S) {
				return max( V(S) , K - S );
			},

			// Spatial grid to interpolate on
			grid
		);
	}

	////////////////////////////////////////////////////////////////////////

	BlackScholes bs(
		grid,

		r, // Interest
		v, // Volatility
		q  // Dividend rate
	);

	ReverseLinearBDFTwo bdf2(grid, bs);
	bdf2.setIteration(stepper);

	Vector initial = (PointwiseMap1(grid))(payoff);

	BiCGSTABSolver solver;

	Vector solutionOnGrid = stepper.iterateUntilDone(
		initial, // Initial iterand
		bdf2,    // Root of linear system tree
		solver   // Linear system solver
	);

	cout << grid.accessor(solutionOnGrid);

	return 0;

}
