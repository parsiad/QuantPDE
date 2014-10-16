////////////////////////////////////////////////////////////////////////////////
// guaranteed_minimum_withdrawal_benefit.cpp
// -----------------------------------------
//
// Computes the price of a Guaranteed Minimum Withdrawal Benefit (GMWB) using an
// implicit, impulse control formulation.
//
// Author: Parsiad Azimzadeh
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>
#include <QuantPDE/Modules/Lambdas>
#include <QuantPDE/Modules/Operators>

////////////////////////////////////////////////////////////////////////////////

#include <algorithm> // max, min
#include <iostream>  // cout
#include <numeric>   // accumulate
#include <tuple>     // get

////////////////////////////////////////////////////////////////////////////////

using namespace QuantPDE;
using namespace QuantPDE::Modules;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

int main() {

	int N = 32; // Initial number of timesteps

	Real T = 10.;
	Real r = .05;
	Real v = .2;

	Real alpha = 0.013886; // Hedging fee

	Real G = 10.; // Contract rate
	Real kappa = 0.1; // Penalty rate

	int E = 10; // Number of events

	int refinement = 5;

	////////////////////////////////////////////////////////////////////////
	// Solution grid
	////////////////////////////////////////////////////////////////////////

	RectilinearGrid2 grid(
		Axis {
			0., 5., 10., 15., 20., 25.,
			30., 35., 40., 45.,
			50., 55., 60., 65., 70., 72.5, 75., 77.5, 80., 82., 84.,
			86., 88., 90., 91., 92., 93., 94., 95.,
			96., 97., 98., 99., 100.,
			101., 102., 103., 104., 105., 106.,
			107., 108., 109., 110., 112., 114.,
			116., 118., 120., 123., 126.,
			130., 135., 140., 145., 150., 160., 175., 200., 225.,
			250., 300., 500., 750., 1000.
		},
		Axis::range(0., 2., 100.)
	);

	int pow2l  = 1; // 2^l
	for(int l = 0; l < refinement; ++l) {

		////////////////////////////////////////////////////////////////
		// Iteration tree
		////////////////////////////////////////////////////////////////

		ReverseConstantStepper stepper(
			0.,              // Initial time
			T,               // Expiry time
			T / (N * pow2l)  // Timestep size
		);

		////////////////////////////////////////////////////////////////
		// Linear system tree
		////////////////////////////////////////////////////////////////

		BlackScholes<2, 0> bs(grid, r, v, alpha);
		ReverseRannacher discretization(grid, bs);
		discretization.setIteration(stepper);

		////////////////////////////////////////////////////////////////
		// Exercise events
		////////////////////////////////////////////////////////////////

		int e = E * pow2l;

		auto withdrawal = [=] (const Interpolant2 &V, Real S, Real W) {
			Real best = V(S, W);

			// You have no money
			if(W < epsilon) {
				return best;
			}

			// Contract withdrawal amount
			const Real Gdt = G * T / e;

			const int partitionSize = 10 * pow2l;

			#if 0
			for(int i = 1; i <= partitionSize; ++i) {
				const Real lambdaW = W * i / partitionSize;
				const Real interp = V(
					max(S - lambdaW, 0.),
					W - lambdaW
				);

				Real cashflow = lambdaW;
				if(lambdaW > Gdt) {
					cashflow -= kappa * (lambdaW - Gdt);
				}

				const Real newValue = interp + cashflow;

				if(newValue > best) {
					best = newValue;
				}
			}
			#endif

			// Nonpenalty
			const Real beta = min(W, Gdt);
			for(int i = 1; i <= partitionSize; ++i) {
				const Real lambdaW = beta * i / partitionSize;
				const Real interp = V(
					max(S - lambdaW, 0.),
					W - lambdaW
				);
				const Real cashflow = lambdaW;
				const Real newValue = interp + cashflow;
				if(newValue > best) {
					best = newValue;
				}
			}

			if(W > Gdt) {
				for(int i = 1; i <= partitionSize; ++i) {
					const Real lambdaW = Gdt + (W - Gdt) * i
							/ partitionSize;
					const Real interp = V(
						max(S - lambdaW, 0.),
						W - lambdaW
					);
					const Real cashflow = lambdaW - kappa
							* (lambdaW - Gdt);
					const Real newValue = interp + cashflow;
					if(newValue > best) {
						best = newValue;
					}
				}
			}

			return best;
		};

		for(int m = 0; m < e; ++m) {
			stepper.add(
				// Time at which the event takes place
				T / e * m,

				withdrawal,

				// Spatial grid to interpolate on
				grid
			);
		}

		////////////////////////////////////////////////////////////////
		// Payoff
		////////////////////////////////////////////////////////////////

		Function2 payoff = [=] (Real S, Real W) {
			return max(S, (1 - kappa) * W);
		};

		////////////////////////////////////////////////////////////////
		// Running
		////////////////////////////////////////////////////////////////

		BiCGSTABSolver solver;

		auto V = stepper.solve(
			grid,    // Domain
			payoff,  // Initial condition
			discretization, // Root of linear system tree
			solver   // Linear system solver
		);

		////////////////////////////////////////////////////////////////
		// Print solution
		////////////////////////////////////////////////////////////////

		RectilinearGrid2 printGrid(
			Axis::range(0., 25., 200.),
			Axis { 100. }
		);
		cout << accessor( printGrid, V ) << endl;

		////////////////////////////////////////////////////////////////
		// Refine Solution grid
		////////////////////////////////////////////////////////////////

		grid.refine( RectilinearGrid2::NewTickBetweenEachPair() );
	}

	return 0;
}

