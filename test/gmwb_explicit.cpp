////////////////////////////////////////////////////////////////////////////////
// gmwb_explicit.cpp
// -----------------
//
// Computes the price of a Guaranteed Minimum Withdrawal Benefit (GMWB) using
// an explicit formulation.
//
// Author: Parsiad Azimzadeh
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>
#include <QuantPDE/Modules/Lambdas>
#include <QuantPDE/Modules/Operators>

#include "gmwb.hpp"

////////////////////////////////////////////////////////////////////////////////

#include <algorithm> // max, min
#include <iostream>  // cout
#include <memory>    // unique_ptr
#include <numeric>   // accumulate
#include <tuple>     // get

////////////////////////////////////////////////////////////////////////////////

using namespace QuantPDE;
using namespace QuantPDE::Modules;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

int main() {

	int N = 32; // Initial number of timesteps

	Real T = 10.; //14.28;
	Real r = .05;
	Real v = .2;

	Real alpha = 0.01389; //0.036; // Hedging fee

	Real G = 10.; //7.; // Contract rate
	Real kappa = 0.1; // Penalty rate

	int refinement = 5;

	bool impulse = false;

	int partitionSize = 10;

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

	for(
		int l = 0;
		l < refinement;
		++l, N *= 2, partitionSize *= 2
	) {

		////////////////////////////////////////////////////////////////
		// Iteration tree
		////////////////////////////////////////////////////////////////

		ReverseConstantStepper stepper(
			0.,    // Initial time
			T,     // Expiry time
			T / N  // Timestep size
		);

		// Tolerance iteration
		unique_ptr<ToleranceIteration> tolerance;
		if(impulse) {
			tolerance = unique_ptr<ToleranceIteration>(
					new ToleranceIteration());
			stepper.setInnerIteration(*tolerance);
		}

		////////////////////////////////////////////////////////////////
		// Linear system tree
		////////////////////////////////////////////////////////////////

		BlackScholes<2, 0> bs(grid, r, v, alpha);
		ReverseLinearBDFOne discretization(grid, bs);
		discretization.setIteration(stepper);

		/*
		unique_ptr<PenaltyMethod> penalty;
		if(impulse) {
			RectilinearGrid1 impulseControls( Axis::range(
							1. / partitionSize,
							1. / partitionSize,
							1.
			) );

			// Impulse withdrawal
			ImpulseWithdrawal impulseWithdrawal(grid, kappa);

			// Impulse withdrawal policy iteration
			MinPolicyIteration2_1 impulsePolicy(
				grid,
				impulseControls,
				impulseWithdrawal
			);
			impulsePolicy.setIteration(tolerance);

			// Penalty method
			penalty = unique_ptr<PenaltyMethod>(
				new PenaltyMethod(
					grid,
					discretization,
					payoff
				)
			);

			penalty->setIteration(*tolerance);
		}
		*/

		////////////////////////////////////////////////////////////////
		// Exercise events
		////////////////////////////////////////////////////////////////

		auto withdrawal = [=] (const Interpolant2 &V, Real S, Real W) {
			Real best = V(S, W);

			// Contract withdrawal amount
			const Real Gdt = G * T / N;

/*
#if   defined(GMWB_CONTRACT_WITHDRAWAL)
			// Constant withdrawal
			const Real gamma = min(W, Gdt);
			const Real interp = V(
				max(S - gamma, 0.),
				W - gamma
			);
			const Real cashflow = gamma;

			best = interp + cashflow;
#elif defined(GMWB_SURRENDER)
			const Real gamma = W;

			const Real interp = V(
				max(S - gamma, 0.),
				W - gamma
			);
			const Real cashflow = gamma  - kappa
					* max(gamma - Gdt, 0.);

			best = interp + cashflow;
#else
*/

			// Nonpenalty
			const Real beta = min(W, Gdt);
			for(int i = 1; i <= partitionSize; ++i) {
				const Real gamma = beta * i / partitionSize;
				const Real interp = V(
					max(S - gamma, 0.),
					W - gamma
				);
				const Real cashflow = gamma;
				const Real newValue = interp + cashflow;
				if(newValue > best) {
					best = newValue;
				}
			}

			// Penalty
			if(!impulse && W > Gdt) {
				for(int i = 1; i <= partitionSize; ++i) {
					const Real gamma = Gdt + (W - Gdt) * i
							/ partitionSize;
					const Real interp = V(
						max(S - gamma, 0.),
						W - gamma
					);
					const Real cashflow = gamma - kappa
							* (gamma - Gdt);
					const Real newValue = interp + cashflow;
					if(newValue > best) {
						best = newValue;
					}
				}
			}
//#endif

			return best;
		};

		for(int m = 0; m < N; ++m) {
			stepper.add(
				// Time at which the event takes place
				T / N * m,

				withdrawal,

				// Spatial grid to interpolate on
				grid

				// Null-event test
				//std::unique_ptr<EventBase>(new NullEvent)
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

