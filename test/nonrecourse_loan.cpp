////////////////////////////////////////////////////////////////////////////////
// nonrecourse_loan.cpp
// --------------------
//
// Author: Parsiad Azimzadeh
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>
#include <QuantPDE/Modules/Lambdas>
#include <QuantPDE/Modules/Operators>

using namespace QuantPDE;
using namespace QuantPDE::Modules;

///////////////////////////////////////////////////////////////////////////////

#include <algorithm> // min
#include <iostream>  // cout, cerr
#include <cmath>     // exp

using namespace std;

///////////////////////////////////////////////////////////////////////////////

Real S_0             = 100.; // Initial stock value

Real r               = 0.04; // Interest rate
Real s               = 0.;   // Spread
Real sigma           = 0.2;  // Volatility

Real lambda          = 0.05; // Jump arrival rate
Real mu_xi           = -.8;  // Mean jump amplitude
Real sigma_xi        = .42;  // Jump amplitude standard deviation

Real L_hat           = 100.; // Representative value of L
                             // Pick \hat{L} = L_0 (initial value of loan)

Real beta_lo         = 0.5;  // Low trigger
Real beta_hi         = 0.75; // High trigger

Real T               = 1.;   // Expiry
int N                = 100;  // Initial number of steps

int borrowerEvents   = 252;  // Daily chance to walk-away or repay the loan
int bankEvents       = 252;  // Daily request to top-up or chance to liquidate
int interestPayments = 12;   // Monthly interest payments

///////////////////////////////////////////////////////////////////////////////

std::vector<Real> interestPaymentDates; // Sorted (ascending) vector of interest
                                        // payment dates

// Payoff for fixed \hat{L}
Real payoff(Real S) { return min(S, L_hat); }

// Across interest payment dates
Real A(Real t) {
	// Binary search for interest payment date
	int lo = 0;
	int hi = interestPaymentDates.size() - 1;
	while(lo < hi - 1) {
		int mid = (lo + hi) / 2;
		if(t <= interestPaymentDates[mid]) {
			hi = mid;
		} else {
			lo = mid;
		}
	}

	// lo = t_{i-1}
	return exp( (r + s) * (t - lo) );
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
	// Constant step-size
	ReverseConstantStepper stepper(
		0., // Initial time
		T,  // Expiry time
		T/N // Timestep size
	);

	// Hand-picked grid, scaled by S_0
	RectilinearGrid1 grid(
		(S_0 / 100.) * Axis {
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
		}
	);

	// Populate vector of interest payment dates
	{
		const Real dt = T / interestPayments;

		interestPaymentDates.reserve(interestPayments);
		for(int i = 1; i <= interestPayments; ++i) {
			interestPaymentDates.push_back(dt * interestPayments);
		}
	}

	////////////////////////////////////////////////////////////////////////

	// Order of events (forward in time)
	// 1. Borrower: walk-away or repay
	// 2. Interest payment and (TODO) dividends
	// 3. Bank: top-up or liquidate

	////////////////////////////////////////////////////////////////////////

	{ // <BorrowerEvent>

	// Time between borrower events
	const Real dt = T / borrowerEvents;

	for(int i = 1; i <= borrowerEvents; ++i) {
		const Real t_i = dt * i;

		// TODO: If the borrower event occurs after an interest
		//       payment, we need to change A(t_i)

		// TODO: Confirm L^{\star}(t+) = L^{\star}(t-)
		//                             = L^{\star}(t ) ?

		stepper.add(
			t_i,
			[=] (const Interpolant1 &U, Real S) {
				return min(
					U(S),                  // Continue
					min(L_hat * A(t_i), S) // Walk away
				);
			},
			grid
		);
	}

	} // </BorrowerEvent>

	////////////////////////////////////////////////////////////////////////

	{ // <InterestPayment>

	// Time between interest payments
	const Real dt = T / interestPayments;

	// Add interest payment events (skip initial time)
	for(int i = 1; i <= interestPayments; ++i) {
		// Time at which interest payment takes place
		const Real t_i = dt * i;

		stepper.add(
			t_i,
			[=] (const Interpolant1 &U, Real S) {
				return U(S) + L_hat * (A(t_i) - 1.);
			},
			grid
		);
	}

	} // </InterestPayment>

	////////////////////////////////////////////////////////////////////////

	{ // <BankEvent>

	// Time between bank events
	const Real dt = T / bankEvents;

	for(int i = 1; i <= bankEvents; ++i) {
		const Real t_i = dt * i;

		stepper.add(
			t_i,
			[=] (const Interpolant1 &U, Real S) {
				// Continuation value
				Real best = U(S);

				// Loan-to-value ratio
				const Real R = L_hat / S;

				// R is below low trigger
				if(R <= beta_lo) {

					// Option to liquidate
					const Real liquidate = min(
						L_hat * A(t_i),
						S
					);
					if(liquidate > best) {
						best = liquidate;
					}

				// R is between low and high triggers
				} else if(R <= beta_hi) {

					// Top-up with shares
					const Real X = U(S_0);
					if(X > best) {
						best = X;
					}

					// Top-up with cash
					// Note: similarity reduction here only
					//       works if \hat{L} = L_0!
					const Real Y = S / S_0 * U(S_0)
							+ L_hat * (1 - S / S_0)
							* A(t_i);
					if(Y > best) {
						best = Y;
					}

				}

				return best;
			},
			grid
		);
	}

	} // </BankEvent>

	////////////////////////////////////////////////////////////////////////

	// Jump-diffusion operator
	BlackScholesJumpDiffusion1 bs(
		grid,

		r,      // Interest
		sigma,  // Volatility
		0.,     // Dividend rate

		lambda, // Mean arrival time
		lognormal(mu_xi, sigma_xi) // Log-normal probability density
	);
	bs.setIteration(stepper);

	// BDF2
	ReverseBDFTwo1 bdf2(grid, bs);
	bdf2.setIteration(stepper);

	// Linear system solver
	BiCGSTABSolver solver;
	auto U = stepper.solve(
		grid,   // Domain
		payoff, // Initial condition
		bdf2,   // Root of linear system tree
		solver  // Linear system solver
	);

	// Solution for L = L_hat
	RectilinearGrid1 printGrid( Axis::range(0., 10., 200.) );
	cout << accessor( printGrid, U );
}
