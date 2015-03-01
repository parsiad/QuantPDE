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
#include <iomanip>   // setw
#include <iostream>  // cout, cerr
#include <cmath>     // exp

using namespace std;

///////////////////////////////////////////////////////////////////////////////

Real S_0             = 100.;  // Initial stock value
Real L_0             = 120.;  // Representative value of L
                              // Pick \hat{L} = L_0 (initial value of loan)

Real r               = 0.04;  // Interest rate
Real s               = 0.;    // Spread
Real sigma           = 0.2;   // Volatility

Real lambda          = 0.05;  // Jump arrival rate
Real mu_xi           = -.8;   // Mean jump amplitude
Real sigma_xi        = .42;   // Jump amplitude standard deviation

Real beta_lo         = 0.7;   // Low trigger
Real beta_hi         = 0.9;   // High trigger

Real T               = 1.;    // Expiry

int borrowerEvents   = 252;   // Daily chance to walk-away or repay the loan
int bankEvents       = 252;   // Daily request to top-up or chance to liquidate
int interestPayments = 12;    // Monthly interest payments

int N_0              = 100;   // Initial number of steps
int maxRefinement    = 8;     // Maximum number of times to refine

bool newton          = false; // Newton iteration to find optimal spread

///////////////////////////////////////////////////////////////////////////////

std::vector<Real> interestPaymentDates; // Sorted (ascending) vector of interest
                                        // payment dates

// Payoff for fixed \hat{L}
Real payoff(Real S) {
	return min(S, L_0);
}

// Accrued interest
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

	// Previous interest payment time
	const Real t_previous = interestPaymentDates[lo];

	const Real A_t = exp( (r + s) * (t - t_previous) );
	return A_t;
}

///////////////////////////////////////////////////////////////////////////////

constexpr int td = 20;

void printHeaders() {
	cout
		<< setw(td) << "Nodes"           << "\t"
		<< setw(td) << "Time Steps"      << "\t"
		<< setw(td) << "Value"           << "\t"
		<< setw(td) << "Change"          << "\t"
	;

	if(!newton) {
		cout << setw(td) << "Ratio";
	}

	cout << endl;
}

int main(int argc, char **argv) {
	cout.precision(12); // Precision
	printHeaders();

	// Populate vector of interest payment dates
	{
		const Real dt = T / interestPayments;

		interestPaymentDates.reserve(interestPayments);
		for(int i = 0; i <= interestPayments; ++i) {
			interestPaymentDates.push_back(dt * i);
		}
	} // Checked 2015-02-28

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

	Real previousValue = nan(""), previousChange = nan("");
	for(
		int N = N_0, refinement = 0;
		refinement < maxRefinement;
		++refinement, N *= 2
	) { // <RefinementLoop>

	// Constant step-size
	ReverseConstantStepper stepper(
		0., // Initial time
		T,  // Expiry time
		T/N // Timestep size
	);

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
				const Real A_t_i = A(t_i - epsilon);
				return min( U(S), min(L_0 * A_t_i, S) );
			},
			grid
		);
	}

	} // </BorrowerEvent>

	////////////////////////////////////////////////////////////////////////

	{ // <InterestPayment>

	// Add interest payment events (skip initial time)
	for(
		auto it = interestPaymentDates.begin() + 1;
		it != interestPaymentDates.end();
		++it
	) {
		// Time at which interest payment takes place
		const Real t_i = *it;

		stepper.add(
			t_i,
			[=] (const Interpolant1 &U, Real S) {
				return U(S) + L_0 * (A(t_i) - 1.);
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

				// Value-to-loan ratio
				const Real Rinv = S / L_0;

				// R is below low trigger
				if(Rinv <= beta_lo) {
					const Real A_t_i = A(t_i + epsilon);

					// Option to liquidate
					const Real liquidate = min(
						L_0 * A_t_i,
						S
					);
					if(liquidate > best) {
						best = liquidate;
					}

				// R is between low and high triggers
				} else if(Rinv <= beta_hi) {
					const Real A_t_i = A(t_i + epsilon);

					// Top-up with shares
					const Real X = U(S_0);
					if(X > best) {
						best = X;
					}

					// Top-up with cash
					// Note: similarity reduction here only
					//       works if \hat{L} = L_0!
					const Real S_ret = S / S_0;
					const Real Y = S_ret * U(S_0) + L_0
							* (1 - S_ret) * A_t_i;
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

	////////////////////////////////////////////////////////////////////////

	Real value = U(S_0);

	Real
		change = value - previousValue,
		ratio = previousChange / change
	;

	cout
		<< setw(td) << grid.size() << "\t"
		<< setw(td) << N           << "\t"
		<< setw(td) << value       << "\t"
		<< setw(td) << change      << "\t"
		<< setw(td) << ratio       << "\t"
		<< endl
	;

	previousChange = change;
	previousValue = value;

	////////////////////////////////////////////////////////////////////////

	// Refine the grid
	grid = grid.refined();

	} // </RefinementLoop>
}
