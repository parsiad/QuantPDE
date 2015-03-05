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

////////////////////////////////////////////////////////////////////////////////

#include <algorithm> // max, min
#include <chrono>    // duration
#include <iomanip>   // setw
#include <iostream>  // cout, cerr
#include <cmath>     // exp

using namespace std;

////////////////////////////////////////////////////////////////////////////////

Real beta_lo         = 0.8;   // Low trigger
Real beta_hi         = 0.9;   // High trigger

Real S_0             = 125.;  // Initial stock value
Real L_0             = 100.;  // Initial loan value
Real L_hat           = 100.;  // Representative value

Real r               = 0.04;  // Interest rate
Real s               = 0.02;  // Spread
Real sigma           = 0.2;   // Volatility

Real lambda          = 0.1;  // Jump arrival rate
Real mu_xi           = -.8;   // Mean jump amplitude
Real sigma_xi        = .42;   // Jump amplitude standard deviation

Real T               = 1.;    // Expiry

int borrowerEvents   = -1;    // Borrower events (-1 for all times)
int bankEvents       = 4;     // Bank events (-1 for all times)
int interestPayments = 4;     // Interest payments (once a quarter)

int N                = 12;    // Initial number of steps
int maxRefinement    = 10;    // Maximum number of times to refine

bool fairSpread      = false; // Calculate fair spread

////////////////////////////////////////////////////////////////////////////////

std::vector<Real> interestPaymentDates; // Sorted (ascending) vector of interest
                                        // payment dates

/**
 * Payoff from the bank's perspective for representative value of L.
 * @param S The stock value at the expiry.
 * @return The payoff min(S, L_hat).
 */
Real payoff(Real S) { return min(S, L_hat); }

/**
 * Accrued interest at a particular time between coupon dates.
 * @param t The particular time.
 * @return e^((r+s)(t-t_p)) where t_p is the previous coupon date.
 */
Real accruedInterest(Real t) {
	// Binary search for interest payment date
	int lo = 0;
	int hi = interestPaymentDates.size() - 1;

	if(t > interestPaymentDates[hi]) {
		lo = hi;
	} else {
		while(lo < hi - 1) {
			int mid = (lo + hi) / 2;
			if(t <= interestPaymentDates[mid]) {
				hi = mid;
			} else {
				lo = mid;
			}
		}
	}

	// Previous interest payment time
	const Real t_previous = interestPaymentDates[lo];

	const Real A_t = exp( (r + s) * (t - t_previous) );
	return A_t;
}

/**
 * Returns a hand-picked axis for the problem, centred at some parameter.
 * @param C The centre.
 * @return An axis with ticks between 0 and 100C (inclusive).
 */
Axis axis(Real C) {
	return (C / 100.) * Axis {
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
	};
}

// Initial grid
RectilinearGrid1 grid( axis(S_0) );

////////////////////////////////////////////////////////////////////////////////

constexpr int td = 20;

/**
 * Print headers to stdout.
 */
void printHeaders() {
	cout
		<< setw(td) << "Nodes"           << "\t"
		<< setw(td) << "Time Steps"      << "\t"
		<< setw(td) << "Value"           << "\t"
		<< setw(td) << "Change"          << "\t"
	;

	if(!fairSpread) {
		cout << setw(td) << "Ratio";
	}

	cout << setw(td) << "Seconds";

	cout << endl;
}

/**
 * Solve the nonrecourse stock loan problem with fixed spread.
 * @param s Spread.
 */
Real solve(Real s) {
	// Constant step-size
	ReverseConstantStepper stepper(
		0.,   // Initial time
		T,    // Expiry time
		T / N // Timestep size
	);

	////////////////////////////////////////////////////////////////////////

	// Order of events (forward in time)
	// 1. Borrower: walk-away or repay
	// 2. Bank: top-up or liquidate
	// 3. Interest payment

	// Apply events at all times
	if(borrowerEvents < 0) { borrowerEvents = N; }
	if(bankEvents     < 0) { bankEvents     = N; }

	////////////////////////////////////////////////////////////////////////

	{ // <BorrowerEvent>

	// Time between borrower events
	const Real dt = T / borrowerEvents;

	for(int i = 0; i < borrowerEvents; ++i) {
		const Real t_i = dt * i;

		// Accrued interest
		const Real A = accruedInterest(t_i);

		stepper.add(
			t_i,
			[=] (const Interpolant1 &U, Real S) {
				return min( U(S), min(L_hat * A, S) );
			},
			grid
		);
	}

	} // </BorrowerEvent>
	// 2015-03-04: checked; respects monotonicity

	////////////////////////////////////////////////////////////////////////

	{ // <BankEvent>

	// Time between bank events
	const Real dt = T / bankEvents;

	for(int i = 0; i < bankEvents; ++i) {
		const Real t_i = dt * i;

		// Accrued interest
		const Real A = accruedInterest(t_i);

		stepper.add(
			t_i,
			[=] (const Interpolant1 &U, Real S) {
				// Continuation value
				Real best = U(S);

				// Value-to-loan ratio
				const Real R = L_hat / S;

				// R is above high trigger
				if(R > beta_hi) {

					// Option to liquidate
					const Real tmp = min(L_hat * A, S);
					if(tmp > best) {
						best = tmp;
					}

				// R is between low and high triggers
				} else if(R >= beta_lo) {
					const Real R_0 = L_0 / S_0;

					{ // Top-up with shares
						const Real tmp = U(L_hat / R_0);
						if(tmp > best) {
							best = tmp;
						}
					}

					// Top-up with cash
					const Real tmp = S * R_0 / L_hat
							* U(L_hat / R_0)
							+ (L_hat - S * R_0) * A;
					;
					if(tmp > best) {
						best = tmp;
					}
				}

				return best;
			},
			grid
		);
	}

	} // </BankEvent>
	// TODO: Can we enforce monotonicity somehow?

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

		// Accrued interest
		const Real A = accruedInterest(t_i);

		stepper.add(
			t_i,
			[=] (const Interpolant1 &U, Real S) {
				return U(S) + L_hat * (A - 1.);
			},
			grid
		);
	}

	} // </InterestPayment>
	// 2015-03-04: checked; respects monotonicity; if spread is zero and
	//             events are off, then U(S_0) -> L_0 as S_0 -> infinity
	//             (i.e. the bank just gets the loan + interest back)

	// TODO: Dividends

	////////////////////////////////////////////////////////////////////////

	// Jump-diffusion operator
	BlackScholesJumpDiffusion1 bs(
		grid,

		r,      // Interest
		sigma,  // Volatility
		0.,     // Continuous dividend rate

		lambda, // Mean arrival time
		lognormal(mu_xi, sigma_xi) // Log-normal probability density
	);
	bs.setIteration(stepper);

	// No jump-diffusion test
	//BlackScholes1 bs(grid, r, sigma, 0.);

	// Discretization method
	typedef ReverseBDFTwo1 Discretization;
	Discretization discretization(grid, bs);
	discretization.setIteration(stepper);

	// Linear system solver
	BiCGSTABSolver solver;
	auto U = stepper.solve(
		grid,           // Domain
		payoff,         // Initial condition
		discretization, // Root of linear system tree
		solver          // Linear system solver
	);

	////////////////////////////////////////////////////////////////////////

	// Similarity reduction
	const Real alpha = L_hat / L_0;
	const Real value = U(alpha * S_0) / alpha;
	// 2015-02-28: checked; without events, exp(-r * T) * L_0 - value is the
	//                      price (at t=0) of a put with strike L_0

	//RectilinearGrid1 printGrid( Axis::range(0., 0.1, 200.) );
	//cout << accessor( printGrid, U );

	return value;
}

int main(int argc, char **argv) {
	cout.precision(12); // Precision
	printHeaders();

	// Populate vector of interest payment dates
	{
		interestPaymentDates.reserve(interestPayments);
		const Real dt = T / interestPayments;

		interestPaymentDates.reserve(interestPayments);
		for(int i = 0; i <= interestPayments; ++i) {
			interestPaymentDates.push_back(dt * i);
		}
	} // 2015-02-28: checked

	Real value, previousValue = nan(""), previousChange = nan("");
	for(
		int refinement = 0;
		refinement < maxRefinement;
		++refinement, N *= 2
	) { // <RefinementLoop>
		auto start = chrono::steady_clock::now();

		// TODO: Fair spread computation

		value = solve(s);

		auto end = chrono::steady_clock::now();
		auto diff = end - start;

		////////////////////////////////////////////////////////////////
		// Print table rows
		////////////////////////////////////////////////////////////////

		Real
			change = value - previousValue,
			ratio = previousChange / change
		;

		cout
			<< setw(td) << grid.size() << "\t"
			<< setw(td) << N           << "\t"
			<< setw(td) << value       << "\t"
			<< setw(td) << change      << "\t"
		;

		if(!fairSpread) {
			cout << setw(td) << ratio << "\t";
		}

		cout
			<< setw(td) << chrono::duration <double> (diff).count()
			<< endl;

		previousChange = change;
		previousValue = value;

		////////////////////////////////////////////////////////////////

		// Refine grid
		grid = grid.refined();

	} // </RefinementLoop>
}
