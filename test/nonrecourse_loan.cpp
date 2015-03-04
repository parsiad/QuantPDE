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

#include <algorithm> // min
#include <iomanip>   // setw
#include <iostream>  // cout, cerr
#include <cmath>     // exp

using namespace std;

////////////////////////////////////////////////////////////////////////////////

Real beta_lo         = 0.8;   // Low trigger
Real beta_hi         = 0.9;   // High trigger

Real S_0             = 125.;  // Initial stock value
Real L_0             = 100.;  // Initial loan value
Real L_hat           = 100.;  // Representative value of L

Real r               = 0.04;  // Interest rate
Real s               = 0.02;  // Spread
Real sigma           = 0.2;   // Volatility

Real lambda          = 0.05;  // Jump arrival rate
Real mu_xi           = -.8;   // Mean jump amplitude
Real sigma_xi        = .42;   // Jump amplitude standard deviation

Real T               = 1.;    // Expiry

int events           = -1;    // Bank and borrower events (-1 for all times)
int interestPayments = 12;    // Interest payments

int N_0              = 12;    // Initial number of steps
int maxRefinement    = 8;     // Maximum number of times to refine

bool newton          = false; // TODO: Newton iteration to find optimal spread

////////////////////////////////////////////////////////////////////////////////

std::vector<Real> interestPaymentDates; // Sorted (ascending) vector of interest
                                        // payment dates

// Payoff for fixed \hat{L}
Real payoff(Real S) {
	return min(S, L_hat);
}

// Accrued interest A(t)
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

////////////////////////////////////////////////////////////////////////////////

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
		interestPaymentDates.reserve(interestPayments);
		const Real dt = T / interestPayments;

		interestPaymentDates.reserve(interestPayments);
		for(int i = 0; i <= interestPayments; ++i) {
			interestPaymentDates.push_back(dt * i);
		}
	} // 2015-02-28: checked

	// Hand-picked grid
	// TODO: Make this more robust; place nodes at S_0 and L_0
	RectilinearGrid1 grid(
		(S_0 / 100.) * Axis {
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
	// 2. Bank: top-up or liquidate
	// 3. Interest payment and (TODO) dividends

	// Apply events at all times
	if(events < 0) { events = N; }

	////////////////////////////////////////////////////////////////////////

	{ // <BorrowerEvent>

	// Time between borrower events
	const Real dt = T / events;

	for(int i = 0; i < events; ++i) {
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

	////////////////////////////////////////////////////////////////////////

	{ // <BankEvent>

	// Time between bank events
	const Real dt = T / events;

	for(int i = 0; i < events; ++i) {
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
					const Real newPos = min(L_hat * A, S);
					if(newPos > best) {
						best = newPos;
					}

				// R is between low and high triggers
				} else if(R >= beta_lo) {

					const Real R_0 = L_0 / S_0;

					{ // Top-up with shares
						const Real newPos =
								U(L_hat / R_0);
						if(newPos > best) {
							best = newPos;
						}
					}

					// Top-up with cash
					if(S < epsilon) {
						const Real newPos = L_hat * A;
						if(newPos > best) {
							best = newPos;
						}
					} else {
						const Real newPos =
							L_hat / (S * R_0)
							* U(L_hat / R_0)
							+ (L_hat - S * R_0) * A
						;
						if(newPos > best) {
							best = newPos;
						}
					}
				}

				return best;
			},
			grid
		);
	}

	} // </BankEvent>

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
	Real value = U(alpha * S_0) / alpha;
	// 2015-02-28: checked; without events, exp(-r * T) * L_0 - value is the
	//                      price (at t=0) of a put with strike L_0

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
