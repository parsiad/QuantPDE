////////////////////////////////////////////////////////////////////////////////
// nonrecourse_stock_loan.cpp
// --------------------------
//
// Author: Parsiad Azimzadeh, Peter Forsyth, Kenneth Vetzal
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>
#include <QuantPDE/Modules/Lambdas>
#include <QuantPDE/Modules/Operators>

using namespace QuantPDE;
using namespace QuantPDE::Modules;

////////////////////////////////////////////////////////////////////////////////

#include <algorithm> // max, min
#include <cmath>     // exp
#include <getopt.h>  // getopt_long
#include <iomanip>   // setw
#include <iostream>  // cout, cerr
#include <limits>    // infinity
#include <string>    // to_string

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Default constants
////////////////////////////////////////////////////////////////////////////////

/**
 * Compute the least common multiple of two integers.
 */
int lcm(int a, int b);

Real beta_lo         = 0.85;  // Low trigger
Real beta_hi         = 0.9;   // High trigger

Real S_0             = 100.;  // Initial stock value
Real L_0             = 80.;   // Initial loan value
Real L_hat;                   // Representative value

Real r               = 0.02;  // Interest rate
Real s_0             = 0.0;   // Spread
Real sigma           = 0.3;   // Volatility
Real q               = 0.;    // Discrete dividend rate

Real p               = 0.;    // Penalty for lapsation

Real lambda          = 0.05;  // Jump arrival rate
Real mu_xi           = -.8;   // Mean jump amplitude
Real sigma_xi        = .42;   // Jump amplitude standard deviation

Real T               = 5.;    // Expiry

// Number of events over course of contract
int borrowerEvents   = -1;    // Borrower events
int bankEvents       = -1;    // Bank events
int interestPayments = 20;    // Interest payments
int dividendPayments = 20;    // Dividend payments

Real firstPrepay     = 0.;    // First time at which the borrower can prepay

bool borrowerAll     = true;  // Borrower events at all times
bool bankAll         = true;  // Bank events at all times

// Initial number of steps
int N                = lcm(interestPayments, dividendPayments);

int maxRefinement    = 5;     // Maximum number of times to refine

bool dividendsToBank = false; // Dividends go to the bank

////////////////////////////////////////////////////////////////////////////////

/**
 * Controls the output of the program.
 */
enum class ProgramOperation {
	PLOT, /**< outputs plotting data */
	PLOT_SPREAD, /**< outputs plotting data with spread as x-axis */
	FAIR_SPREAD, /**< computes the fair spread */
	FIXED_SPREAD /**< convergence test for a fixed spread */
};
ProgramOperation op = ProgramOperation::FIXED_SPREAD;

////////////////////////////////////////////////////////////////////////////////

std::vector<Real> interestPaymentDates; // Sorted (ascending) vector of interest
                                        // payment dates; initialized in main()

/**
 * Similarity reduction.
 * @param U One-dimensional solution (for fixed L = L_hat)
 * @param S The stock price.
 * @param L The loan value (assumed > 0).
 * @return The value of U.
 */
Real simU(const Interpolant1 &U, Real S, Real L) {
	assert(L > 0.);
	const Real alpha = L_hat / L;
	return U(alpha * S) / alpha;
}

/**
 * Accrued interest at a particular time between coupon dates.
 * @param s The spread.
 * @param t The particular time.
 * @return e^((r+s)(t-t_p)) where t_p is the previous coupon date.
 */
Real accruedInterest(Real s, Real t) {
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
 * Payoff from the bank's perspective for representative value of L.
 * @param S The stock price.
 * @return The payoff.
 */
Real payoff(Real S) {
	return min(S, L_hat);
}

////////////////////////////////////////////////////////////////////////////////

/**
 * Prints help to stderr.
 */
void help();

/**
 * Processes command line options.
 * @param argc Number of arguments.
 * @param argv Arguments.
 * @return 0 on success.
 */
int processOptions(int argc, char **argv);

/**
 * Print selected options to stderr.
 * @param The initial grid before any refinement.
 */
void printOptions(const RectilinearGrid1 &initialGrid);

/**
 * Solve the nonrecourse stock loan problem with fixed spread.
 * @param grid The stock grid.
 * @param s A fixed value of the spread.
 * @return The solution U(S, L) as a lambda function that can be queried at any
 *         point.
 */
Function2 solve(RectilinearGrid1 &grid, Real s);

/**
 * Main code.
 */
int main(int argc, char **argv) {

	// Process options
	{
		int error;
		if( (error = processOptions(argc, argv)) != 0 ) {
			return error;
		}
	}

	////////////////////////////////////////////////////////////////////////

	// Adjust number of timesteps accordingly
	{
		int newN = interestPayments;

		if(bankEvents     > 0) { newN = lcm(newN, bankEvents      ); }
		if(borrowerEvents > 0) { newN = lcm(newN, borrowerEvents  ); }

		if(q > 0 && dividendPayments > 0) {
			newN = lcm(newN, dividendPayments);
		}

		int tmp = newN;
		while(tmp < N) { tmp += newN; }

		if(N != tmp) {
			N = tmp;
			cerr <<
"warning: changed initial number of timesteps for better convergence";
			cerr << endl << endl;
		}
	}


	// Initial grid with a node at S_0
	RectilinearGrid1 initialGrid( S_0 * Axis::special );

	// Representative value of L_hat is S_0 (ideal for choice of grid)
	L_hat = S_0;

	////////////////////////////////////////////////////////////////////////

	printOptions(initialGrid);

	////////////////////////////////////////////////////////////////////////

	// Populate vector of interest payment dates
	{
		//interestPaymentDates.reserve(interestPayments);
		const Real dt = T / interestPayments;

		// i = 0 is not an interest payment date!
		// It is added in the list but not processed (do not remove)

		for(int i = 0; i <= interestPayments; ++i) {
			interestPaymentDates.push_back(dt * i);
		}
	} // 2015-02-28: checked

	cout.precision(12); // Precision
	constexpr int td = 20; // Spacing used to print

	if(op != ProgramOperation::FIXED_SPREAD) {

		// Double the number of timesteps as many times as necessary
		for(int i = 0; i < maxRefinement; ++i) { N *= 2; }

		// Refine grid ref times
		auto grid = initialGrid.refined( maxRefinement );

		if(op == ProgramOperation::PLOT) {

			// Plot

			// Fixed spread computation
			auto U = solve(grid, s_0);

			// Print U(S, L_0) (fixed L_0) at all grid nodes
			cout << accessor(
				grid,
				[&] (Real S) { return U(S, L_0); }
			);

		} else if(op == ProgramOperation::PLOT_SPREAD) {

			// Plot spread

			Real hi = r; // Initial guess for upper bound is r

			// Find upper bound
			while(1) {
				auto U = solve(grid, hi);
				if( U(S_0, L_0) >= L_0 ) {
					break;
				}
				hi *= 2; // Double upper bound
			}

			const Real ds = 0.001; // Make this an option

			RectilinearGrid1 spreads(Axis::range(0, ds, hi));

			// U as a function of the spread
			auto U_spread = [&] (Real spread) {
				auto U = solve(grid, spread);
				return U(S_0, L_0);
			};

			// Print U_spread for all spreads on the grid
			cout << accessor(spreads, U_spread);

		} else {
			// Fair spread

			// Print headers
			cout
				<< setw(td) << "Spread"       << "\t"
				<< setw(td) << "Value at S_0" << "\t"
				<< endl
			;

			// Initial brackets
			Real lo = 0.;
			Real hi = r; // Initial guess for upper bound

			// Find upper bound
			while(1) {

				auto U = solve(grid, hi);
				const Real value = U(S_0, L_0);

				cout
					<< setw(td) << hi    << "\t"
					<< setw(td) << value << "\t"
				;

				if( value >= L_0 ) {
					cerr << "\t" << "bracketed in ["
							<< lo << ", " << hi
							<< "]; starting hybrid "
							"Newton" << endl;
					break;
				}

				cout << endl;

				hi *= 2; // Double upper bound
			}

			// Hybrid Newton's method
			Real s_0 = (hi + lo) / 2., s_new;
			while(1) {
				auto U_0 = solve(grid, s_0);
				const Real f_0 = U_0(S_0, L_0);

				// Update brackets
				if(f_0 >= L_0 && s_0 < hi) {
					hi = s_0;
				} else if(f_0 < L_0 && s_0 > lo) {
					lo = s_0;
				}

				auto U_1 = solve(grid, s_0
						+ QuantPDE::epsilon);
				const Real f_1 = U_1(S_0, L_0);

				cout
					<< setw(td) << s_0 << "\t"
					<< setw(td) << f_0 << "\t"
					<< endl
				;

				if(f_1 - f_0 < QuantPDE::epsilon) {
					// Binary search
					s_new = (lo + hi) / 2;
				} else {
					// Newton iteration
					s_new = s_0 - QuantPDE::epsilon
							* (f_0 - L_0)
							/ (f_1 - f_0);
				}

				const Real err = fabs(s_new - s_0);
				if(err < QuantPDE::tolerance) {
					// Found fair fee; break
					break;
				}

				s_0 = s_new;
			}
		}

	} else {

		// Print headers
		cout
			<< setw(td) << "Nodes"           << "\t"
			<< setw(td) << "Time Steps"      << "\t"
			<< setw(td) << "Value"           << "\t"
			<< setw(td) << "Change"          << "\t"
			<< setw(td) << "Ratio"
			<< endl
		;

		Real value, previousValue = nan(""), previousChange = nan("");
		for(int ref = 0; ref <= maxRefinement; ++ref, N *= 2) {
			// Refine grid ref times
			auto grid = initialGrid.refined( ref );

			// Fixed spread computation
			auto U = solve(grid, s_0);
			value = U(S_0, L_0);

			////////////////////////////////////////////////////////
			// Print table rows
			////////////////////////////////////////////////////////

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
				<< endl;

			previousChange = change;
			previousValue = value;
		}

	}

	return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Solution for fixed spread
////////////////////////////////////////////////////////////////////////////////

Function2 solve(RectilinearGrid1 &grid, Real s) {
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
	if(borrowerAll) { borrowerEvents = N; }
	if(bankAll)     { bankEvents     = N; }

	////////////////////////////////////////////////////////////////////////

	{ // <BorrowerEvent>

	// Time between borrower events
	const Real dt = T / borrowerEvents;

	for(int i = 0; i < borrowerEvents; ++i) {
		const Real t_i = dt * i;

		// Accrued interest
		const Real A = accruedInterest(s, t_i);

		stepper.add(
			t_i,
			[=] (const Interpolant1 &U, Real S) {
				Real newValue = S;
				if(t_i >= firstPrepay) {
					const Real prepay = (1+p) * L_hat * A;
					if(prepay < newValue) {
						newValue = prepay;
					}
				}

				return min( U(S), newValue );
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
		const Real A = accruedInterest(s, t_i);

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
				}

				// R is above low trigger
				if(R > beta_lo) {
					const Real R_0 = L_0 / S_0;

					// Top-up with shares
					/*
					{
						const Real tmp = U(L_hat / R_0);
						if(tmp > best) {
							best = tmp;
						}
					}
					*/

					// Top-up with cash
					// (similarity reduction)
					{
						const Real tmp =
							simU(U, S, S * R_0)
							+ (L_hat - S * R_0) * A
						;
						if(tmp > best) {
							best = tmp;
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
			const Real A = accruedInterest(s, t_i);

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

	////////////////////////////////////////////////////////////////////////

	if(q > 0.) { // <DividendPayments>

		// Add dividend payment events (skip initial time)
		const Real dt = T / dividendPayments;

		for(int i = 1; i <= dividendPayments; ++i) {
			// Time at which interest payment takes place
			const Real t_i = dt * i;

			// Accrued interest
			const Real A = accruedInterest(s, t_i);

			if(dividendsToBank) {
				stepper.add(
					t_i,
					[=] (const Interpolant1 &U, Real S) {
						// Bank gets dividends
						Real Lp = L_hat - q * S / A;
						if(Lp < QuantPDE::epsilon) {
							Lp = QuantPDE::epsilon;
							// Saves us from solving
							// another 1d PDE
						}

						// Similarity reduction
						return simU(U, (1-q) * S, Lp)
								+ q * S;
					},
					grid
				);
			} else {
				stepper.add(
					t_i,
					[=] (const Interpolant1 &U, Real S) {
						// Borrower gets dividends
						return U( (1-q) * S );
					},
					grid
				);
			}
		}

	} // </DividendPayments>

	////////////////////////////////////////////////////////////////////////

	// Jump-diffusion operator
	BlackScholesJumpDiffusion1 bs(
		grid,

		r,      // Interest
		sigma,  // Volatility
		0.,     // No continuous dividends

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
	return [=] (Real S, Real L) {
		if(L < QuantPDE::epsilon) { L = QuantPDE::epsilon; }
		const Real alpha = L_hat / L;
		const Real value = U(alpha * S) / alpha;
		return value;
	};
	// 2015-02-28: checked; without events, exp(-r * T) * L_0 - value is the
	//                      price (at t=0) of a put with strike L_0
}

////////////////////////////////////////////////////////////////////////////////
// Least common multiple
////////////////////////////////////////////////////////////////////////////////

int lcm(int a, int b)
{
	int gcd;
	{
		int a_ = a, b_ = b;
		for (;;)
		{
			if (a_ == 0) {
				gcd = b_;
				break;
			}

			b_ %= a_;

			if (b_ == 0) {
				gcd = a_;
				break;
			}

			a_ %= b_;
		}
	}

	return gcd ? (a / gcd * b) : 0;
}

////////////////////////////////////////////////////////////////////////////////
// Help and options
////////////////////////////////////////////////////////////////////////////////

void help() {
	cerr <<
"nonrecourse_stock_loan [OPTIONS]" << endl << endl <<
"Used to analyze a nonrecourse stock loan." << endl <<
endl <<
"-g REAL" << endl <<
endl <<
"    Sets the jump amplitude standard deviation (default is 0.42)." << endl <<
endl <<
"-l NONNEGATIVE_REAL" << endl <<
endl <<
"    Sets the mean arrival time for the jump process (default is 0.1)." << endl <<
endl <<
"-L NONNEGATIVE_REAL" << endl <<
endl <<
"    Sets the initial loan value (default is 80.)." << endl <<
endl <<
"-m REAL" << endl <<
endl <<
"    Sets the mean jump amplitude (default is -0.8)." << endl <<
endl <<
"-N POSITIVE_INTEGER" << endl <<
endl <<
"    Sets the initial number of timesteps (default is 12)." << endl <<
endl <<
"-p PROPORTION" << endl <<
endl <<
"    Sets the penalty rate for borrower lapsation (default is 0.)." << endl <<
endl <<
"-q PROPORTION" << endl <<
endl <<
"    Sets the dividend rate (default is 0.)." << endl <<
endl <<
"-r REAL" << endl <<
endl <<
"    Sets the interest rate (default is 0.02)." << endl <<
endl <<
"-R NONNEGATIVE_INTEGER" << endl <<
endl <<
"    Controls the coarseness of the grid, with 0 being coarsest (default is 5)." << endl <<
endl <<
"-s REAL" << endl <<
endl <<
"    Sets the spread (default is 0.)." << endl <<
endl <<
"-S NONNEGATIVE_REAL" << endl <<
endl <<
"    Sets the initial stock price (default is 100.)." << endl <<
endl <<
"-T POSITIVE_REAL" << endl <<
endl <<
"    Sets the expiry time (default is 5.)." << endl <<
endl <<
"-v REAL" << endl <<
endl <<
"    Sets the volatility (default is 0.3)." << endl <<
endl <<
"--bank-lo-trigger POSITIVE_REAL" << endl <<
endl <<
"    Sets the bank top-up trigger (default is 0.85)." << endl <<
endl <<
"--bank-hi-trigger POSITIVE_REAL" << endl <<
endl <<
"    Sets the bank liquidation trigger (default is 0.9)." << endl <<
endl <<
"--bank-events INTEGER" << endl <<
endl <<
"    Number of bank events; -1 for at all time steps (default is -1)." << endl <<
endl <<
"--borrower-events INTEGER" << endl <<
endl <<
"    Number of borrower events; -1 for at all time steps (default is -1)." << endl <<
endl <<
"--coupons INTEGER" << endl <<
endl <<
"    Number of coupon payments (default is 20)." << endl <<
endl <<
"--dividend-payments INTEGER" << endl <<
endl <<
"    Number of dividend payments (default is 20)." << endl <<
endl <<
"--bank-gets-dividends" << endl <<
endl <<
"    Specifies that the bank gets the dividends (default is off)." << endl <<
endl <<
"--no-prepay-until" << endl <<
endl <<
"    Specifies the first time at which the borrower can prepay (default is 0.)." << endl <<
endl <<
"--plot" << endl <<
endl <<
"    Plot data for the solution vs. the initial stock price (default is off)." << endl <<
endl <<
"--plot-spread" << endl <<
endl <<
"    Plot data for the solution vs. the spread (default is off)." << endl <<
endl <<
"--fair-spread" << endl <<
endl <<
"    Computes the fair spread (default is off)." << endl <<
endl;
}

int processOptions(int argc, char **argv) {

	{
		// Long option names
		static struct option opts[] = {
			{ "bank-events",         required_argument, 0, 0 },
			{ "borrower-events",     required_argument, 0, 0 },
			{ "coupons",             required_argument, 0, 0 },
			{ "dividend-payments",   required_argument, 0, 0 },
			{ "bank-gets-dividends", no_argument,       0, 0 },
			{ "bank-lo-trigger",     required_argument, 0, 0 },
			{ "bank-hi-trigger",     required_argument, 0, 0 },
			{ "plot",                no_argument,       0, 0 },
			{ "plot-spread",         no_argument,       0, 0 },
			{ "fair-spread",         no_argument,       0, 0 },
			{ "no-prepay-until",     required_argument, 0, 0 },
			{ nullptr, 0, 0, 0 }
		};

		int c;
		int index;
		while(
			(
				c = getopt_long(
					argc,
					argv,
					"g:hl:L:m:N:p:q:r:R:s:S:T:v:",
					opts,
					&index
				)
			) != -1
		) {
			switch(c) {
				// Long options
				case 0:
					switch(index)
					{
						case 0:
						bankEvents = atoi(optarg);
						if(bankEvents < 0) {
							bankAll = true;
						} else {
							bankAll = false;
						}
						break;

						case 1:
						borrowerEvents = atoi(optarg);
						if(borrowerEvents < 0) {
							borrowerAll = true;
						} else {
							borrowerAll = false;
						}
						break;

						case 2:
						interestPayments = atoi(optarg);
						if(interestPayments < 1) {
							cerr <<
"error: there must be at least one coupon payment" << endl;
							return 1;
						}
						break;

						case 3:
						dividendPayments = atoi(optarg);
						if(dividendPayments < 0) {
							cerr <<
"error: number of dividend payments must be nonnegative" << endl;
							return 1;
						}
						break;

						case 4:
						dividendsToBank = true;
						break;

						case 5:
						beta_lo = atof(optarg);
						break;

						case 6:
						beta_hi = atof(optarg);
						break;

						case 7:
						op = ProgramOperation
								::PLOT;
						break;

						case 8:
						op = ProgramOperation
								::PLOT_SPREAD;
						break;

						case 9:
						op = ProgramOperation
								::FAIR_SPREAD;
						break;

						case 10:
						firstPrepay = atof(optarg);
						break;

						default:
							break;
					}
				break;

				case 'g':
				sigma_xi = atof(optarg);
				break;

				case 'h':
				help();
				return 1;

				case 'l':
				if((lambda = atof(optarg)) < 0) {
					cerr <<
"error: the mean arrival time must be nonnegative" << endl;
					return 1;
				}
				break;

				case 'L':
				if((L_0 = atof(optarg)) < 0.) {
					cerr <<
"error: the initial loan value must be nonnegative" << endl;
					return 1;
				}
				break;

				case 'm':
				mu_xi = atof(optarg);
				break;

				case 'N':
				if((N = atoi(optarg)) <= 0) {
					cerr <<
"error: the number of steps must be positive" << endl;
					return 1;
				}
				break;

				case 'p':
				p = atof(optarg);
				if(p < 0. || p > 1.) {
					cerr <<
"error: the penalty rate must be a proportion" << endl;
					return 1;
				}
				break;

				case 'q':
				q = atof(optarg);
				if(q < 0. || q > 1.) {
					cerr <<
"error: the dividend rate must be a proportion" << endl;
					return 1;
				}
				break;

				case 'r':
				r = atof(optarg);
				break;

				case 'R':
				if((maxRefinement = atoi(optarg)) < 0) {
					cerr <<
"error: the maximum level of refinement must be nonnegative" << endl;
					return 1;
				}
				break;

				case 's':
				s_0 = atof(optarg);
				break;

				case 'S':
				if((S_0 = atof(optarg)) < 0.) {
					cerr <<
"error: the initial stock price must be nonnegative" << endl;
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
				sigma = atof(optarg);
				break;

				case ':':
				case '?':
				cerr << "specify -h for a list of options"
						<< endl;
				return 1;
			}
		}
	}

	return 0;
}

void printOptions(const RectilinearGrid1 &initialGrid) {
	cerr
		<< "Initial stock price:     \t" << S_0
		<< endl
		<< "Initial loan value:      \t" << L_0
		<< endl
		<< endl

		<< "Interest rate:           \t" << r
		<< endl
		<< "Spread:                  \t" << s_0
		<< endl
		<< "Volatility:              \t" << sigma
		<< endl
		<< endl

		<< "Mean jump arrival time:  \t" << lambda
		<< endl
		<< "Mean jump amplitude:     \t" << mu_xi
		<< endl
		<< "Jump amplitude std:      \t" << sigma_xi
		<< endl
		<< endl

		<< "Expiry time:             \t" << T
		<< endl
		<< endl

		<< "Penalty rate:            \t" << p
		<< endl
		<< "Dividend rate:           \t" << q
		<< endl
		<< endl

		<< "Bank top-up trigger:     \t" << beta_lo
		<< endl
		<< "Bank liquidation trigger:\t" << beta_hi
		<< endl
		<< endl

		<< "No prepay until:         \t" << firstPrepay
		<< endl
		<< endl

		<< "Bank events:             \t" << ( bankEvents < 0
				? "at all times" : to_string(bankEvents) )
		<< endl
		<< "Borrower events:         \t" << ( borrowerEvents < 0
				? "at all times" : to_string(borrowerEvents) )
		<< endl
		<< "Coupon payments:         \t" << interestPayments
		<< endl
		<< "Dividend payments:       \t" << dividendPayments
		<< endl
		<< "Bank gets dividends:     \t" <<
				(dividendsToBank ? "true" : "false")
		<< endl
		<< endl

		<< "Initial # of time steps: \t" << N
		<< endl
		<< "Initial grid:            \t" << initialGrid
		<< endl
		<< "Maximum refinement:      \t" << maxRefinement
		<< endl
		<< endl
	;
}

