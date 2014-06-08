#include <QuantPDE/Core>

// TODO: Change these includes; shouldn't include src directory explicitly
#include <QuantPDE/src/Modules/Payoffs.hpp>
#include <QuantPDE/src/Modules/BlackScholesOperator.hpp>

using namespace QuantPDE;
using namespace QuantPDE::Modules;

///////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <unistd.h>

using namespace std;

///////////////////////////////////////////////////////////////////////////////

// TODO: Remove

template <size_t Lookback = 1>
class AmericanPutPenalty : public Linearizer<Lookback> {

	const Domain1 *domain;
	LinearizerBase *left;

	Real strike, tolerance, large;
	Matrix P;

public:

	template <typename D>
	AmericanPutPenalty(D &domain, LinearizerBase &left, Real strike,
			Real tolerance = 1e-6) noexcept : domain(&domain),
			left(&left), strike(strike), tolerance(tolerance),
			P( domain.size(), domain.size() ) {
		P.reserve(IntegerVector::Constant(domain.size(), 1));
	}

	virtual void onIterationStart() {
		const Vector &v = std::get<1>(this->iterands()[0]);

		P.setZero();

		for(Index i = 0; i < domain->size(); i++) {
			auto x = domain->coordinates(i);
			if( v(i) < strike - x[0] - tolerance ) {
				P.insert(i, i) = 1. / tolerance;
			}
		}
	}

	virtual Matrix A() {
		return left->A() + P;
	}

	virtual Vector b() {
		return left->b() + P * domain->image(
				QUANT_PDE_MODULES_PUT_PAYOFF_FIXED_STRIKE(
				strike));
	}

};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {

	// Default options
	Real expiry         = 1.;
	Real interest       = 0.04;
	Real volatility     = 0.2;
	Real dividends      = 0.;
	Real stock          = 100.;
	Real strike         = 100.;
	unsigned refinement = 5;
	unsigned steps      = 25;
	bool call           = true;

	// Setting options with getopt
	{ char c;
	while((c = getopt(argc, argv, "d:hK:pr:R:s:S:T:v:")) != -1) {
		switch(c) {
			case 'd':
				dividends = atof(optarg);
				break;
			case 'h':
				cerr <<
"ConvergenceTableEuropean [OPTIONS]" << endl << endl <<
"Outputs the rate of convergence for computing the price of a European call or" << endl <<
"put using an upwind, fully implicit discretization of the Black-Scholes partial" << endl <<
"differential equation." << endl <<
endl <<
"-d REAL" << endl <<
endl <<
"    sets the dividend rate (default is 0.)" << endl <<
"-K REAL" << endl <<
endl <<
"    sets the strike price (default is 100.)" << endl <<
endl <<
"-p" << endl <<
endl <<
"    computes the price of a European put (default is call)" << endl <<
endl <<
"-r REAL" << endl <<
endl <<
"    sets interest rate (default is 0.04)" << endl <<
endl <<
"-R POSITIVE_INTEGER" << endl <<
endl <<
"    sets the maximum number of refinement steps in the computation (default is" << endl <<
"    5)" << endl <<
endl <<
"-s POSITIVE_INTEGER" << endl <<
endl <<
"    sets the initial number of steps to take in time (default is 25)" << endl <<
endl <<
"-S REAL" << endl <<
endl <<
"    sets the initial stock price (default is 100.)" << endl <<
endl <<
"-T POSITIVE_REAL" << endl <<
endl <<
"    sets the expiry time (default is 1.)" << endl <<
endl <<
"-v REAL" << endl <<
endl <<
"    sets the volatility" << endl << endl;
				return 0;
			case 'K':
				strike = atof(optarg);
				break;
			case 'p':
				call = false;
				break;
			case 'r':
				interest = atof(optarg);
				break;
			case 'R':
				refinement = atoi(optarg);
				break;
			case 's':
				steps = atoi(optarg);
				if(steps <= 0) {
					cerr <<
"error: the number of steps must be positive" << endl;
					return 1;
				}
				break;
			case 'S':
				stock = atof(optarg);
				break;
			case 'T':
				if((expiry = atof(optarg)) <= 0.) {
					cerr <<
"error: expiry time must be positive" << endl;
					return 1;
				}
				break;
			case 'v':
				volatility = atof(optarg);
				break;
		}
	} }

	// Setting up the table
	const int td = 10;
	cout
		<< setw(td) << "Nodes"  << "\t"
		<< setw(td) << "Steps"  << "\t"
		<< setw(td) << "Value"  << "\t"
		<< setw(td) << "Change" << "\t"
		<< setw(td) << "Ratio"
		<< endl
	;
	cout.precision(6);
	Real previousValue = nan(""), previousChange = nan("");

	// Initial discretization
	// TODO: Create grid based on initial stock price and strike price
	Axis S {
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
	RectilinearGrid1 R(S);

	// Payoff function
	Function1 payoff = bind(
		call ? callPayoff : putPayoff,
		placeholders::_1,
		strike
	);

	// Alternatively, we could have used...
	//auto payoff = QUANT_PDE_MODULES_CALL_PAYOFF_FIXED_STRIKE(strike);
	//auto payoff = QUANT_PDE_MODULES_PUT_PAYOFF_FIXED_STRIKE(strike);

	for(unsigned l = 0; l < refinement; l++, steps *= 2) {

		///////////////////////////////////////////////////////////////
		// Build spatial grid
		///////////////////////////////////////////////////////////////

		// Refine the grid
		R.refine( RectilinearGrid1::NewTickBetweenEachPair() );

		///////////////////////////////////////////////////////////////
		// Build problem
		///////////////////////////////////////////////////////////////

		/*
		EuropeanOption europeanOption(
			payoff,
			[interest]   (Real, Real) { return interest;   },
			[volatility] (Real, Real) { return volatility; },
			[dividends]  (Real, Real) { return dividends;  },
			0., expiry
		);
		*/

		///////////////////////////////////////////////////////////////
		// Solve problem
		///////////////////////////////////////////////////////////////

		BlackScholesOperator blackScholes(
			R,
			[interest]   (Real, Real) { return interest;   },
			[volatility] (Real, Real) { return volatility; },
			[dividends]  (Real, Real) { return dividends;  }
		);

		ConstantStepper<6> stepper(
			0., // Initial time
			expiry,
			steps
		);
		//ToleranceIteration<> tolerance;
		//stepper.setChildIteration(tolerance);

		LinearBDFSix<> bdf(R, blackScholes);
		bdf.setIteration(stepper);

		//AmericanPutPenalty<> penalty(R, bdf, strike);
		//penalty.setIteration(tolerance);

		BiCGSTABSolver solver;
		Vector v = stepper.iterateUntilDone(
			R.image(payoff),
			bdf, //penalty,
			solver
		);

		///////////////////////////////////////////////////////////////
		// Table
		///////////////////////////////////////////////////////////////

		// Solution at S = 100.
		Real value = R.accessor( v )(stock);

		// Change and ratio between successive solutions
		Real
			change = value - previousValue,
			ratio = previousChange / change
		;

		// Print out row of table
		cout
			<< scientific
			<< setw(td) << R.size() << "\t"
			<< setw(td) << steps    << "\t"
			<< setw(td) << value    << "\t"
			<< setw(td) << change   << "\t"
			<< setw(td) << ratio
			<< endl
		;

		previousChange = change;
		previousValue = value;

	}

	return 0;
}

