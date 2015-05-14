////////////////////////////////////////////////////////////////////////////////
// fex_rate.cpp
// ------------
//
// Optimal combined control of the exchange rate, as formulated in [1].
//
// [1] Mundaca, Gabriela, and Bernt Øksendal. "Optimal stochastic intervention
// control with application to the exchange rate." Journal of Mathematical
// Economics 29.2 (1998): 225-243.
//
////////////////////////////////////////////////////////////////////////////////

//#define QUANT_PDE_FEX_RATE_EVEN

#include <QuantPDE/Core>

#include <cmath>    // fabs
#include <fstream>  // ofstream
#include <getopt.h> // getopt_long
#include <iomanip>  // setw
#include <iostream> // cout, cerr
#include <limits>   // numeric_limits
#include <memory>   // unique_ptr
#include <numeric>  // accumulate
#include <string>   // to_string

////////////////////////////////////////////////////////////////////////////////

using namespace QuantPDE;
using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Methods
////////////////////////////////////////////////////////////////////////////////

constexpr int SEMI_LAGRANGIAN_WITHDRAWAL_CONTINUOUS = 1 << 0;
constexpr int EXPLICIT_IMPULSE = 1 << 1;

constexpr int EXPLICIT =
		  EXPLICIT_IMPULSE
		| SEMI_LAGRANGIAN_WITHDRAWAL_CONTINUOUS;

constexpr int IMPLICIT = 0;

////////////////////////////////////////////////////////////////////////////////

int method = 0;

// The problem is solved on [-boundary, boundary]
Real boundary = 6.; // log(1000) ~= 7

// Interest rate differential in [q_min, q_max]
Real q_min = 0.;
Real q_max = 0.05;

// Number of points in space (initial discretization)
int gridPoints = 32;

// Number of points for interest rate control (initial discretization)
int controlPoints = 16;

// Constants
Real c_1 = 1.;
Real c_2 = 1.; // lambda
Real c_3 = QuantPDE::epsilon; // c
Real c_4 = .25; // a
Real c_5 = 1.;
Real c_6 = 3.; // b

// Discount
Real rho = 0.02;

// Volatility
Real v = 0.3;

// Central parity
Real m = 0.;

// Initial (log of the) exchange rate
Real x_0 = 0.;

// Finite horizon
Real T = -1;

// Initial number of timesteps
int TN = 32;

// Max/min refinement level
int Rmin = 0;
int Rmax = 6;

////////////////////////////////////////////////////////////////////////////////

/**
 * Controls the output of the program.
 */
enum class ProgramOperation {
	CONVERGENCE_TEST, /**< convergence test */
	PLOT, /**< prints plotting data for the solution at time 0 */
	STOCHASTIC, /**< prints stochastic control data at time 0 */
	IMPULSE /** < prints impulse control data at time 0 */
};
ProgramOperation op = ProgramOperation::CONVERGENCE_TEST;

////////////////////////////////////////////////////////////////////////////////

// Spacing
constexpr int td = 20;

////////////////////////////////////////////////////////////////////////////////
// Functions
////////////////////////////////////////////////////////////////////////////////

// gamma(zeta) = c_1 * zeta
// L(zeta)     = c_2 * |zeta| + c_3
// F(y)        = c_4 * y
// M(y)        = c_5 * y^2
// N(y)        = c_6 * y^2

/**
* Cost of an impulse.
*/
inline Real interventionCost(Real t, Real x, Real x_j) {
	const Real zeta = (x_j - x) / c_1;
	return -(c_2 * fabs(zeta) + c_3);
}

/**
* State of the exchange rate after an impulse.
*/
inline Real xplus(Real t, Real x, Real x_j) {
	return x_j;
}

/**
 * Effect of the interest rate differential on the exchange rate.
 *
 */
inline Real F(Real y) {
	return c_4 * y;
}

/**
* Cost of exchange rate differential.
*/
inline Real M(Real y) {
	#ifdef QUANT_PDE_FEX_RATE_EVEN
	return c_5 * y * y;
	#else
	return (y > 0.) ? (c_5 * y * y) : 0.;
	#endif
}

/**
* Cost of interest rate differential.
*/
inline Real N(Real y) {
	return c_6 * y * y;
}

////////////////////////////////////////////////////////////////////////////////

class Discretizee final : public RawControlledLinearSystem1_1 {

	const RectilinearGrid1 &grid;

	const Function1 F, M, N;
	const Real rho, v, m;

	const bool controlled;
	const Vector zero;

	// TODO: Nonconstant coefficients

public:

	template <typename G, typename F1, typename F2, typename F3>
	Discretizee(
		G &grid,
		F1 &&interestDifferentialEffect,
		F2 &&exchangeRateCost,
		F3 &&interestDifferentialCost,
		Real discount,
		Real volatility,
		Real centralParity = 0.,
		bool controlled = true
	) noexcept :
		grid(grid),
		F( std::forward<F1>(interestDifferentialEffect) ),
		M( std::forward<F3>(exchangeRateCost) ),
		N( std::forward<F2>(interestDifferentialCost) ),
		rho( discount ),
		v( volatility ),
		m( centralParity ),
		controlled( controlled ),
		zero( grid.zero() )
	{
	}

	virtual Matrix A(Real time) {
		Matrix A = grid.matrix();
		A.reserve( IntegerVector::Constant(grid.size(), 3) );

		// x axis
		const Axis &x = grid[0];
		const int n = x.size();

		// Control as a vector
		const Vector &raw = controlled ? control(0) : zero;

		// Left boundary (derivatives disappear)
		A.insert(0, 0) = rho;

		// Interior point
		for(Index i = 1; i < n - 1; ++i) {
			const Real
				dxb = x[i]     - x[i - 1],
				dxc = x[i + 1] - x[i - 1],
				dxf = x[i + 1] - x[i]
			;

			const Real tmp1 = v * v;
			const Real tmp2 = -F( raw(i) );

			const Real alpha_common = tmp1 / dxb / dxc;
			const Real  beta_common = tmp1 / dxf / dxc;

			// Central
			Real alpha_i = alpha_common - tmp2 / dxc;
			Real beta_i  =  beta_common + tmp2 / dxc;
			if(alpha_i < 0.) {
				// Forward
				alpha_i = alpha_common;
				beta_i  =  beta_common + tmp2 / dxf;
			} else if(beta_i < 0.) {
				// Backward
				alpha_i = alpha_common - tmp2 / dxb;
				beta_i  =  beta_common;
			}

			assert(alpha_i >= 0.);
			assert( beta_i >= 0.);

			A.insert(i, i - 1) = -alpha_i;
			A.insert(i, i    ) =  alpha_i + beta_i + rho;
			A.insert(i, i + 1) = -beta_i;
		}

		#ifdef QUANT_PDE_FEX_RATE_EVEN
		{
			const Real dx = x[n - 1] - x[n - 2];

			const Real tmp1 = v * v;
			const Real tmp2 = -F( raw(n - 1) );

			const Real alpha_common = tmp1 / dx / (2 * dx);
			const Real  beta_common = tmp1 / dx / (2 * dx);

			// Central
			Real alpha_i = alpha_common - tmp2 / (2 * dx);
			Real beta_i  =  beta_common + tmp2 / (2 * dx);
			if(alpha_i < 0.) {
				// Forward
				alpha_i = alpha_common;
				beta_i  =  beta_common + tmp2 / dx;
			} else if(beta_i < 0.) {
				// Backward
				alpha_i = alpha_common - tmp2 / dx;
				beta_i  =  beta_common;
			}

			assert(alpha_i >= 0.);
			assert( beta_i >= 0.);

			A.insert(n - 1, n - 2) = -(alpha_i + beta_i);
			A.insert(n - 1, n - 1) =   alpha_i + beta_i + rho;
		}
		#else
		A.insert(n - 1, n - 1) = rho;
		#endif

		A.makeCompressed();
		return A;
	}

	virtual Vector b(Real time) {
		Vector b = grid.vector();

		const Axis &x = grid[0];

		// Control as a vector
		const Vector &raw = controlled ? control(0) : zero;

		// Interior and boundary points
		for(Index i = 0; i < x.size(); ++i) {
			b(i) = -( M( x[i] - m ) + N( raw(i) ) );
		}

		return b;
	}

	virtual bool isATheSame() const {
		// TODO: Check the coefficients
		return true;
	}

};

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {

	////////////////////////////////////////////////////////////////////////
	// Options
	////////////////////////////////////////////////////////////////////////

	{
		// Long option names
		static struct option opts[] = {
			{ "mixed"           , no_argument,       0, 0 },
			{ "explicit"        , no_argument,       0, 0 },
			{ "plot"            , no_argument,       0, 0 },
			{ "stochastic"      , no_argument,       0, 0 },
			{ "impulse"         , no_argument,       0, 0 },
			{ "min-refinement"  , required_argument, 0, 0 },
			{ "max-refinement"  , required_argument, 0, 0 },
			{ "expiry"          , required_argument, 0, 0 },
			{ "fixed-cost"      , required_argument, 0, 0 },
			{ nullptr           , 0,                 0, 0 }
		};

		int c;
		int index;
		while(
			(
				c = getopt_long(
					argc,
					argv,
					"",
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
							method =
					SEMI_LAGRANGIAN_WITHDRAWAL_CONTINUOUS;
							break;
						case 1:
							method = EXPLICIT;
							break;
						case 2:
							op = ProgramOperation::
									PLOT;
							break;
						case 3:
							op = ProgramOperation::
								STOCHASTIC;
							break;
						case 4:
							op = ProgramOperation::
									IMPULSE;
							break;
						case 5:
							Rmin = atoi(optarg);
							break;
						case 6:
							Rmax = atoi(optarg);
							break;
						case 7:
							T = atof(optarg);
							if(T <= 0.) {
								cerr <<
"error: expiry time must be positive" << endl;
								return 1;
							}
							break;
						case 8:
							c_3 = atof(optarg);
							if(c_3 < 0.) {
								cerr <<
"error: fixed transaction cost must be nonnegative" << endl;
								return 1;
							}
							break;
						default:
							break;
					}
				break;
			}
		}

		if(Rmin < 0 || Rmax < Rmin) {
			cerr << "error: the minimum level of refinement must be"
					" nonnegative and less than or equal to"
					" the maximum level of refinement"
					<< endl;
			return 1;
		}

		if(method != IMPLICIT && T <= 0.) {
			cerr << "error: an implicit method must be used for the"
					" infinite horizon solution" << endl;
			return 1;
		}

		if( (method & SEMI_LAGRANGIAN_WITHDRAWAL_CONTINUOUS)
				&& op == ProgramOperation::STOCHASTIC) {
			cerr << "error: querying the stochastic control for the"
					" semi-Lagrangian scheme is not"
					" implemented" << endl;
		}

		if( (method & EXPLICIT_IMPULSE)
				&& op == ProgramOperation::IMPULSE) {
			cerr << "error: querying the impulse control for the"
					" explicit impulse scheme is not"
					" implemented" << endl;
		}
	}

	// Finite-horizon or infinite-horizon?
	const bool finite_horizon = T > 0.;

	Real previousValue = nan(""), previousChange = nan("");

	RectilinearGrid1 initialGrid(
		// Uniform grid
		/*
		Axis::uniform(
			-boundary,
			#ifdef QUANT_PDE_FEX_RATE_EVEN
			m,
			#else
			+boundary,
			#endif
			gridPoints
		)
		*/

		// Clustered grid
		Axis::cluster(
			-boundary,  // Left-hand boundary
			m,         // Feature to cluster around
			#ifdef QUANT_PDE_FEX_RATE_EVENT
			m,
			#else
			+boundary,  // Right-hand boundary
			#endif
			gridPoints, // Number of points
			10.         // Clustering intensity
		)
	);

	//cerr << initialGrid << endl;

	// Table headers
	if(op != ProgramOperation::CONVERGENCE_TEST) {
		Rmin = Rmax;
	} else {
		cout
			<< setw(td) << "Nodes"         << "\t"
			<< setw(td) << "Control Nodes" << "\t"
		;

		if(finite_horizon) { cout << setw(td) << "Timesteps" << "\t"; }

		cout
			<< setw(td) << "Value at x=0"  << "\t"
			<< setw(td) << "Iterations"    << "\t"
			<< setw(td) << "Change"        << "\t"
			<< setw(td) << "Ratio"
			<< endl
		;
	}

	for(
		int ref = 0;
		ref <= Rmax;
		++ref, controlPoints*=2, TN*=2
	) {

		if(ref < Rmin) { continue; }

		// Refine the grid ref times
		auto grid = initialGrid.refined( ref );

		RectilinearGrid1 stochasticControls(
			Axis::uniform(
				q_min,
				q_max,
				controlPoints
			)
		);

		ToleranceIteration toleranceIteration;

		////////////////////////////////////////////////////////////////
		// Linear system tree
		////////////////////////////////////////////////////////////////

		Discretizee discretizee(
			grid, // Domain
			F,    // Effect of interest rate differential on the
			      // exchange rate
			M,    // Cost of interest rate differential
			N,    // Cost of exchange rate differential
			rho,  // Discount factor
			v,    // Volatility
			m,    // Central parity

			// Is the q control handled implicitly?
			!( method & SEMI_LAGRANGIAN_WITHDRAWAL_CONTINUOUS )
		);

		// Policy iteration for stochastic control
		MinPolicyIteration1_1 stochasticPolicy(
			grid,               // Domain
			stochasticControls, // Control grid
			discretizee         // Lu - rho*u + f
		);
		stochasticPolicy.setIteration(toleranceIteration);

		// Policy iteration for impulse control
		Impulse1_1 impulseWithdrawal(
			grid,             // Spatial grid
			interventionCost, // Cost of an intervention
			xplus             // State of exchange rate after
			                  // intervention
		);
		MinPolicyIteration1_1 impulsePolicy(
			grid,             // Domain
			grid,             // Control grid
			impulseWithdrawal // Impulse
		);
		impulsePolicy.setIteration(toleranceIteration);

		////////////////////////////////////////////////////////////////
		// Stepper
		////////////////////////////////////////////////////////////////

		typedef ReverseBDFOne1 Discretization;

		ReverseConstantStepper stepper(
			0.,   // Initial time
			T,    // Expiry time
			T/TN  // Timestep size
		);

		if(method != EXPLICIT) {
			stepper.setInnerIteration(toleranceIteration);
		}

		// What to discretize
		LinearSystem *discretize;
		if(method & SEMI_LAGRANGIAN_WITHDRAWAL_CONTINUOUS) {
			discretize = &discretizee;
		} else {
			discretize = &stochasticPolicy;
		}

		// Discretize
		Discretization discretization(grid, *discretize);
		discretization.setIteration(stepper);

		////////////////////////////////////////////////////////////////
		// Horizon selection
		////////////////////////////////////////////////////////////////

		Iteration *iteration;
		IterationNode *penalized;

		if(finite_horizon) {
			iteration = &stepper;
			penalized = &discretization;
		} else {
			iteration = &toleranceIteration;
			penalized = &stochasticPolicy;
		}

		////////////////////////////////////////////////////////////////

		// Penalty method
		PenaltyMethod penalty(
			grid,         // Domain
			*penalized,   // Stochastic control component
			impulsePolicy // Impulse control component
		);
		penalty.setIteration(toleranceIteration);

		IterationNode *root;
		if(method & EXPLICIT_IMPULSE) {
			// No impulse root
			root = penalized;
		} else {
			// Impulse root
			root = &penalty;
		}

		////////////////////////////////////////////////////////////////
		// Exercise events
		////////////////////////////////////////////////////////////////

		auto setDifferential = [=,&stepper] (
			const Interpolant1 &u,
			Real x
		) {

			Real best = -numeric_limits<Real>::infinity();

			const Real dt = stepper.timestep();

			if(method & SEMI_LAGRANGIAN_WITHDRAWAL_CONTINUOUS) {
				const Real dq = (q_max-q_min)/(controlPoints-1);
				for(int i = 0; i < controlPoints; ++i) {
					const Real q = q_min + i * dq;
					const Real newValue =
						u( x - F(q) * dt )
						- N(q) * dt
					;
					if(newValue > best) {
						best = newValue;
					}
				}
			}

			if(method & EXPLICIT_IMPULSE) {
				const Axis &X = grid[0];
				for(int i = 0; i < X.size(); ++i) {
					const Real newValue =
						u( X[i] )
						+ interventionCost(0., x, X[i])
					;
					if(newValue > best) {
						best = newValue;
					}
				}
			}

			return best;

		};

		if(method != IMPLICIT) {
			for(int m = 0; m < TN; ++m) {
				stepper.add(
					// Time at which the event takes place
					T / TN * m,

					// Domestic gov't sets differential
					setDifferential,

					// Spatial grid to interpolate on
					grid
				);
			}
		}

		////////////////////////////////////////////////////////////////
		// Run
		////////////////////////////////////////////////////////////////

		BiCGSTABSolver solver;
		auto u = iteration->solve(
			grid,                        // Domain
			[=] (Real x) { return 0.; }, // Payoff
			*root,                       // Root of linear system
			solver                       // Linear system solver
		);

		if(op == ProgramOperation::PLOT) {

			RectilinearGrid1 printGrid(
				Axis::uniform(
					-boundary,
					#ifdef QUANT_PDE_FEX_RATE_EVEN
					m,
					#else
					+boundary,
					#endif
					100
				)
			);

			cout << accessor( printGrid, u );

		} else if(op == ProgramOperation::STOCHASTIC) {

			// Print stochastic control

			Vector beta = discretizee.control(0);
			auto mask = penalty.constraintMask();
			for(int i = 0; i < grid.size(); i++) {
				if(mask[i]) {
					beta(i) = numeric_limits<Real>
							::infinity();
				}
			}

			for(Index i = 0; i < grid[0].size(); ++i) {
				cout << grid[0][i] << "\t" << beta(i) << endl;
			}

		} else if(op == ProgramOperation::IMPULSE) {

			// Print impulse control

			Vector beta = impulseWithdrawal.control(0);
			auto mask = penalty.constraintMask();
			for(int i = 0; i < grid.size(); i++) {
				if(!mask[i]) {
					beta(i) = numeric_limits<Real>
							::infinity();
				}
			}

			for(Index i = 0; i < grid[0].size(); ++i) {
				cout << grid[0][i] << "\t" << beta(i) << endl;
			}

		} else {

			////////////////////////////////////////////////////////
			// Print solution at x_0
			////////////////////////////////////////////////////////

			Real adjusted = x_0;
			#ifdef QUANT_PDE_FEX_RATE_EVEN
			// If x_0 > m, use the fact that u(.) is symmetric about
			// x=m
			if(x_0 > m) {
				adjusted = m - (x_0 - m);
			}
			#endif

			Real value = u( adjusted );
			Real
				change = value - previousValue,
				ratio = previousChange / change
			;
			previousValue = value;
			previousChange = change;

			cout.precision(12);

			auto its = toleranceIteration.iterations();
			const Real mean = accumulate(its.begin(),its.end(),0.)
					/ its.size();

			// Print row of table
			cout
				<< setw(td) << grid.size()               << "\t"
				<< setw(td) << stochasticControls.size() << "\t"
			;

			if(finite_horizon) { cout << setw(td) << TN << "\t"; }

			cout
				<< setw(td) << value  << "\t"
				<< setw(td) << mean   << "\t"
				<< setw(td) << change << "\t"
				<< setw(td) << ratio  << "\t"
				<< endl
			;

		}

	}

	return 0;

}

