////////////////////////////////////////////////////////////////////////////////
// gmwb.cpp
// --------
//
// Computes the price of a guaranteed minimum withdrawal benefit (GMWB) variable
// annuity (insurance contract) using several formulations.
//
// The ITERATED_OPTIMAL_STOPPING macro uses iterated optimal stopping to solve
// the GMWB problem. QuantPDE does not normally support iterated optimal
// stopping and hence the implementation is a bit of a hack, and will most
// likely not compile in future versions.
//
// Author: Parsiad Azimzadeh
////////////////////////////////////////////////////////////////////////////////

#ifdef ITERATED_OPTIMAL_STOPPING
	#define private public
	#define protected public
	#define QUANT_PDE_DO_EVENT_PUBLIC
#endif

#include <QuantPDE/Core>

#ifdef ITERATED_OPTIMAL_STOPPING
	#undef private
	#undef protected
#endif

#include <QuantPDE/Modules/Operators>

////////////////////////////////////////////////////////////////////////////////

#include <algorithm> // max, min
#include <chrono>    // duration
#include <climits>   // INT_MAX
#include <cmath>     // sqrt
#include <cstdlib>   // abs
#include <getopt.h>  // getopt_long
#include <iomanip>   // setw
#include <iostream>  // cout
#include <memory>    // unique_ptr
#include <numeric>   // accumulate
#include <thread>    // thread
#include <tuple>     // get

////////////////////////////////////////////////////////////////////////////////

using namespace QuantPDE;
using namespace QuantPDE::Modules;

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Methods
////////////////////////////////////////////////////////////////////////////////

constexpr int SEMI_LAGRANGIAN_WITHDRAWAL_CONTINUOUS = 1 << 0;
constexpr int SEMI_LAGRANGIAN_WITHDRAWAL_IMPULSE    = 1 << 1;

constexpr int EXPLICIT =
		  SEMI_LAGRANGIAN_WITHDRAWAL_IMPULSE
		| SEMI_LAGRANGIAN_WITHDRAWAL_CONTINUOUS;

constexpr int IMPLICIT = 0;

////////////////////////////////////////////////////////////////////////////////
// Options
////////////////////////////////////////////////////////////////////////////////

int method = 0;

Real T = 10.; // 14.28;
Real r = .05;
Real v = .2;

//Real alpha = 0.0138457882249;
Real alpha = 0.013891; // 0.036; // Hedging fee

Real G = 10.; // 7.; // Contract rate
Real kappa = 0.1; // Penalty rate

Real w0 = 100.; // Initial value of the account

int N = 32; // Initial number of timesteps
int M = 2; // Initial control set partition
int Mmax = INT_MAX; //16; // Maximum control set partition size

int Rmin = 0;
int Rmax = 6; // Maximum level of refinement

int realizedN = -1; // Realized timesteps
Real target = 1.; // Target relative error for variable timestepping

bool newton = false;
bool variable = false; // Variable stepping does not work for anything other
                       // than fully implicit!
bool quarter = false; // Quarter timestep on each refinement

////////////////////////////////////////////////////////////////////////////////
// Solution grid
////////////////////////////////////////////////////////////////////////////////

// Grid
RectilinearGrid2 grid(
	// Hand-picked axis, scaled by w0
	(w0 / 100.) * Axis {
		0., 5., 10., 15., 20., 25.,
		30., 35., 40., 45.,
		50., 55., 60., 65., 70., 72.5, 75., 77.5, 80., 82., 84.,
		86., 88., 90.,91., 92., 93., 94., 95.,
		96., 97., 98., 99., 100.,
		101., 102., 103., 104., 105., 106.,
		107., 108., 109., 110., 112., 114.,
		116., 118., 120., 123., 126.,
		130., 135., 140., 145., 150., 160., 175., 200., 225.,
		250., 300., 500., 750., 1000.
	},

	// 51 points evenly spaced on [0, w0]
	Axis::uniform(0., w0, 51)
);

// Automatically generated grid
/*
constexpr int points1 = 64;
constexpr int points2 = 50;

RectilinearGrid2 grid(
	Axis::cluster(0., 1000., points1, w0, 5.),
	Axis::cluster(0.,  100., points2, w0, 5.)
);
*/

////////////////////////////////////////////////////////////////////////////////
// Impulse functions
////////////////////////////////////////////////////////////////////////////////

/**
 * Penalized cashflow.
 */
inline Real penalizedCashflow(Real t, Real W, Real A, Real q) {
	return (1 - kappa) * q * A;
}

/**
 * Regular cashflow. Penalized only if the amount withdrawn exceeds Gdt.
 */
inline Real cashflow(Real t, Real W, Real A, Real q, Real Gdt) {
	return q * A - kappa * max( q * A - Gdt, 0. );
}

/**
 * State of the investment account after withdrawal.
 */
inline Real Wplus(Real t, Real W, Real A, Real q) {
	return max( W - q * A, 0. );
}

/**
 * State of the withdrawal account after withdrawal.
 */
inline Real Aplus(Real t, Real W, Real A, Real q) {
	return A - q * A;
}

/**
 * Value of the contract after withdrawal.
 */
inline Real Vminus(const Interpolant2 &V, Real t, Real W, Real A, Real q) {
	return V( Wplus(t, W, A, q), Aplus(t, W, A, q) );
}

////////////////////////////////////////////////////////////////////////////////
// Operator to discretize
// i.e. Lu + f, where L is the infinitesimal generator of process (W_t, A_t)
////////////////////////////////////////////////////////////////////////////////

class Discretizee final : public RawControlledLinearSystem2_1 {

	const RectilinearGrid2 &grid;
	const Real r, v, q;

	const bool controlled;
	const Vector zero;

public:

	template <typename G1>
	Discretizee(
		G1 &grid,
		Real interest,
		Real volatility,
		Real dividends,
		bool controlled = true
	) noexcept :
		grid( grid ),
		r( interest ),
		v( volatility ),
		q( dividends ),
		controlled( controlled ),
		zero( grid.zero() )
	{
	}

	virtual Matrix A(Real) {
		Matrix M = grid.matrix();
		M.reserve( IntegerVector::Constant(grid.size(), 4) );

		// Axes
		const Axis &W = grid[0];
		const Axis &A = grid[1];

		// Control as a vector
		Index k = 0;
		const Vector &raw = controlled ? control(0) : zero;

		// A
		for(Index j = 0; j < A.size(); ++j) {

			// W = 0
			if(j > 0) {
				const Real dA  = A[j] - A[j-1];
				const Real tmp = raw(k) / dA;

				M.insert(k, k           ) =  tmp + r;
				M.insert(k, k - W.size()) = -tmp;
			} else {
				M.insert(k, k) = r;
			}
			++k;

			// 0 < W < W_max
			for(Index i = 1; i < W.size() - 1; ++i) {

				// Deltas
				const Real
					dWb = W[i  ] - W[i-1],
					dWc = W[i+1] - W[i-1],
					dWf = W[i+1] - W[i  ]
				;

				const Real tmp1 = v * v * W[i] * W[i];
				const Real tmp2 = (r - q) * W[i] - raw(k);

				const Real alpha_common = tmp1 / dWb / dWc ;
				const Real  beta_common = tmp1 / dWf / dWc ;

				// Central
				Real alpha_i = alpha_common - tmp2 / dWc ;
				Real  beta_i =  beta_common + tmp2 / dWc ;
				if(alpha_i < 0.) {
					// Forward
					alpha_i = alpha_common;
					 beta_i = beta_common + tmp2 / dWf ;
				} else if(beta_i < 0.) {
					// Backward
					alpha_i = alpha_common - tmp2 / dWb;
					 beta_i =  beta_common;
				}

				M.insert(k, k - 1) = -alpha_i;
				M.insert(k, k + 1) = - beta_i;

				// Branching is slow, but peanuts compared to
				// the other work that we have to do
				const Real base = alpha_i + beta_i + r;
				if(j > 0) {
					const Real dA  = A[j] - A[j-1];
					const Real tmp = raw(k) / dA;

					M.insert(k, k           ) =  tmp + base;
					M.insert(k, k - W.size()) = -tmp;
				} else {
					M.insert(k, k) = base;
				}

				++k;

			}

			// W = W_max
			if(j > 0) {
				const Real dA  = A[j] - A[j-1];
				const Real tmp = raw(k) / dA;

				M.insert(k, k           ) =  tmp + q;
				M.insert(k, k - W.size()) = -tmp;
			} else {
				M.insert(k, k) = q;
			}
			++k;

		}

		M.makeCompressed();
		return M;
	}

	virtual Vector b(Real t) {
		Vector b = grid.vector();

		const Axis &W = grid[0];
		const Axis &A = grid[1];

		// Control as a vector
		Index k = 0;
		const Vector &raw = controlled ? control(0) : zero;

		// A = 0 (no withdrawal)
		for(Index i = 0; i < W.size(); ++i) {
			b(k) = 0.;
			++k;
		}

		//const Real Wmax = W[ W.size() - 1 ];

		// A > 0
		for(Index j = 1; j < A.size(); ++j) {
			// 0 <= W < W_max
			for(Index i = 0; i < W.size() - 1; ++i) {
				b(k) = raw(k);
				++k;
			}

			// W = W_max
			b(k) = raw(k);
			++k;
		}

		return b;
	}

	virtual bool isATheSame() const {
		return !controlled;
	}

};

////////////////////////////////////////////////////////////////////////////////

constexpr int progressWidth = 70;
constexpr int progressSleep = 100;

std::tuple<Real, Real, Real, int> solve(Real alpha) {

	////////////////////////////////////////////////////////////////////////
	// Iteration tree
	////////////////////////////////////////////////////////////////////////

	unique_ptr<ReverseTimeIteration> stepper;

	if(variable) {
		stepper = unique_ptr<ReverseTimeIteration>(
			(ReverseTimeIteration *)
			new ReverseVariableStepper(
				0.,    // Initial time
				T,     // Expiry time
				T / N, // Initial timestep size
				target // Target error
			)
		);
	} else {
		stepper = unique_ptr<ReverseTimeIteration>(
			(ReverseTimeIteration *)
			new ReverseConstantStepper(
				0.,    // Initial time
				T,     // Expiry time
				T / N  // Timestep size
			)
		);
	}

	// Tolerance iteration
	ToleranceIteration toleranceIteration;
	if(method != EXPLICIT) {
		stepper->setInnerIteration(toleranceIteration);
	}

	////////////////////////////////////////////////////////////////////////
	// Linear system tree
	////////////////////////////////////////////////////////////////////////

	typedef ReverseBDFOne2 Discretization;

	////////////////////////////////////////////////////////////////////////

	// Black-Scholes
	Discretizee discretizee(
		grid, r, v, alpha,
		!( method & SEMI_LAGRANGIAN_WITHDRAWAL_CONTINUOUS )
	);

	RectilinearGrid1 stochasticControls( Axis { 0., G } );
	MinPolicyIteration2_1 stochasticPolicy(
		grid,
		stochasticControls,
		discretizee
	);
	stochasticPolicy.setIteration(toleranceIteration);

	// What to discretize
	LinearSystem *discretize;
	if(method & SEMI_LAGRANGIAN_WITHDRAWAL_CONTINUOUS) {
		discretize = &discretizee;
	} else {
		discretize = &stochasticPolicy;
	}

	// Discretize
	Discretization discretization(grid, *discretize);
	discretization.setIteration(*stepper);

	// Impulse withdrawal
	RectilinearGrid1 impulseControls( Axis::range( 1. / M, 1. / M, 1. ) );
	Impulse2_1 impulseWithdrawal(
		grid,

		penalizedCashflow, // Cash received from impulse withdrawal

		Wplus,             // State of investment account after impulse
		                   // withdrawal

		Aplus              // State of withdrawal account after impulse
		                   // withdrawal
	);
	MinPolicyIteration2_1 impulsePolicy(
		grid,
		impulseControls,
		impulseWithdrawal
	);
	impulsePolicy.setIteration(toleranceIteration);

	// Penalty method
	PenaltyMethod penalty(grid, discretization, impulsePolicy);
	penalty.setIteration(toleranceIteration);

	// Root
	IterationNode *root;
	if(method & SEMI_LAGRANGIAN_WITHDRAWAL_IMPULSE) {
		// No impulse root
		root = &discretization;
	} else {
		// Impulse root
		root = &penalty;
	}

	////////////////////////////////////////////////////////////////////////
	// Exercise events
	////////////////////////////////////////////////////////////////////////

	auto withdrawal = [=,&stepper] (const Interpolant2 &V, Real W, Real A) {

		// No withdrawal
		Real best = V(W, A);

		const Real dt = stepper->timestep();

		// Contract withdrawal amount
		const Real Gdt = G * dt;

		if(method & SEMI_LAGRANGIAN_WITHDRAWAL_CONTINUOUS) {
			// Nonpenalty

			const Real beta = min(A, Gdt);

			for(int i = 1; i <= M; ++i) {
				const Real gamma = beta * i/M;
				//const Real gamma = beta;

				const Real newValue =
					 Vminus(V, 0., W, A, gamma / A     )
					+ cashflow(0., W, A, gamma / A, Gdt)
				;

				if(newValue > best) {
					best = newValue;
				}
			}

		}

		// Penalty
		if(
			(method & SEMI_LAGRANGIAN_WITHDRAWAL_IMPULSE)
			&& A > Gdt
		) {

			for(int i = 1; i <= M; ++i) {

				const Real gamma = Gdt + (A-Gdt) * i/M;

				const Real newValue =
					 Vminus(V, 0., W, A, gamma / A     )
					+ cashflow(0., W, A, gamma / A, Gdt)
				;

				if(newValue > best) {
					best = newValue;
				}

			}

		}

		return best;

	};

	if(method != IMPLICIT) {
		for(int m = 0; m < N; ++m) {
			stepper->add(
				// Time at which the event takes place
				T / N * m,

				withdrawal,

				// Spatial grid to interpolate on
				grid
			);
		}
	}

	////////////////////////////////////////////////////////////////////////
	// Payoff
	////////////////////////////////////////////////////////////////////////

	Function2 payoff = [=] (Real W, Real A) {
		return max(W, (1 - kappa) * A);
	};

	////////////////////////////////////////////////////////////////////////
	// Running
	////////////////////////////////////////////////////////////////////////

	Real value;
	BiCGSTABSolver solver;

	#ifdef ITERATED_OPTIMAL_STOPPING
	{ // Start iterated optimal stopping test

		// Solutions at each time
		Vector *current  = new Vector[N+1];
		Vector *previous = new Vector[N+1];

		// Solution at the expiry
		previous[0] = current[0] = grid.image(payoff);

		// Initialize degenerate circular buffers
		toleranceIteration.history = new Iteration::CB(1);
		stepper          ->history = new Iteration::CB(1);

		// Tolerance loop
		Event<2> event(withdrawal, grid);
		toleranceIteration.its.push_back(0);
		bool converged, first = true;
		do {
			// Timestep loop
			converged = true;
			for(int n = 0; n < N; ++n) {
				// Explicit and implicit times
				const Real texp = T - T *  n    / N;
				const Real t    = T - T * (n+1) / N;

				// What to use as the previous iterand?
				const Vector *const initial =
					first ?
						  &current [n  ]
						: &previous[n+1]
				;
				const Real initialTime = first ? texp : t;

				// Previous iterand
				toleranceIteration.history->clear();
				toleranceIteration.history->push( make_tuple(
					initialTime,
					*initial
				));
				stepper->history->clear();
				stepper->history->push( make_tuple(
					texp,
					current[n]
				));

				toleranceIteration.implicitTime = t;
				stepper          ->implicitTime = t;

				// Start nodes
				toleranceIteration.startNodes();
				stepper          ->startNodes();

				// Solve Ax=b
				((LinearSolver *) &solver)->initialize(
					root->A(t)
				);
				current[n+1] = solver.solve(
					root->b(t),
					*initial
				);

				// End nodes
				stepper          ->endNodes();
				toleranceIteration.endNodes();

				if(method != IMPLICIT) {
					// Apply event
					current[n+1] = event.doEvent(
							current[n+1]);
				}

				if(!first && converged) {
					// Compare current[n + 1] and
					// previous[n + 1]

					const Vector *const a = &current [n+1];
					const Vector *const b = &previous[n+1];

					const Real tmp = relativeError(*a, *b);

					if(tmp > QuantPDE::tolerance) {
						converged = false;
					}
				}
			}

			// No longer on the first outer iteration
			if(first) {
				converged = false;
				first = false;
			}

			// Swap solutions
			Vector *const tmp = current;
			current  = previous;
			previous = tmp;

			// Increment number of iterations
			++toleranceIteration.its.back();
		} while(!converged);

		// Linear interpolate to get V(w0, w0)
		PiecewiseLinear2 V(grid, previous[N]);
		value = V(w0, w0);

		// Housekeeping

		delete toleranceIteration.history;
		toleranceIteration.history = nullptr;

		delete stepper->history;
		stepper->history = nullptr;

		delete [] current;
		current = nullptr;

		delete [] previous;
		previous = nullptr;

	} // End iterated optimal stopping test
	#else
	{ // Start policy iteration test

		#ifdef _OPENMP
		// Progress bar
		thread progress([&] () {
			chrono::milliseconds duration( progressSleep );

			while(stepper->nextTime() < 0.) {
				// Idle state

				this_thread::sleep_for( duration );
			}

			while(stepper->nextTime() != 0.) {
				// Busy state

				const Real progress = (T-stepper->nextTime())/T;

				cerr << "[";
				int pos = progressWidth * progress;

				for(int i = 0; i < pos; ++i) {
					cerr << "=";
				}

				if(pos < progressWidth) {
					cerr << ">";
				}

				for(int i = pos + 1; i < progressWidth; ++i) {
					cerr << " ";
				}

				cerr << "] " << int(progress * 100.) << "%  \r";
				cerr.flush();

				this_thread::sleep_for( duration );
			}
		});
		#endif

		auto V = stepper->solve(
			grid,   // Domain
			payoff, // Initial condition
			*root,  // Root of linear system tree
			solver  // Linear system solver
		);

		#ifdef _OPENMP
		progress.join();
		#endif

		value = V(w0, w0);

	} // End policy iteration test
	#endif

	// Actual number of timesteps
	#ifdef ITERATED_OPTIMAL_STOPPING
	realizedN = N;
	#else
	realizedN = stepper->iterations().back();
	#endif

	////////////////////////////////////////////////////////////////////////
	// Statistics
	////////////////////////////////////////////////////////////////////////

	Real mean = 1., var = 0.;
	int max = 1;

	if(method != EXPLICIT) {
		auto its = toleranceIteration.iterations();

		mean = accumulate(its.begin(),its.end(),0.)/its.size();

		var = 0.;
		for(auto x : its) { var += (x - mean) * (x - mean); }

		max = ( *max_element(its.begin(), its.end()) );
	}

	return make_tuple( value, mean, var, max );

}

////////////////////////////////////////////////////////////////////////////////

constexpr int td = 20;

void printHeaders() {
	cout
		<< setw(td) << "Nodes"           << "\t"
		<< setw(td) << "Control Nodes"   << "\t"
		<< setw(td) << "Time Steps"      << "\t"
		<< setw(td) << "Value"           << "\t"
		<< setw(td) << "Mean Iterations" << "\t"
		<< setw(td) << "Std Iterations"  << "\t"
		<< setw(td) << "Max Iterations"  << "\t"
		<< setw(td) << "Change"          << "\t"
	;

	if(!newton) {
		cout << setw(td) << "Ratio";
	}

	cout << setw(td) << "Seconds";

	cout << endl;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {

	#ifdef ITERATED_OPTIMAL_STOPPING
	cerr << "warning: compiled with iterated optimal stopping" << endl
			<< endl;
	#endif

	////////////////////////////////////////////////////////////////////////
	// Options
	////////////////////////////////////////////////////////////////////////

	{
		// Long option names
		static struct option opts[] = {
			{ "semi-implicit"   , no_argument,       0, 0 },
			{ "explicit"        , no_argument,       0, 0 },
			{ "newton"          , no_argument,       0, 0 },
			{ "variable"        , no_argument,       0, 0 },
			{ "quarter-timestep", no_argument,       0, 0 },
			{ "fair-fee"        , required_argument, 0, 0 },
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
							newton = true;
							break;
						case 3:
							variable = true;
							break;
						case 4:
							quarter = true;
							break;
						case 5:
							alpha = atof(optarg);
							break;
						default:
							break;
					}
				break;
			}
		}

		if(variable && method != IMPLICIT) {
			cerr << "error: variable timestepping cannot be used "
					"with a nonimplicit method" << endl;
			return 1;
		}

		#ifdef ITERATED_OPTIMAL_STOPPING
		if(method == SEMI_LAGRANGIAN_WITHDRAWAL_CONTINUOUS) {
			cerr << "error: cannot use iterated optimal stopping "
					"with semi-implicit method" << endl;
			return 1;
		}
		#endif
	}

	////////////////////////////////////////////////////////////////////////
	// Table headers
	////////////////////////////////////////////////////////////////////////

	cout.precision(12);
	Real previousValue = nan(""), previousChange = nan("");

	if(!newton) {
		printHeaders();
	}

	////////////////////////////////////////////////////////////////////////
	// Refinement loop
	////////////////////////////////////////////////////////////////////////

	for(
		int l = 0;
		l <= Rmax;
		++l, N *= (quarter ? 4 : 2), M *= 2, target /= 2.
	) {

		if( l < Rmin ) {
			// Refine grid
			grid = grid.refined();
			continue;
		}

		M = min(M, Mmax);

		////////////////////////////////////////////////////////////////
		// Outermost Newton iteration to find fair fee
		////////////////////////////////////////////////////////////////

		Real value, mean, var;
		int max;

		auto start = chrono::steady_clock::now();

		if(newton) {
			// Spacing
			cout << endl;

			// Newton iteration
			cout
				<< setw(td) << "\t"
				<< setw(td) << "Fair Fee"   << "\t"
				<< setw(td) << "Value"      << "\t"
				<< setw(td) << "Residual"   << "\t"
				<< endl
			;

			while(true) {
				// f(alpha)
				std::tie(value, mean, var, max) = solve(alpha);
				const Real f0 = value - w0;

				cout
					<< setw(td) << "\t"
					<< setw(td) << alpha << "\t"
					<< setw(td) << value << "\t"
					<< setw(td) << f0    << "\t"
					<< endl
				;

				if(abs(f0) < tolerance) {
					break;
				}

				// f(alpha + epsilon)
				auto tmp = solve(alpha + epsilon);
				const Real f1 = std::get<0>(tmp) - w0;

				// f'(alpha)
				const Real fp = (f1 - f0) / epsilon;

				// Next iterand
				alpha -= f0 / fp;
			}

			// Spacing
			cout << endl;

			// Print headers
			printHeaders();
		} else {
			// No Newton iteration
			std::tie(value, mean, var, max) = solve(alpha);
		}

		auto end = chrono::steady_clock::now();
		auto diff = end - start;

		////////////////////////////////////////////////////////////////
		// Print table rows
		////////////////////////////////////////////////////////////////

		Real
			change = value - previousValue,
			ratio = previousChange / change
		;

		int tmp = M;
		if(method == EXPLICIT) {
			tmp *= 2;
		}

		cout
			<< setw(td) << grid.size() << "\t"
			<< setw(td) << tmp         << "\t"
			<< setw(td) << realizedN   << "\t"
			<< setw(td) << value       << "\t"
			<< setw(td) << mean        << "\t"
			<< setw(td) << sqrt(var)   << "\t"
			<< setw(td) << max         << "\t"
			<< setw(td) << change      << "\t"
		;

		if(!newton) {
			cout << setw(td) << ratio << "\t";
		}

		cout
			<< setw(td) << chrono::duration <double> (diff).count()
			<< endl;

		previousChange = change;
		previousValue = value;

		// Refine grid
		grid = grid.refined();

	}

	return 0;
}

#if 0

////////////////////////////////////////////////////////////////////////////////
// Legacy (no longer used)
////////////////////////////////////////////////////////////////////////////////

class ImpulseWithdrawal final : public RawControlledLinearSystem2_1 {

	const RectilinearGrid2 &grid;
	Noncontrollable2 kappa;

public:

	template <typename G, typename F1>
	ImpulseWithdrawal(
		G &grid,
		F1 &&kappa
	) noexcept :
		grid(grid),
		kappa(kappa)
	{
	}

	virtual Matrix A(Real t) {
		Matrix M = grid.matrix();
		M.reserve(IntegerVector::Constant(grid.size(), 4));

		Index k = 0;
		const Vector &raw = control(0);

		for(auto node : grid) {
			const Real W = node[0]; // Investment
			const Real A = node[1]; // Withdrawal

			// Amount withdrawn pre-penalty
			const Real gamma = raw(k) * A;

			const Real Wplus = max(W - gamma, 0.);
			const Real Aplus = A - gamma;

			// Interpolation data
			auto data = linearInterpolationData(grid,
					{{Wplus,Aplus}});

			const Index i0 = std::get<0>( data[0] );
			const Index i1 = std::get<0>( data[1] );
			const Real  w0 = std::get<1>( data[0] );
			const Real  w1 = std::get<1>( data[1] );

			assert( (grid[0][i0+1] - Wplus)
					/ (grid[0][i0+1] - grid[0][i0]) == w0 );
			assert( (grid[1][i1+1] - Aplus)
					/ (grid[1][i1+1] - grid[1][i1]) == w1 );

			const Index j = grid.index(i0, i1);

			M.insert(k, j                     ) =    w0  *    w1 ;
			M.insert(k, j     + grid[0].size()) =    w0  * (1-w1);
			M.insert(k, j + 1                 ) = (1-w0) *    w1 ;
			M.insert(k, j + 1 + grid[0].size()) = (1-w0) * (1-w1);

			++k;
		}

		M.makeCompressed();
		return grid.identity() - M;
	}

	virtual Vector b(Real t) {
		Vector b = grid.vector();

		Index k = 0;
		const Vector &raw = control(0);

		for(auto node : accessor(grid, b)) {
			auto WA = &node;
			const Real W = WA[0]; // Investment
			const Real A = WA[1]; // Withdrawal

			// Amount withdrawn, pre-penalty
			const Real gamma = raw(k) * A;

			// Cashflow minus adjustment
			*node = (1 - kappa(t, W, A)) * gamma;

			++k;
		}

		return b;
	}

};

////////////////////////////////////////////////////////////////////////////////

// This results in a nonomonotone scheme and should not be used
// (testing purposes only)

class ContinuousWithdrawal final : public RawControlledLinearSystem2_1 {

	const RectilinearGrid2 &grid;
	Noncontrollable2 contractRate;

public:

	template <typename G, typename F1>
	ContinuousWithdrawal(
		G &grid,
		F1 &&contractRate
	) noexcept :
		grid(grid),
		contractRate(contractRate)
	{
	}

	virtual Matrix A(Real t) {
		Matrix M = grid.matrix();
		M.reserve(IntegerVector::Constant(grid.size(), 4));

		const Axis &W = grid[0];
		const Axis &A = grid[1];

		// Control as a vector
		Index k = W.size();
		const Vector &raw = control(0);

		// A > 0
		for(Index j = 1; j < A.size(); ++j) {
			// W = 0
			{
				const Real G  = contractRate(t, W[0], A[j]);
				const Real g = raw(k) * G;

				const Real tA = g / (A[j] - A[j-1]);

				M.insert(k, k           ) = +tA;
				M.insert(k, k - W.size()) = -tA;

				++k;
			}

			/*
			// W > 0
			for(Index i = 1; i < W.size(); ++i) {
				const Real G  = contractRate(t, W[i], A[j]);
				const Real g = raw(k) * G;

				const Real tA = g / (A[j] - A[j-1]);
				const Real tW = g / (W[i] - W[i-1]);

				M.insert(k, k           ) = + tA + tW;
				M.insert(k, k - W.size()) = - tA     ;
				M.insert(k, k - 1       ) =      - tW;

				++k;
			}
			*/

			//#if 0
			// 0 < W < W_max
			for(Index i = 1; i < W.size() - 1; ++i) {
				const Real G  = contractRate(t, W[i], A[j]);
				const Real g = raw(k) * G;

				const Real tA = g / (A[j  ] - A[j-1]);
				const Real tW = g / (W[i+1] - W[i-1]);

				M.insert(k, k + 1       ) =      + tW;
				M.insert(k, k           ) = + tA     ;
				M.insert(k, k - W.size()) = - tA     ;
				M.insert(k, k - 1       ) =      - tW;

				++k;
			}

			// W = W_max
			{
				const Index i = W.size() - 1;

				const Real G  = contractRate(t, W[i], W[j]);
				const Real g = raw(k) * G;

				const Real tA = g / (A[j] - A[j-1]);
				const Real tW = g / (W[i] - W[i-1]);

				M.insert(k, k           ) = + tA + tW;
				M.insert(k, k - W.size()) = - tA     ;
				M.insert(k, k - 1       ) =      - tW;

				++k;
			}
			//#endif
		}

		M.makeCompressed();
		return M;
	}

	virtual Vector b(Real t) {
		Vector b = grid.vector();

		const Axis &W = grid[0];
		const Axis &A = grid[1];

		// Control as a vector
		Index k = 0;
		const Vector &raw = control(0);

		// A = 0 (no withdrawal)
		for(Index i = 0; i < W.size(); ++i) {
			b(k) = 0.;
			++k;
		}

		// A > 0
		for(Index j = 1; j < A.size(); ++j) {
			// W >= 0
			for(Index i = 0; i < W.size(); ++i) {
				const Real G = contractRate(t, W[i], A[j]);
				const Real g = raw(k) * G;
				b(k) = g;
				++k;
			}
		}

		return b;
	}

};

////////////////////////////////////////////////////////////////////////////////

// Use this in lieu of linear interpolation for a smoother convergence rate for
// the first few iterates

class WithdrawalEvent final : public EventBase {

	const RectilinearGrid2 &grid;
	const Real Gdt, kappa;

	template <typename V>
	Vector _doEvent(V &&Vp) const {
		// Interpolant
		PiecewiseLinear2 Vplus(grid, Vp);

		// Solution before withdrawal
		Vector Vm = grid.vector();

		// Axes
		const Axis &W = grid[0];
		const Axis &A = grid[1];

		assert(W[0] == 0.);
		assert(A[0] == 0.);

		// Minimum withdrawal amount (exclusive)
		Real L;
		if(method & SEMI_LAGRANGIAN_WITHDRAWAL_CONTINUOUS) {
			// Continuous withdrawal is handled here
			L = 0.;
		} else {
			// Continuous withdrawal is not handled here
			L = Gdt;
		}

		Index k = 0;
		for(Index j = 0; j < A.size(); ++j) {

			// Maximum withdrawal amount (inclusive)
			Real U;
			if(method & SEMI_LAGRANGIAN_WITHDRAWAL_IMPULSE) {
				// Impulse withdrawal is handled here
				U = A[j];
			} else {
				// Impulse withdrawal is not handled here
				U = min(A[j], Gdt);
			}

			for(Index i = 0; i < W.size(); ++i) {

				// Find the optimal control at this node

				// No withdrawal
				Real best = Vp(k);

				// Withdrawal at the contract rate
				if(
					method
					& SEMI_LAGRANGIAN_WITHDRAWAL_CONTINUOUS
				) {
					const Real beta = min(A[j], Gdt);

					const Real newValue =
							Vminus(
								Vplus,
								0.,
								W[i], A[j],
								beta / A[j]
							)
							+ cashflow(
								0.,
								W[i], A[j],
								beta / A[j],
								Gdt
							)
					;

					if(newValue > best) {
						best = newValue;
					}
				}

				// Full surrender at penalty
				if(
					(A[j] > Gdt) && (method
					& SEMI_LAGRANGIAN_WITHDRAWAL_IMPULSE)
				) {
					const Real newValue =
							Vminus(
								Vplus,
								0.,
								W[i], A[j],
								1.
							)
							+ cashflow(
								0.,
								W[i], A[j],
								1.,
								Gdt
							)
					;

					if(newValue > best) {
						best = newValue;
					}
				}

				// Withdraw <= W[i] (along W axis)
				for(Index ii = i; ii >= 0; --ii) {
					// Amount withdrawn (pre-penalty)
					const Real gamma = W[i] - W[ii];

					// Skip anything outside the bounds
					// TODO: Binary search to find starting
					//       point instead
					if( gamma <= L ) { continue; }
					if( gamma >  U ) {    break; }

					// Interpolate on the A axis
					const Real Ap = Aplus(
						0.,
						W[i], A[j],
						gamma / A[j]
					);
					auto data = linearInterpolationData(A,
							Ap);

					const Index jj = get<0>(data);
					const Real   w = get<1>(data);

					const Real newValue =
						     w *Vp(ii+ jj   *W.size())
						+ (1-w)*Vp(ii+(jj+1)*W.size())
						+ cashflow(
							0.,
							W[ii], A[jj],
							gamma / A[jj],
							Gdt
						);
					;

					if(newValue > best) {
						best = newValue;
					}
				}

				// Withdraw > W[i]
				/*
				if(A[j] > W[i]) {
					const Real Alast = A[j] - W[i];

					for(Index jj = 0; A[jj] < Alast; ++jj) {
						// Amount withdrawn
						const Real gamma = A[j] - A[jj];

						const Real newValue =
							  Vp(jj * W.size())
							+ cashflow(
								0.,
								W[ii], A[jj]
								gamma / A[jj],
								Gdt
							)
						;

						if(newValue > best) {
							best = newValue;
						}
					}
				}
				*/

				// Withdraw <= A[j] (along A axis)
				for(Index jj = j; jj >= 0; --jj) {
					// Amount withdrawn pre-penalty
					const Real gamma = A[j] - A[jj];

					// Skip anything outside the bounds
					// TODO: Binary search to find starting
					//       point instead
					if( gamma <= L ) { continue; }
					if( gamma >  U ) {    break; }

					// Interpolate on the W axis
					const Real Wp = Wplus(
						0.,
						W[i], A[j],
						gamma / A[j]
					);
					auto data = linearInterpolationData(W,
							Wp);

					const Index ii = get<0>(data);
					const Real   w = get<1>(data);

					const Real newValue =
						     w *Vp( ii   +jj*W.size())
						+ (1-w)*Vp((ii+1)+jj*W.size())
						+ cashflow(
							0.,
							W[ii], A[jj],
							gamma / A[jj],
							Gdt
						);
					;

					if(newValue > best) {
						best = newValue;
					}
				}

				Vm(k++) = best;
			}
		}

		return Vm;
	}

	virtual Vector doEvent(const Vector &vector) const {
		return _doEvent(vector);
	}

	virtual Vector doEvent(Vector &&vector) const {
		return _doEvent(move(vector));
	}

public:

	template <typename G>
	WithdrawalEvent(
		G &&grid,
		Real Gdt,
		Real kappa
	) noexcept :
		grid(grid),
		Gdt(Gdt),
		kappa(kappa)
	{
	}

};

#endif
