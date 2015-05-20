////////////////////////////////////////////////////////////////////////////////
// optimal_consumption.cpp
// -----------------------
//
// The optimal consumption problem subject to transaction costs described in
// [1].
//
// [1] Chancelier, Jean-Philippe, Bernt Øksendal, and Agnès Sulem. "Combined
// stochastic control and optimal stopping, and application to numerical
// approximation of combined stochastic and impulse control." Proceedings of the
// Steklov Institute of Mathematics. VA Steklov. 237.0 (2002): 149-172.
//
// Author: Parsiad Azimzadeh
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>

////////////////////////////////////////////////////////////////////////////////

#include <cmath>    // fabs
#include <getopt.h>  // getopt_long
#include <iomanip>  // setw
#include <iostream> // cout, cerr
#include <limits>   // numeric_limits

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
// Parameters from Chancelier
////////////////////////////////////////////////////////////////////////////////

int method = 0;

// Stochastic control
Real q_max = 100.;
int stochasticControlPoints = 11;
int impulseControlPoints = 11;

// Volatility
Real v = 0.3;

// Risk-aversion
Real aver = 0.3;

// Drift
Real mu = 0.11;

// Interest rate
Real r = 0.07;

// Discount factor
Real rho = 0.1;

// Transaction costs
Real c = 0.05;
Real lambda = 0.1;

// Expiry
Real T = 10.;

// Initial value of stock and bank account
Real w_0 = 100.;
Real b_0 = 100.;

// Initial number of timesteps
int N = 32;

// Max/min refinement level
int Rmin = 0;
int Rmax = 6;

////////////////////////////////////////////////////////////////////////////////

/**
 * Controls the output of the program.
 */
enum class ProgramOperation {
	CONVERGENCE_TEST, /**< convergence test */
	PLOT, /**< outputs value function */
	PLOT_CONTROLS /**< outputs control data */
};
ProgramOperation op = ProgramOperation::CONVERGENCE_TEST;

////////////////////////////////////////////////////////////////////////////////
// Solution grid
////////////////////////////////////////////////////////////////////////////////

Axis axis {
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
};

////////////////////////////////////////////////////////////////////////////////
// Impulse functions
////////////////////////////////////////////////////////////////////////////////

/**
 * Amount moved to (result is positive) or from (result is negative) stock.
 */
inline Real zeta(Real W, Real B, Real zeta_frac) {
	// Withdrawal boundaries
	const Real lo = -W;
	const Real hi = (B - c)/(1 + lambda);

	if(lo > hi) {
		// No admissible strategies available
		return 0.;
	}

	// zeta > 0: Move money to stock
	// zeta < 0: Move money to bank
	return zeta_frac * (hi - lo) + lo;
}

/**
 * State of the stock after withdrawal.
 */
inline Real Wplus(Real t, Real W, Real B, Real zeta_frac) {
	const Real z = zeta(W, B, zeta_frac);
	return W + z;
}

/**
 * State of the bank account after withdrawal.
 */
inline Real Bplus(Real t, Real W, Real B, Real zeta_frac) {
	const Real z = zeta(W, B, zeta_frac);
	return B - z - c - lambda * fabs(z);
}

////////////////////////////////////////////////////////////////////////////////
// Operator to discretize
// i.e. Lu - rho + f, where L is the infinitesimal generator of process
//      (B_t, W_t)
////////////////////////////////////////////////////////////////////////////////

class Discretizee final : public RawControlledLinearSystem2_1 {

	const RectilinearGrid2 &grid;

	const Real r, v, mu, rho, aver;

	const bool controlled;
	const Vector zero;

public:

	template <typename G>
	Discretizee(
		G &grid,
		Real interest,
		Real volatility,
		Real drift,
		Real discount,
		Real aversion,
		bool controlled = true
	) noexcept :
		grid( grid ),
		r( interest ),
		v( volatility ),
		mu( drift ),
		rho( discount ),
		aver( aversion ),
		controlled( controlled ),
		zero( grid.zero() )
	{
	}

	virtual Matrix A(Real) {
		Matrix M = grid.matrix();
		M.reserve( IntegerVector::Constant(grid.size(), 4) );

		// Axes
		const Axis &W = grid[0]; const Index Ws = W.size();
		const Axis &B = grid[1]; const Index Bs = B.size();

		// Control as vector
		Index k = 0;
		const Vector &raw = controlled ? control(0) : zero;

		for(Index j = 0; j < Bs; ++j) {
			for(Index i = 0; i < Ws; ++i) {

				////////////////////////////////////////////////

				const Real
					dWb = W[i  ] - W[i-1],
					dWc = W[i+1] - W[i-1],
					dWf = W[i+1] - W[i  ]
				;

				// Neumann boundary condition set to zero as in
				// Chancelier's paper
				const Real W_i = W[i];

				const Real tmp1 = v * v * W_i * W_i;
				const Real tmp2 = mu * W_i;

				const Real alpha_common = tmp1 / dWb / dWc;
				const Real  beta_common = tmp1 / dWf / dWc;

				Real alpha_i = 0.;
				Real  beta_i = 0.;
				if(i > 0 && i < Ws - 1) {
					// 0 < W < W_max

					// Central
					alpha_i = alpha_common - tmp2 / dWc;
					 beta_i =  beta_common + tmp2 / dWc;

					if(alpha_i < 0.) {
						// Forward
						alpha_i = alpha_common;
						 beta_i =  beta_common
						 		+ tmp2 / dWf;
					} else if(beta_i < 0.) {
						// Backward
						alpha_i = alpha_common
								- tmp2 / dWb;
						 beta_i =  beta_common;
					}

					M.insert(k, k - 1) = -alpha_i;
					M.insert(k, k + 1) = - beta_i;
				}

				////////////////////////////////////////////////

				const Real B_j = B[j];

				const Real
					dBb = B[j  ] - B[j-1],
					dBf = B[j+1] - B[j  ]
				;

				// Upwind (cannot use central for monotonicity)
				Real alpha_j = 0.;
				Real  beta_j = 0.;
				if(j > 0 && j < Bs - 1) {
					// 0 < B < B_max

					const Real tmp = r * B_j - raw(k);
					if(tmp < 0.) {
						alpha_j = -tmp / dBb;
						M.insert(k, k - Ws) = -alpha_j;
					} else {
						 beta_j =  tmp / dBf;
						M.insert(k, k + Ws) = - beta_j;
					}
				}

				////////////////////////////////////////////////

				M.insert(k, k) =
					  alpha_i + beta_i
					+ alpha_j + beta_j
					+ rho;

				// Check positive coefficient condition
				assert(alpha_i >= 0.);
				assert( beta_i >= 0.);
				assert(alpha_j >= 0.);
				assert( beta_j >= 0.);

				++k;

			}
		}

		M.makeCompressed();
		return M;
	}

	virtual Vector b(Real) {

		Vector b = grid.vector();

		const Axis &W = grid[0];
		const Axis &B = grid[1];

		// Control as a vector
		Index k = 0;
		const Vector &raw = controlled ? control(0) : zero;

		// B = 0 (no consumption)
		for(Index i = 0; i < W.size(); ++i) {
			b(k) = 0.;
			++k;
		}

		// 0 < B < B_max
		for(Index j = 1; j < B.size(); ++j) {
			// 0 <= W <= W_max
			for(Index i = 0; i < W.size(); ++i) {
				b(k) = pow( raw(k) , aver ) / aver;
				++k;
			}
		}

		return b;
	}

	virtual bool isATheSame() const {
		return !controlled;
	}

};

// Spacing
constexpr int td = 20;

int main(int argc, char **argv) {

	////////////////////////////////////////////////////////////////////////
	// Options
	////////////////////////////////////////////////////////////////////////

	{
		// Long option names
		static struct option opts[] = {
			{ "mixed"           , no_argument,       0, 0 },
			{ "explicit"        , no_argument,       0, 0 },
			{ "min-refinement"  , required_argument, 0, 0 },
			{ "max-refinement"  , required_argument, 0, 0 },
			{ "plot-controls"   , no_argument,       0, 0 },
			{ "plot"            , no_argument,       0, 0 },
			{ nullptr           , 0,                 0, 0 }
		};

		int c;
		int index;
		while(
			(
				c = getopt_long(
					argc,
					argv,
					"v:",
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
							Rmin = atoi(optarg);
							break;
						case 3:
							Rmax = atoi(optarg);
							break;
						case 4:
							op = ProgramOperation::
								PLOT_CONTROLS;
							break;
						case 5:
							op = ProgramOperation::
									PLOT;
							break;
						default:
							break;
					}
				break;
				case 'v':
					v = atof(optarg);
					break;
				case ':':
				case '?':
					return 1;
			}
		}

		if(Rmin < 0 || Rmax < Rmin) {
			cerr << "error: the minimum level of refinement must be"
					" nonnegative and less than or equal to"
					" the maximum level of refinement"
					<< endl;
			return 1;
		}

		if(op == ProgramOperation::PLOT_CONTROLS && method
				== SEMI_LAGRANGIAN_WITHDRAWAL_CONTINUOUS) {
			cerr << "error: cannot output controls unless an "
					"implicit/explicit method is used"
					<< endl;
		}
	}

	const bool plot = (op == ProgramOperation::PLOT_CONTROLS)
			|| (op == ProgramOperation::PLOT);

	// Grid
	RectilinearGrid2 initialGrid( (w_0/100.) * axis, (b_0/100.) * axis );

	Real previousValue = nan(""), previousChange = nan("");

	if(plot) {
		Rmin = Rmax;
	} else {
		cout
			<< setw(td) << "Nodes"              << "\t"
			<< setw(td) << "q Control Nodes"    << "\t"
			<< setw(td) << "zeta Control Nodes" << "\t"
			<< setw(td) << "Timesteps"          << "\t"
			<< setw(td) << "Value"              << "\t"
		;

		if(method != EXPLICIT) {
			cout << setw(td) << "Iterations"    << "\t";
		}

		cout
			<< setw(td) << "Change"             << "\t"
			<< setw(td) << "Ratio"              << "\t"
			<< endl
		;
	}

	for(
		int ref = 0;
		ref <= Rmax;
		++ref, N *= 2,
		stochasticControlPoints *= 2, impulseControlPoints *= 2
	) {

		if(ref < Rmin) { continue; }

		// Refine the grid ref times
		auto grid = initialGrid.refined( ref );

		RectilinearGrid1 stochasticControls(
			Axis::uniform(
				0.,                     // Control lower bound
				q_max,                  // Control upper bound
				stochasticControlPoints // Number of nodes
			)
		);

		RectilinearGrid1 impulseControls(
			Axis::uniform(
				0.,                  // Control lower bound
				1.,                  // Control upper bound
				impulseControlPoints // Number of nodes
			)
		);

		ToleranceIteration toleranceIteration;

		////////////////////////////////////////////////////////////////
		// Linear system tree
		////////////////////////////////////////////////////////////////

		Discretizee discretizee(
			grid,  // Domain
			r,     // Interest rate
			v,     // Volatility
			mu,    // Drift
			rho,   // Discount factor
			aver,  // Risk-aversion

			// Is the q control handled implicitly?
			!( method & SEMI_LAGRANGIAN_WITHDRAWAL_CONTINUOUS )
		);

		// Policy iteration for stochastic control
		MinPolicyIteration2_1 stochasticPolicy(
			grid,               // Domain
			stochasticControls, // Control grid
			discretizee         // Lu - rho*u + f
		);
		stochasticPolicy.setIteration(toleranceIteration);

		// Policy iteration for impulse control
		Impulse2_1 impulseTransfer(
			// Spatial grid
			grid,

			// The intervention does not cost anything; the cost
			// is taken out of the bank account
			[] (Real, Real, Real, Real) { return 0.; },

			// State after an intervention
			Wplus, Bplus
		);
		MinPolicyIteration2_1 impulsePolicy(
			grid,            // Domain
			impulseControls, // Control grid
			impulseTransfer  // Impulse
		);
		impulsePolicy.setIteration(toleranceIteration);

		////////////////////////////////////////////////////////////////
		// Stepper
		////////////////////////////////////////////////////////////////

		typedef ReverseBDFOne2 Discretization;

		ReverseConstantStepper stepper(
			0.,  // Initial time
			T,   // Expiry time
			T/N  // Timestep size
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

		// Penalty method
		PenaltyMethod penalty(
			grid,           // Domain
			discretization, // Stochastic control component
			impulsePolicy   // Impulse control component
		);
		penalty.setIteration(toleranceIteration);

		IterationNode *root;
		if(method & EXPLICIT_IMPULSE) {
			// No impulse root
			root = &discretization;
		} else {
			// Impulse root
			root = &penalty;
		}

		////////////////////////////////////////////////////////////////
		// Exercise events
		////////////////////////////////////////////////////////////////

		auto ctSubroutine = [=,&stepper] (const Interpolant2 &u, Real W,
				Real B) {

			// Do nothing
			Real best = u(W, B);

			// Optimal controls
			Real sc = numeric_limits<Real>::infinity();
			Real ic = numeric_limits<Real>::infinity();

			const Real dt = stepper.timestep();

			// Consume
			if(method & SEMI_LAGRANGIAN_WITHDRAWAL_CONTINUOUS) {
				const Real dq = q_max
						/ (stochasticControlPoints - 1);

				for(
					int i = 1;
					i < stochasticControlPoints;
					++i
				) {
					const Real q = i * dq;

					// Don't have enough money to consume
					// that much!
					if(B < q * dt) { continue; }

					const Real newValue =
						u(W, B - q * dt)
						+ pow( q, aver ) / aver * dt
					;

					if(newValue > best) {
						best = newValue;
						sc = q;
					}
				}
			}

			// Transfer money from accounts
			const Real lo = -W;
			const Real hi = (B - c)/(1 + lambda);
			if(
				(method & EXPLICIT_IMPULSE)
				&& lo <= hi
			) {
				const Real dz = 1. / (impulseControlPoints - 1);

				for(
					int i = 1;
					i < impulseControlPoints;
					++i
				) {
					const Real zeta_frac = i * dz;

					const Real newValue = u(
						Wplus(0., W, B, zeta_frac),
						Bplus(0., W, B, zeta_frac)
					);

					if(newValue > best) {
						best = newValue;
						sc = numeric_limits<Real>
								::infinity();
						ic = zeta(W, B, zeta_frac);
					}
				}
			}

			return make_tuple(best, sc, ic);

		};

		auto ct = [=,&stepper] (const Interpolant2 &u, Real W, Real B) {
			return get<0>( ctSubroutine(u, W, B) );
		};

		auto ctPrintControl = [=,&stepper] (const Interpolant2 &u,
				Real W, Real B) {
			Real best, sc, ic;
			tie(best, sc, ic) = ctSubroutine(u, W, B);
			cout << W << "\t" << B << "\t" << sc << "\t" << ic
					<< endl;
			return best;
		};

		if(method != IMPLICIT) {
			// Special routine for printing
			const int mstart = (
				op == ProgramOperation::PLOT_CONTROLS
			) ? 1 : 0;
			if(mstart) { stepper.add(0., ctPrintControl, grid); }

			for(int m = mstart; m < N; ++m) {
				stepper.add(
					// Time at which the event takes place
					T / N * m,

					// Consume or transfer
					ct,

					// Spatial grid to interpolate on
					grid
				);
			}
		}

		////////////////////////////////////////////////////////////////
		// Running
		////////////////////////////////////////////////////////////////

		BiCGSTABSolver solver;
		auto u = stepper.solve(
			grid,                            // Domain
			[=] (Real, Real) { return 0.; }, // Payoff
			*root,                           // Root
			solver                           // Linear system solver
		);

		////////////////////////////////////////////////////////////////
		// Print solution at (w_0, b_0)
		////////////////////////////////////////////////////////////////

		if(op == ProgramOperation::PLOT) {

			RectilinearGrid2 printGrid(
				Axis::uniform( 0., 200., 200 ),
				Axis::uniform( 0., 200., 200 )
			);
			cout << accessor(printGrid, u);

		} else if(op == ProgramOperation::PLOT_CONTROLS) {

			// Print stochastic and impulse control

			auto mask = penalty.constraintMask();

			Vector ctrl_s = discretizee.control(0);
			for(int i = 0; i < grid.size(); i++) {
				if(mask[i]) {
					ctrl_s(i) = numeric_limits<Real>
							::infinity();
				}
			}

			Vector ctrl_i = impulseTransfer.control(0);
			for(int i = 0; i < grid.size(); i++) {
				if(!mask[i]) {
					ctrl_i(i) = numeric_limits<Real>
							::infinity();
				}
			}

			int k = 0;
			for(auto node : grid) {
				cout
					<< node[0]   << "\t"
					<< node[1]   << "\t"
					<< ctrl_s(k) << "\t"
					<< ctrl_i(k) << "\t"
					<< endl
				;
				++k;
			}

		} else {

			Real value = u( w_0, b_0 );
			Real
				change = value - previousValue,
				ratio = previousChange / change
			;
			previousValue = value;
			previousChange = change;

			cout.precision(12);

			// Print row of table
			cout
				<< setw(td) << grid.size()               << "\t"
				<< setw(td) << stochasticControls.size() << "\t"
				<< setw(td) << impulseControls.size()    << "\t"
				<< setw(td) << N                         << "\t"
				<< setw(td) << value                     << "\t"
			;

			if(method != EXPLICIT) {
				auto its = toleranceIteration.iterations();
				const Real mean = accumulate(its.begin(),
						its.end(), 0.) / its.size();
				cout << setw(td) << mean << "\t";
			}

			cout
				<< setw(td) << change << "\t"
				<< setw(td) << ratio  << "\t"
				<< endl
			;

		}

	}

	return 0;
}
