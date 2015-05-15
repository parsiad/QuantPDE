////////////////////////////////////////////////////////////////////////////////
// optimal_consumption.cpp
// -----------------------
//
// TODO: Description
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>

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

// Stochastic control
Real q_max = 1.;
int controlPoints = 16;

// Interest rate
Real r = 0.02;

// Volatility
Real v = 0.3;

// Drift
Real mu = 0.1;

// Risk-aversion
Real gamma = 1.;

// Discount factor
Real rho = 0.02;

// Finite horizon
Real T = 10.;

// Initial number of timesteps
int N = 32;

// Max/min refinement level
int Rmin = 0;
int Rmax = 6;

////////////////////////////////////////////////////////////////////////////////
// Solution grid
////////////////////////////////////////////////////////////////////////////////

Axis axis = (w0 / 100.) * Axis {
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

// Grid
RectilinearGrid2 initialGrid( axis, axis );

////////////////////////////////////////////////////////////////////////////////
// Operator to discretize
// i.e. Lu - rho + f, where L is the infinitesimal generator of process
//      (B_t, W_t)
////////////////////////////////////////////////////////////////////////////////

class Discretizee final : public RawControlledLinearSystem2_1 {

	const RectilinearGrid2 &grid;

	const Real r, v, mu, rho, gamma;

	const bool controlled;
	const Vector zero;

	template <typename G>
	Discretizee(
		G1 &grid,
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
		gamma( aversion ),
		controlled( controlled ),
		zero( grid.zero() )
	{
	}

	virtual Matrix A(Real) {
		Matrix M = grid.matrix();
		M.reserve( IntegerVector::Constant(grid.size(), 4) );

		// Axes
		const Axis &W = grid[0];
		const Axis &B = grid[1];

		// Control as vector
		Index k = 0;
		const Vector &raw = controlled ? control(0) : zero;

		for(Index j = 0; j < B.size(); ++j) {
			for(Index i = 0; i < W.size(); ++i) {

				const Real
					dWb = W[i  ] - W[i-1],
					dWc = W[i+1] - W[i-1],
					dWf = W[i+1] - W[i  ],
					dBb = B[j  ] - B[j-1],
					dBc = B[j+1] - B[j-1],
					dBf = B[j+1] - B[j  ]
				;

				const Real tmp1 = v * v * W[i] * W[i];
				const Real tmp2 = mu * W[i];

				const Real alpha_common = tmp1 / dWb / dWc;
				const Real  beta_common = tmp1 / dWf / dWc;

				// Central
				Real alpha_i = alpha_common - tmp2 / dWc;
				Real  beta_i =  beta_common + tmp2 / dWc;

				if(alpha_i < 0.) {
					// Forward
					alpha_i = alpha_common;
					 beta_i =  beta_common + tmp2 / dWf;
				} else if(beta_i < 0.) {
					// Backward
					alpha_i = alpha_common - tmp2 / dWb;
					 beta_i =  beta_common;
				}

				if(j == 0) {
					if(i < W.size()) {
						// (B=0, 0<=W<W_max)
					} else {
						// (B=0, W=W_max)
					}
				} else if(j < B.size()) {
					if(i < W.size()) {
						// (0<B<B_max, 0<=W<W_max)
					} else {
						// (0<B<B_max, W=W_max)
					}
				} else {
					if(i < W.size()) {
						// (B=B_max, 0<=W<W_max)
					} else {
						// (B=B_max, W=W_max)
					}
				}

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

		// B > 0
		for(Index j = 1; j < B.size(); ++j) {
			// 0 <= W <= W_max
			for(Index i = 0; i < W.size(); ++i) {
				b(k) = raw(k);
				++k;
			}
		}

		return b;
	}

	virtual bool isATheSame() const {
		return !controlled;
	}

};

int main() {

	for(int ref = 0; ref <= Rmax; ++ref, N *= 2, controlPoints *= 2) {

		if(ref < Rmin) { continue; }

		// Refine the grid ref times
		auto grid = initialGrid.refined( ref );

		RectilinearGrid1 stochasticControls(
			Axis::uniform(
				0.,           // Control lower bound
				q_max,        // Control upper bound
				controlPoints // Number of nodes
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
			gamma, // Risk-aversion

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
		Impulse2_1 impulseWithdrawal(
			// Spatial grid
			grid,

			// The intervention does not cost anything; the cost
			// is taken out of the bank account
			[] (Real t, Real w, Real b, Real w_j, Real b_j) {
				return 0.;
			},

			// State after an intervention
			xplus
		);
		MinPolicyIteration2_1 impulsePolicy(
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

	}

	return 0;
}
