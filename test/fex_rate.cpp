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

#include <QuantPDE/Core>

#include <cstdlib>  // abs
#include <iomanip>  // setw
#include <iostream> // cout, cerr

////////////////////////////////////////////////////////////////////////////////

using namespace QuantPDE;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

// The problem is solved on [-boundary, boundary]
Real boundary = 7.; // log(1000) ~= 7

// The stochastic control is truncated to [-betaMax, betaMax]
Real betaMin = 0.;
Real betaMax = 0.2;

// Number of points in space
int N = 128;

// Number of points for control
int M = 4;

// Cost of intervention is lambda * |zeta| + c (zeta is the intervention size)
Real lambda = 1.;
Real c = 0.;

// Effect of an intervention
Real a = 1.;

// Discount
Real rho = 0.04;

// Volatility
Real v = 0.2;

// Foreign interest rate
Real rbar = 0.04;

// Central parity
Real m = 0.;

int maxRefinement = 10;

constexpr int td = 20;

////////////////////////////////////////////////////////////////////////////////

/**
* Cost of an impulse.
*/
inline Real interventionCost(Real t, Real x, Real x_j) {
	const Real zeta = (x_j - x) / a;
	return -(lambda * abs(zeta) + c);
}

/**
* State of the exchange rate after an impulse.
*/
inline Real xplus(Real t, Real x, Real x_j) {
	return x_j;
}

////////////////////////////////////////////////////////////////////////////////

class Discretizee final : public RawControlledLinearSystem1_1 {

	const RectilinearGrid1 &grid;

	const Function1 F, N, M;
	const Real rho, v, rbar, m;

	// TODO: Nonconstant coefficients

public:

	template <typename G, typename F1, typename F2, typename F3>
	Discretizee(
		G &grid,
		F1 &&interestDifferentialEffect,
		F2 &&interestDifferentialCost,
		F3 &&exchangeRateCost,
		Real discount,
		Real volatility,
		Real foreignInterestRate = 0.,
		Real centralParity = 0.
	) noexcept :
		grid(grid),
		F( std::forward<F1>(interestDifferentialEffect) ),
		N( std::forward<F2>(interestDifferentialCost) ),
		M( std::forward<F3>(exchangeRateCost) ),
		rho( discount ),
		v( volatility ),
		rbar( foreignInterestRate ),
		m( centralParity )
	{
	}

	virtual Matrix A(Real time) {
		Matrix A = grid.matrix();
		A.reserve( IntegerVector::Constant(grid.size(), 3) );

		// x axis
		const Axis &x = grid[0];

		// Control as a vector
		const Vector &raw = control(0);

		// Left boundary (derivatives disappear)
		A.insert(0, 0) = rho;

		// Interior point
		for(Index i = 1; i < x.size() - 1; ++i) {
			const Real
				dxb = x[i]     - x[i - 1],
				dxc = x[i + 1] - x[i - 1],
				dxf = x[i + 1] - x[i]
			;

			const Real tmp1 = v * v;
			const Real tmp2 = -F(raw(i) - rbar);

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

			A.insert(i, i - 1) = -alpha_i;
			A.insert(i, i    ) =  alpha_i + beta_i + rho;
			A.insert(i, i + 1) = -beta_i;
		}

		// Right boundary (derivatives disappear)
		A.insert(x.size() - 1, x.size() - 1) = rho;

		A.makeCompressed();
		return A;
	}

	virtual Vector b(Real time) {
		Vector b = grid.vector();

		const Axis &x = grid[0];

		// Control as a vector
		const Vector &raw = control(0);

		// Interior and boundary points
		for(Index i = 0; i < x.size(); ++i) {
			b(i) = -( M(x[i] - m) + N(raw(i) - rbar) );
		}

		return b;
	}

	virtual bool isATheSame() const {
		// TODO: Check the coefficients
		return true;
	}

};

////////////////////////////////////////////////////////////////////////////////

int main() {

	RectilinearGrid1 grid(
		// Uniform grid
		Axis::uniform(
			-boundary,
			+boundary,
			N
		)
	);

	RectilinearGrid1 stochasticControls(
		Axis::uniform(
			betaMin,
			betaMax,
			M
		)
	);

	for(int ref = 0; ref < maxRefinement; ++ref) {

		ToleranceIteration toleranceIteration;

		////////////////////////////////////////////////////////////////
		// Linear system tree
		////////////////////////////////////////////////////////////////

		Discretizee discretizee(
			grid,
			[] (Real y) { return y  ; }, // F(y)
			[] (Real y) { return y  ; }, // N(y)
			[] (Real y) { return y*y; }, // M(y)
			rho,
			v,
			rbar,
			m
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

		// Penalty method
		PenaltyMethod penalty(
			grid,             // Domain
			stochasticPolicy, // Stochastic control component
			impulsePolicy     // Impulse control component
		);
		penalty.setIteration(toleranceIteration);

		////////////////////////////////////////////////////////////////
		// Run
		////////////////////////////////////////////////////////////////

		BiCGSTABSolver solver;
		auto u = toleranceIteration.solve(
			grid,                        // Domain
			[=] (Real x) { return 0.; }, // Initial guess (zero
			                             // everywhere)
			penalty,                     // Root of linear system
			                             // tree
			solver                       // Linear system solver
		);

		////////////////////////////////////////////////////////////////
		// Print solution
		////////////////////////////////////////////////////////////////

		cout
			<< setw(td) << -u(0.) << "\t"
			<< setw(td) << toleranceIteration.iterations().back()
			<< endl
		;

		// Refine grid and control set
		grid = grid.refined();
		stochasticControls = stochasticControls.refined();

	}

	return 0;

}
