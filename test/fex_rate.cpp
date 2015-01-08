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

#include <cstdlib> // abs

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
int n = 32;

// Number of points for control
int m = 4;

// Cost of intervention is lambda * |zeta| + c (zeta is the intervention size)
Real lambda = 1.;
Real c = 0.;

// Effect of an intervention
Real a = 1.;

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

	const Real rho, v, rbar, m;

public:

	template <typename G>
	Discretizee(
		G &grid,
		Real discount,
		Real volatility,
		Real foreignInterestRate = 0.,
		Real centralParity = 0.
	) noexcept :
		grid(grid),
		rho(discount),
		v(volatility),
		rbar(foreignInterestRate),
		m(centralParity)
	{
	}

	virtual Matrix A(Real time) {
		Matrix M = grid.matrix();
		M.reserve( IntegerVector::Constant(grid.size(), 3) );

		// x axis
		const Axis &x = grid[0];

		// Control as a vector
		const Vector &raw = controlled ? control(0) : zero;

		// Left boundary
		// TODO

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

			// TODO
			M.insert(idx, idx - offset) = -alpha_i;
			M.insert(idx, idx)          =  alpha_i + beta_i + rho;
			M.insert(idx, idx + offset) = -beta_i;
		}

		// Right boundary
		// TODO

		M.makeCompressed();
		return M;
	}

	virtual Vector b(Real time) {
		Vector b = grid.vector();

		const Axis &x = grid[0];

		// Control as a vector
		const Vector &raw = control(0);

		// Left boundary
		// TODO

		// Interior point
		for(Index i = 1; i < x.size() - 1; ++i) {
			b(i) = -( M(x[i] - m) + N(raw(i) - rbar) );
			++k;
		}

		// Right boundary
		// TODO

		return b;
	}

};

////////////////////////////////////////////////////////////////////////////////

int main() {

	RectilinearGrid1 grid(
		// Uniform grid
		/*
		Axis::uniform(
			-boundary,
			+boundary,
			n
		)
		*/

		// Automagic grid
		Axis::cluster(
			-boundary,
			+boundary,
			n,
			0., // Feature
			5.  // Intensity
		)
	);

	RectilinearGrid1 stochasticControls(
		Axis::uniform(
			betaMin,
			betaMax,
			m
		)
	);

	ToleranceIteration toleranceIteration;

	////////////////////////////////////////////////////////////////////////
	// Linear system tree
	////////////////////////////////////////////////////////////////////////

	Discretizee discretizee;

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
		xplus             // State of exchange rate after intervention
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

	////////////////////////////////////////////////////////////////////////
	// Running
	////////////////////////////////////////////////////////////////////////

	BiCGSTABSolver solver;
	auto psi = toleranceIteration.solve(
		grid,                        // Domain
		[=] (Real X) { return 0.; }, // Initial guess (zero everywhere)
		penalty,                     // Root of linear system tree
		solver                       // Linear system solver
	);

	////////////////////////////////////////////////////////////////////////
	// Print solution
	////////////////////////////////////////////////////////////////////////

	RectilinearGrid1 printGrid( Axis::uniform(-boundary, +boundary, 10) );
	cout << accessor( printGrid, psi );

	return 0;

}
