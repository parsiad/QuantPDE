////////////////////////////////////////////////////////////////////////////////
// stock_loan.cpp
// --------------
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
#include <cassert>   // assert
#include <cmath>     // exp
#include <limits>    // numeric_limits

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/**
 * Interpolation function for similarity reduction.
 * @param v The value of the stock loan for fixed q_hat.
 * @param s The value of the collateral.
 * @param q The outstanding loan.
 * @param q_hat A fixed, positive number.
 */
inline Real similarity_value(
	const Interpolant1 &v,
	Real s, Real q,
	Real q_hat
) {
	assert(q > 0.);
	const Real alpha = q_hat / q;
	return v(alpha * s) / alpha;
}

int main(int argc, char **argv) {

	////////////////////////////////////////////////////////////////////////
	// Constants
	////////////////////////////////////////////////////////////////////////

	const Real q_hat = 80.;      // Special loan value used in computation
	const Real T = 5.;           // Expiry time
	const Real r = .02;          // Interest rate
	const Real v = .3;           // Volatility
	const Real xi = 0.01;        // Spread
	const Real rho = 1.;         // Penalty scaling (>= 0)
	const Real gamma = r + xi;   // Loan interest rate
	const Real P = 1.;           // Lockout time
	const Real beta = 80. / 89.; // Liquidation trigger
	const Real eta = 80. / 94.;  // Margin call trigger
	const Real lambda = 1;       // Jump arrival rate
	const Real mu_xi = -.8;      // Mean jump amplitude
	const Real sigma_xi = .42;   // Jump amplitude standard deviation
	const Real theta = .8;       // Post-margin call loan-to-value ratio

	const int N = 10000;         // Timesteps
	const int R = 3;             // Refinement

	////////////////////////////////////////////////////////////////////////
	// Spatial grid
	////////////////////////////////////////////////////////////////////////

	RectilinearGrid1 initialGrid( q_hat * Axis::special );
	auto grid = initialGrid.refined(R); // Refine grid R times

	////////////////////////////////////////////////////////////////////////
	// Payoffs and impulses
	////////////////////////////////////////////////////////////////////////

	auto repay = [=] (Real s) {
		const Real grow = exp(gamma * T);
		return max(s - grow * q_hat, 0.);
	};

	auto prepay = [=] (Real t, Real s) {
		const Real grow = exp(gamma * t);
		if(t >= P) { return max(s - rho * grow * q_hat, 0.); }
		else { return -numeric_limits<Real>::infinity(); }
	};

	auto liquidation = [=] (Real t, Real s) {
		const Real grow = exp(gamma * t);
		if(grow * q_hat >= beta * s) {
			return max(s - grow * q_hat, 0.);
		} else { return numeric_limits<Real>::infinity(); }
	};

	/*
	auto margin_call = [=] (const Interpolant1 &v, Real t, Real s, Real z) {
		const Real grow = exp(gamma * t);
		if(grow * q_hat >= eta * s) {
			return similarity_value(v, s, z * s / grow, q_hat)
					- (grow * q_hat - z * s);
		} else { return numeric_limits<Real>::infinity(); }
	};
	*/

	////////////////////////////////////////////////////////////////////////
	// Iteration tree
	////////////////////////////////////////////////////////////////////////

	const Real dt = T / N;
	ReverseConstantStepper stepper(0., T, dt);
	ToleranceIteration tolerance;
	stepper.setInnerIteration(tolerance);

	// Lender's events (handled explicitly)
	for(int n = 1; n < N; ++n) {
		const Real t = n * dt;
		stepper.add(
			n * dt,
			[=] (const Interpolant1 &v, Real s) {
				const Real v0 = v(s);
				const Real v1 = liquidation(t, s);
				//const Real v2 = margin_call(v, t, s, theta);
				//return min(v0, min(v1, v2));
				return min(v0, v1);
			},
			grid
		);
	}

	////////////////////////////////////////////////////////////////////////
	// Linear system tree
	////////////////////////////////////////////////////////////////////////

	BlackScholesJumpDiffusion1 bs(
		grid,

		r,  // Interest
		v,  // Volatility
		0., // Dividend rate

		lambda, // Mean arrival time
		lognormal(mu_xi, sigma_xi) // Log-normal probability density
	);
	bs.setIteration(stepper);

	typedef ReverseBDFOne Discretization;
	Discretization discretization(grid, bs);
	discretization.setIteration(stepper);

	PenaltyMethodDifference1 penalty(grid, discretization, prepay);
	penalty.setIteration(tolerance);

	////////////////////////////////////////////////////////////////////////
	// Running
	////////////////////////////////////////////////////////////////////////

	BiCGSTABSolver solver;

	auto V = stepper.solve(
		grid,    // Domain
		repay,   // Initial condition
		penalty, // Root of linear system tree
		solver   // Linear system solver
	);

	////////////////////////////////////////////////////////////////////////
	// Print solution
	////////////////////////////////////////////////////////////////////////

	RectilinearGrid1 printGrid( Axis::range(0., .1, 300.) );
	cout << accessor( printGrid, V );

	return 0;

}

