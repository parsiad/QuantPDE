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
#include <fstream>  // ofstream
#include <iomanip>  // setw
#include <iostream> // cout, cerr
#include <string>   // to_string

////////////////////////////////////////////////////////////////////////////////

using namespace QuantPDE;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

// The problem is solved on [-boundary, boundary]
Real boundary = 7.; // log(1000) ~= 7

// Interest rate differential in [-betaMax, +betaMax]
Real betaMax = 0.5;

// Number of points in space (initial discretization)
int gridPoints = 32;

// Number of points for interest rate control (initial discretization)
int controlPoints = 2;

// Constants
Real c_1 = 1.;
Real c_2 = 1.;
Real c_3 = QuantPDE::epsilon;
Real c_4 = 1.;
Real c_5 = 1.;
Real c_6 = 1.;

// Effect of an intervention
Real a = 1.;

// Discount
Real rho = 0.04;

// Volatility
Real v = 0.2;

// Central parity
Real m = 0.;

// Initial (log of the) exchange rate
Real x_0 = 0.;

int maxRefinement = 10;

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
	return -(c_2 * abs(zeta) + c_3);
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
* Cost of exchange ate differential.
*/
inline Real M(Real y) {
	return c_5 * y * y;
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
		Real centralParity = 0.
	) noexcept :
		grid(grid),
		F( std::forward<F1>(interestDifferentialEffect) ),
		M( std::forward<F3>(exchangeRateCost) ),
		N( std::forward<F2>(interestDifferentialCost) ),
		rho( discount ),
		v( volatility ),
		m( centralParity )
	{
	}

	virtual Matrix A(Real time) {
		Matrix A = grid.matrix();
		A.reserve( IntegerVector::Constant(grid.size(), 3) );

		// x axis
		const Axis &x = grid[0];
		const int n = x.size();

		// Control as a vector
		const Vector &raw = control(0);

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

int main() {

	Real previousValue = nan(""), previousChange = nan("");

	RectilinearGrid1 grid(
		// Uniform grid
		Axis::uniform(
			-boundary,
			m,
			gridPoints
		)

		// Clustered grid
		/*Axis::cluster(
			-boundary,  // Left-hand boundary
			0.,         // Feature to cluster around
			0.,         // Right-hand boundary
			gridPoints, // Number of points
			10.         // Clustering intensity
		)*/
	);

	// Table headers
	cout
		<< setw(td) << "Nodes"         << "\t"
		<< setw(td) << "Control Nodes" << "\t"
		<< setw(td) << "Value at x=0"  << "\t"
		<< setw(td) << "Iterations"    << "\t"
		<< setw(td) << "Change"        << "\t"
		<< setw(td) << "Ratio"
		<< endl
	;

	for(int ref = 0; ref < maxRefinement; ++ref, controlPoints *= 2) {

		RectilinearGrid1 stochasticControls(
			Axis::uniform(
				-betaMax,
				+betaMax,
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
			m     // Central parity
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

		IterationNode *root;
		//root = &stochasticPolicy;
		root = &penalty;

		////////////////////////////////////////////////////////////////
		// Run
		////////////////////////////////////////////////////////////////

		BiCGSTABSolver solver;
		auto u = toleranceIteration.solve(
			grid,                        // Domain
			[=] (Real x) { return 0.; }, // Initial guess (zero
			                             // everywhere)
			*root,                       // Root of linear system
			                             // tree
			solver                       // Linear system solver
		);

		////////////////////////////////////////////////////////////////
		// Print solution at x_0
		////////////////////////////////////////////////////////////////

		// If x_0 > m, use the fact that u(.) is symmetric about x=m
		Real adjusted;
		if(x_0 > m) {
			adjusted = m - (x_0 - m);
		} else {
			adjusted = x_0;
		}

		Real value = -u( adjusted );
		Real
			change = value - previousValue,
			ratio = previousChange / change
		;
		previousValue = value;
		previousChange = change;

		cout.precision(12);

		const int its = toleranceIteration.iterations().back();

		// Print row of table
		cout
			<< setw(td) << grid.size()               << "\t"
			<< setw(td) << stochasticControls.size() << "\t"
			<< setw(td) << value                     << "\t"
			<< setw(td) << its                       << "\t"
			<< setw(td) << change                    << "\t"
			<< setw(td) << ratio                     << "\t"
			<< endl
		;

		////////////////////////////////////////////////////////////////
		// Write to file
		////////////////////////////////////////////////////////////////

		// Write to file
		/*
		RectilinearGrid1 printGrid(
			Axis::uniform(
				-boundary,
				0.,
				100
			)
		);
		ofstream file;
		file.open("/tmp/fex_rate_ref_" + to_string(ref) + ".txt");
		file << accessor( printGrid, u );
		file.close();
		*/

		////////////////////////////////////////////////////////////////
		// Refinements
		////////////////////////////////////////////////////////////////

		grid = grid.refined();

	}

	return 0;

}
