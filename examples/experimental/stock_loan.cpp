////////////////////////////////////////////////////////////////////////////////
// stock_loan.cpp
// --------------
//
// Author: Parsiad Azimzadeh
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>
#include <QuantPDE/Modules/Lambdas>
#include <QuantPDE/Modules/Operators>
#include <QuantPDE/Modules/Utilities>

using namespace QuantPDE;
using namespace QuantPDE::Modules;

////////////////////////////////////////////////////////////////////////////////

#include <algorithm>     // max, min
#include <cassert>       // assert
#include <cmath>         // exp
#include <limits>        // numeric_limits
#include <numeric>       // accumulate
#include <tuple>         // tie
#include <unordered_map> // unordered_map

using namespace std;

////////////////////////////////////////////////////////////////////////////////

Real
	T, r, divs, vol, S_0,
	q_hat, xi, rho,
	lockout_time,
	liquidation_trigger, margin_call_trigger, post_margin_call_ratio,
	jump_arrival, up_probability, up_mean_r, down_mean_r
;
bool
	explicit_impulses,
	margin_call_enabled, liquidation_enabled,
	event_at_zero
;
int N;
RectilinearGrid1 *grid;

////////////////////////////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////////////////////////////

constexpr Real inf = numeric_limits<Real>::infinity();

/**
 * Interpolation function for similarity reduction.
 * @param v The value of the stock loan for fixed q_hat.
 * @param s Collateral value.
 * @param q The outstanding loan.
 * @param q_hat A fixed, positive number.
 * @return The interpolated value of v at (s, q).
 */
inline Real v_similarity_reduction(
	const Interpolant1 &v,
	Real s, Real q, Real q_hat
) {
	const Real alpha = q_hat / q;
	return v(alpha * s) / alpha;
}

unordered_map<Real, Real> growth_cache;

/**
 * Computes the loan growth factor for some time.
 * @param t A time.
 * @return The growth factor.
 */
inline Real growth(Real t) {
	if(growth_cache.count(t) > 0) {
		return growth_cache[t];
	}
	growth_cache[t] = exp((r+xi) * t);
}

/**
 * Computes the loan-to-value ratio.
 * @param t A time.
 * @param s Collateral value.
 * @return The ratio.
 */
inline Real loan_to_value_ratio(Real t, Real s) {
	const Real g = growth(t);
	return g * q_hat / s;
}

/**
 * Computes the payoff for prepayment.
 * @param t A time.
 * @param s Collateral value.
 * @return The payoff or -infinity if t < lockout time.
 */
inline Real client_obstacle(Real t, Real s) {
	if(t < lockout_time) { return -inf; }
	return max(s - rho * growth(t) * q_hat, 0. );
}

/**
 * Computes cash flow to the client when the lender performs an action.
 * @param t A time.
 * @param s Collateral value.
 * @param z The control. Zero indicates liquidation, positive indicates a margin
 *          call.
 */
inline Real flow_to_client(Real t, Real s, Real z) {
	const Real g = growth(t);

	// Loan-to-value ratio
	const Real ratio = g * q_hat / s;

	if(
		liquidation_enabled
		&& z <= epsilon && ratio >= liquidation_trigger
	) {
		// Liquidate
		return max(s - g * q_hat, 0.);
	} else if(
		margin_call_enabled
		&& z > epsilon && ratio >= margin_call_trigger && s > 0.
	) {
		// Margin call
		return -(g * q_hat - z * s);
	}

	// Cannot issue a margin call or liquidate
	return inf;
}

/**
 * Computes the post margin call loan value.
 * @param t A time.
 * @param s Collateral value.
 * @param z The control.
 */
inline Real post_margin_call_loan(Real t, Real s, Real z) {
	const Real g = growth(t);
	return z * s / g;
}

/**
 * Computes the payoff at expiry.
 * @param s Collateral value.
 * @return Payoff.
 */
Real payoff(Real s) {
	return max(s - growth(T) * q_hat, 0.);
}

////////////////////////////////////////////////////////////////////////////////
// Lender intervention
////////////////////////////////////////////////////////////////////////////////

class LenderIntervention final : public RawControlledLinearSystem<1, 1> {

	const RectilinearGrid1 &grid;

	virtual Matrix A(Real t) {
		// S axis
		const Axis &S = grid[0];

		// Reserve matrix (at most 2 elements per row)
		Matrix M = grid.matrix();
		M.reserve(IntegerVector::Constant( grid.size(), 2 ));

		Index k = 0;
		for(auto node : grid) {
			// s-Coordinate
			const Real s = node[0];

			// Get control
			const Real z = (this->control(0))(k);

			// Loan-to-value ratio
			const Real ratio = loan_to_value_ratio(t, s);

			// Either liquidating or not past margin call trigger
			if(
				!(
					margin_call_enabled
					&& z > epsilon
					&& ratio >= margin_call_trigger
				)
			) {
				++k;
				continue;
			}

			// Loan value after margin call
			const Real q_plus = post_margin_call_loan(t, s, z);

			// Similarity reduction factor
			const Real alpha = q_hat / q_plus;

			// Approximate
			//   V(alpha * s) / alpha
			// via linear interpolation

			Index j;
			Real w;
			tie(j, w) = linearInterpolationData(S, alpha * s);

			M.insert(k, j  ) =    w  / alpha;
			M.insert(k, j+1) = (1-w) / alpha;

			++k;
		}

		M.makeCompressed();

		return grid.identity() - M;
	}

	virtual Vector b(Real t) {
		Vector b(grid.vector());

		Index k = 0;
		for(auto node : grid) {
			// s-Coordinate
			const Real s = node[0];

			// Get control
			const Real z = (this->control(0))(k);

			// Flow to client due to lender action
			b(k) = flow_to_client(t, s, z);

			++k;
		}

		return b;
	}

public:

	template <typename G>
	LenderIntervention(G &grid) : grid(grid) {}

};

////////////////////////////////////////////////////////////////////////////////

ResultsTuple1 run(int k) {

	// 2^k
	int factor = 1;
	for(int i = 0; i < k; ++i) {
		factor *= 4;
	}
	const int Nck = N * factor;

	////////////////////////////////////////////////////////////////////////
	// Spatial grid
	////////////////////////////////////////////////////////////////////////

	auto refined_grid = grid->refined(k); // Refine grid R times

	////////////////////////////////////////////////////////////////////////
	// Iteration tree
	////////////////////////////////////////////////////////////////////////

	const Real dt = T / Nck;
	ReverseConstantStepper stepper(0., T, dt);
	ToleranceIteration tolerance;

	// Explicit events
	const int n0 = event_at_zero ? 0 : 1;
	for(int n = n0; n < Nck; ++n) { // Skip n = Nck as it has no effect
		const Real t = n * dt;

		// Client
		stepper.add(
			t,
			[=] (const Interpolant1 &v, Real s) {
				// Do nothing
				const Real v0 = v(s);

				// Prepay
				const Real v1 = client_obstacle(t, s);

				// Maximum
				return max(v0, v1);
			},
			refined_grid
		);

		// Lender
		if(
			explicit_impulses
			&& (liquidation_enabled || margin_call_enabled)
		) {
			stepper.add(
				t,
				[=] (const Interpolant1 &v, Real s) {
					// Do nothing
					const Real v0 = v(s);

					// Liquidate
					const Real v1 = flow_to_client(
						t, s,
						0.
					);

					// Margin call
					const Real v2 = v_similarity_reduction(
						v, s,
						post_margin_call_loan(
							t,
							s,
							post_margin_call_ratio
						),
						q_hat
					) + flow_to_client(
						t, s,
						post_margin_call_ratio
					);

					// Minimum
					return min(v0, min(v1, v2));
				},
				refined_grid
			);
		}
	}

	////////////////////////////////////////////////////////////////////////
	// Linear system tree
	////////////////////////////////////////////////////////////////////////

	IterationNode *root;
	LinearSystem *bsptr;

	// Black-Scholes operator
	BlackScholesJumpDiffusion1 bsj(
		refined_grid, r, vol, divs,
		jump_arrival,
		doubleExponential(
			up_probability,
			up_mean_r,
			down_mean_r
		)
	);
	BlackScholes1 bs(refined_grid, r, vol, divs);
	if(jump_arrival < epsilon) {
		bsptr = &bs;
	} else {
		bsptr = &bsj;
		bsj.setIteration(stepper);
	}

	// Discretization
	typedef ReverseBDFOne Discretization;
	Discretization discretization(refined_grid, *bsptr);
	discretization.setIteration(stepper);

	// Policy iteration on lender control
	LenderIntervention lender_intervention(refined_grid);
	auto control_grid = RectilinearGrid1(
		Axis { 0., post_margin_call_ratio }
	);
	MaxPolicyIteration1_1 policy(
		refined_grid,
		control_grid,
		lender_intervention
	);
	policy.setIteration(tolerance);
	MaxPenaltyMethod penalty(
		refined_grid,
		discretization,
		policy
	);
	penalty.setIteration(tolerance);

	// Pick root element
	if(explicit_impulses) {
		root = &discretization;
	} else {
		stepper.setInnerIteration(tolerance);
		root = &penalty;
	}

	////////////////////////////////////////////////////////////////////////
	// Running
	////////////////////////////////////////////////////////////////////////

	SparseLUSolver solver;

	auto solution = stepper.solve(
		refined_grid, // Domain
		payoff,       // Initial condition
		*root,        // Root of linear system tree
		solver        // Linear system solver
	);

	////////////////////////////////////////////////////////////////////////

	// Timesteps
	unsigned timesteps = stepper.iterations()[0];

	// Average number of policy iterations
	Real policy_its = nan("");
	auto its = tolerance.iterations();
	policy_its = accumulate(its.begin(), its.end(), 0.) / its.size();

	return ResultsTuple1(
		{(Real) refined_grid.size(), (Real) timesteps, policy_its },
		solution, S_0
	);
}

int main(int argc, char **argv) {
	// Parse configuration file
	Configuration configuration = getConfiguration(argc, argv);

	// Get options
	int kn, k0;
	Real S_max, S_min, dS;
	kn = getInt(configuration, "maximum_refinement", 7);
	k0 = getInt(configuration, "minimum_refinement", 0);
	T = getReal(configuration, "time_to_expiry", 10.);
	r = getReal(configuration, "interest_rate", .05);
	divs = getReal(configuration, "dividend_rate", 0.02);
	vol = getReal(configuration, "volatility", .15);
	S_0 = getReal(configuration, "asset_price", 100.);
	q_hat = getReal(configuration, "loan_value", 80.);
	xi = getReal(configuration, "spread", 0.02);
	rho = getReal(configuration, "penalty_scaling", 1.);
	lockout_time = getReal(configuration, "lockout_time", 0.);
	liquidation_trigger = getReal(configuration, "liquidation_trigger", 80. / 90.);
	margin_call_trigger = getReal(configuration, "margin_call_trigger", 80. / 95.);
	jump_arrival = getReal(configuration, "jump_arrival_rate", 0.5);
	up_probability = getReal(configuration, "jump_up_probability", 0.09);
	up_mean_r = getReal(configuration, "jump_up_mean_reciprocal", 2.3);
	down_mean_r = getReal(configuration, "jump_down_mean_reciprocal", 1.8);
	post_margin_call_ratio = getReal(configuration, "post_margin_call_ratio", 80. / 100.);
	S_min = getReal(configuration, "print_asset_price_minimum", 0.);
	S_max = getReal(configuration, "print_asset_price_maximum", 300.);
	dS = getReal(configuration, "print_asset_price_step_size", 10.);
	N = getInt(configuration, "initial_number_of_timesteps", 2);
	explicit_impulses = getInt(configuration, "fully_explicit", false);
	event_at_zero = getInt(configuration, "apply_event_at_time_zero", true);
	margin_call_enabled = getInt(configuration, "margin_call_enabled", false);
	liquidation_enabled = getInt(configuration, "liquidation_enabled", true);
	RectilinearGrid1 default_grid( S_0 * Axis::special );
	RectilinearGrid1 tmp = getGrid(configuration, "initial_grid", default_grid);
	grid = &tmp;

	// Print configuration file
	cerr << configuration << endl << endl;

	// Run and print results
	ResultsBuffer1 buffer(
		run,
		{ "Nodes", "Steps", "Mean Policy Iterations" },
		kn, k0
	);
	buffer.addPrintGrid( RectilinearGrid1(Axis::range(S_min, dS, S_max)) );
	buffer.stream();

	return 0;
}

