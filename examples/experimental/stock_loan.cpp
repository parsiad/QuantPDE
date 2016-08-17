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
#include <cmath>         // ceil, exp
#include <fstream>       // ofstream
#include <limits>        // numeric_limits
#include <memory>        // unique_ptr
#include <numeric>       // accumulate
#include <sstream>       // stringstream
#include <tuple>         // tie
#include <unordered_map> // unordered_map
#include <vector>        // vector

using namespace std;

////////////////////////////////////////////////////////////////////////////////

#define QUANT_PDE_STOCK_LOAN_EXPLICIT_ONLY

////////////////////////////////////////////////////////////////////////////////

Real
	T, r, divs, vol, S_0,
	q_hat, spread,
	prepayment_penalty, liquidation_penalty,
	lockout_time,
	liquidation_trigger, margin_call_trigger, post_margin_call_ratio,
	jump_arrival, up_probability, up_mean_r, down_mean_r
;
bool
	explicit_impulses,
	margin_call_enabled, liquidation_enabled,
	use_dirty_loan_to_value_ratio
;
int
	kn, k0,
	N,
	coupon_payments
;
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

/**
 * Computes the loan growth factor minus unity for some fixed increment in time.
 * @param delta A time increment.
 * @return The growth factor.
 */
inline Real growth_delta(Real delta) {
	return (r+spread) * delta;
}

/**
 * Computes the loan growth factor for some time.
 * @param t A time.
 * @return The growth factor.
 */
inline Real growth(Real t) {
	// Get last coupon payment time (or initial time)
	const Real d = (t / T) * (coupon_payments + 1);
	int n = floor(d);
	if(n > 0 && abs(n - d) < epsilon) { --n; }
	const Real t_n = n * (T / (coupon_payments + 1));
	return 1. + growth_delta(t - t_n);
}

/**
 * Computes the loan-to-value ratio.
 * @param t A time.
 * @param s Collateral value.
 * @return The ratio.
 */
inline Real loan_to_value_ratio(Real t, Real s) {
	if(use_dirty_loan_to_value_ratio) {
		return growth(t) * q_hat / s;
	}
	return q_hat / s;
}

/**
 * Computes the payoff for prepayment.
 * @param t A time.
 * @param s Collateral value.
 * @return The payoff or -infinity if t < lockout time.
 */
inline Real client_obstacle(Real t, Real s) {
	if(t < lockout_time) { return -inf; }
	return max(s - prepayment_penalty * growth(t) * q_hat, 0. );
}

/**
 * Computes cash flow to the client when the lender performs an action.
 * @param t A time.
 * @param s Collateral value.
 * @param control The control. Zero indicates liquidation, positive indicates a
 *                margin call.
 */
inline Real flow_to_client(Real t, Real s, Real control) {

	const Real g = growth(t);

	// Loan-to-value ratio
	const Real ratio = loan_to_value_ratio(t, s);

	if(
		liquidation_enabled
		&& control <= epsilon && ratio >= liquidation_trigger
	) {
		// Liquidate
		return max(s - liquidation_penalty * g * q_hat, 0.);
	} else if(
		margin_call_enabled
		&& control > epsilon && ratio >= margin_call_trigger
	) {
		// Margin call
		if(use_dirty_loan_to_value_ratio) {
			return -(g * q_hat - control * s);
		}
		return -g * (q_hat - control * s);
	}

	// Cannot issue a margin call or liquidate
	return inf;
}

/**
 * Computes the post margin call loan value.
 * @param t A time.
 * @param s Collateral value.
 * @param control The control.
 */
inline Real post_margin_call_loan(Real t, Real s, Real control) {
	if(use_dirty_loan_to_value_ratio) {
		return control * s / growth(t);
	}
	return control * s;
}

/**
 * Computes the payoff at expiry.
 * @param s Collateral value.
 * @return Payoff.
 */
Real payoff(Real s) {
	const Real coupon_dt = T / (coupon_payments + 1);
	return max(s - (1. + growth_delta(coupon_dt)) * q_hat, 0.);
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
			const Real control = (this->control(0))(k);

			// Loan-to-value ratio
			const Real ratio = loan_to_value_ratio(t, s);

			// Either liquidating or not past margin call trigger
			if(
				!(
					margin_call_enabled
					&& control > epsilon
					&& ratio >= margin_call_trigger
				)
			) {
				++k;
				continue;
			}

			// Loan value after margin call
			const Real q_plus = post_margin_call_loan(
				t,
				s,
				control
			);

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
			const Real control = (this->control(0))(k);

			// Flow to client due to lender action
			b(k) = flow_to_client(t, s, control);

			++k;
		}

		return b;
	}

public:

	template <typename G>
	LenderIntervention(G &grid) noexcept : grid(grid) {}

};

class ExplicitEvent : public EventBase {

	const Real t;
	const RectilinearGrid1 &grid;
	Vector &control;

	template <typename V>
	Vector _doEvent(V &&vector) const {
		const Axis &S = grid[0];

		Vector best = grid.vector();

		for(int i = 0; i < S.size(); ++i) {
			const Real S_i = S[i];

			// Interpolant
			PiecewiseLinear1 interpolant(grid, vector);

			// Do nothing
			const Real continue_i = vector(i);

			// Prepay
			const Real prepay_i = client_obstacle(t, S_i);

			// Liquidate
			const Real liquidate_i = flow_to_client(t, S_i, 0.);

			// Margin call
			const Real margin_call_i = v_similarity_reduction(
				interpolant, S_i,
				post_margin_call_loan(
					t, S_i,
					post_margin_call_ratio
				),
				q_hat
			) + flow_to_client(t, S_i, post_margin_call_ratio);

			// Saddle point
			best(i)  = continue_i;
			control(i) = 0;
			if(margin_call_i < best(i)) {
				best(i) = margin_call_i;
				control(i) = 1;
			}
			if(liquidate_i < best(i)) {
				best(i) = liquidate_i;
				control(i) = 2;
			}
			if(prepay_i > best(i)) {
				best(i) = prepay_i;
				control(i) = 3;
			}
		}

		return best;
	}

	virtual Vector doEvent(const Vector &vector) const {
		return _doEvent(vector);
	}

	virtual Vector doEvent(Vector &&vector) const {
		return _doEvent(std::move(vector));
	}

public:

	template <typename G>
	ExplicitEvent(Real t, G &grid, Vector &control) noexcept :
			t(t), grid(grid), control(control) {
		control.resize(grid.size()); // Size vector
	}

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

	Vector controls[Nck];

	// Explicit events
	for(int n = 0; n < Nck; ++n) { // Skip n = Nck as it has no effect
		const Real t = n * dt;

		#ifdef QUANT_PDE_STOCK_LOAN_EXPLICIT_ONLY
		stepper.add(
			t,
			unique_ptr<EventBase>(
				new ExplicitEvent(
					t,
					refined_grid,
					controls[n]
				)
			)
		);
		#else
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
		#endif
	}

	// Coupon payments
	const Real coupon_dt = T / (coupon_payments + 1);
	for(int n = 1; n <= coupon_payments; ++n) {
		const Real t = n * coupon_dt;

		stepper.add(
			t,
			[=] (const Interpolant1 &v, Real s) {
				return v(s) - q_hat * growth_delta(coupon_dt);
			},
			refined_grid
		);
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
		if(explicit_impulses) {
			bsj.setIteration(tolerance);
		} else {
			// Not sure why this works...
			bsj.setIteration(stepper);
		}
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
		policy,
		1e-2 * dt
	);
	penalty.setIteration(tolerance);

	// Penalty on client control
	/*
	MinPenaltyMethodDifference1 penalty(
		refined_grid,
		discretization,
		client_obstacle,
		1e-2 * dt
	);
	penalty.setIteration(tolerance);
	*/

	// Pick root element
	stepper.setInnerIteration(tolerance);
	if(explicit_impulses) {
		root = &discretization;
	} else {
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

	/*
	// TODO: Incorporate printing controls into ResultsTuple
	#ifdef QUANT_PDE_STOCK_LOAN_EXPLICIT_ONLY
	if(k == kn) {
		ofstream fout("control.txt");
		for(int n = 0; n < Nck; ++n) {
			const Real t = n * dt;
			const Axis &S = refined_grid[0];
			for(int i = 0; i < S.size(); ++i) {
				const Real S_i = S[i];
				Real ctrl = controls[n](i) == 3;
				fout << t << "\t" << S_i << "\t" << ctrl << endl;
			}
		}
	}
	#endif
	*/

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
	Real S_max, S_min, dS;
	kn = getInt(configuration, "maximum_refinement", 6);
	k0 = getInt(configuration, "minimum_refinement", 0);
	T = getReal(configuration, "time_to_expiry", 3.);
	r = getReal(configuration, "interest_rate", .05);
	divs = getReal(configuration, "dividend_rate", 0);
	vol = getReal(configuration, "volatility", .15);
	S_0 = getReal(configuration, "asset_price", 100.);
	q_hat = getReal(configuration, "loan_value", 70.);
	spread = getReal(configuration, "spread", 0.02);
	prepayment_penalty = getReal(configuration, "prepayment_penalty", 1.);
	liquidation_penalty = getReal(configuration, "liquidation_penalty", 1.);
	lockout_time = getReal(configuration, "lockout_time", 0.);
	liquidation_trigger = getReal(configuration, "liquidation_trigger", 70. / 75.);
	margin_call_trigger = getReal(configuration, "margin_call_trigger", 70. / 90.);
	jump_arrival = getReal(configuration, "jump_arrival_rate", 0.5);
	up_probability = getReal(configuration, "jump_up_probability", 0.09);
	up_mean_r = getReal(configuration, "jump_up_mean_reciprocal", 2.3);
	down_mean_r = getReal(configuration, "jump_down_mean_reciprocal", 1.8);
	post_margin_call_ratio = getReal(configuration, "post_margin_call_ratio", 70. / 100.);
	S_min = getReal(configuration, "print_asset_price_minimum", 0.);
	S_max = getReal(configuration, "print_asset_price_maximum", 300.);
	dS = getReal(configuration, "print_asset_price_step_size", 10.);
	N = getInt(configuration, "initial_number_of_timesteps", 2);
	coupon_payments = getInt(configuration, "coupon_payments", 4*T);
	margin_call_enabled = getBool(configuration, "margin_call_enabled", true);
	liquidation_enabled = getBool(configuration, "liquidation_enabled", true);
	use_dirty_loan_to_value_ratio = getBool(configuration, "use_dirty_loan_to_value_ratio", true);
	RectilinearGrid1 default_grid( S_0 * (Axis { 1e-4 } + Axis::special +  Axis { 1e4 }) );
	RectilinearGrid1 tmp = getGrid(configuration, "initial_grid", default_grid);
	grid = &tmp;

	#ifdef QUANT_PDE_STOCK_LOAN_EXPLICIT_ONLY
	explicit_impulses = true;
	#else
	explicit_impulses = getBool(configuration, "explicit_impulses", true);
	#endif

	// Ensure that N is picked appropriately
	N = max(N, coupon_payments + 1);

	// Print configuration file
	cerr << configuration << endl << endl;

	// Run and print results
	ResultsBuffer1 buffer(
		run,
		{ "Nodes", "Steps", "Mean Fixpt. Iterations" },
		kn, k0
	);
	buffer.setPrintGrid( RectilinearGrid1(Axis::range(S_min, dS, S_max)) );
	buffer.stream();

	return 0;
}

