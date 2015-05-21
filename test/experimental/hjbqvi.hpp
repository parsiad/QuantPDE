////////////////////////////////////////////////////////////////////////////////
//
// Code to solve HJBQVI problems.
//
// Author: Parsiad Azimzadeh
//
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>

using namespace QuantPDE;

////////////////////////////////////////////////////////////////////////////////

#include <array>      // array
#include <chrono>     // chrono::duration
#include <cmath>      // nan
#include <functional> // function
#include <getopt.h>   // getopt_long
#include <iostream>   // cerr, cout
#include <limits>     // numeric_limits
#include <memory>     // unique_ptr
#include <vector>     // vector

using namespace std;

////////////////////////////////////////////////////////////////////////////////

template <unsigned N, unsigned M>
using ArrayFunction = function<
	NaryFunctionSignature<
		array<Real, M>,
		N,
		Real
	>
>;

template <Index Dimension>
struct Result {

	const InterpolantWrapper<Dimension> solution; // u(t,x)

	const InterpolantWrapper<Dimension> stochastic_control;
	const InterpolantWrapper<Dimension> impulse_control;

	const double execution_time_seconds;

	Result(
		const InterpolantWrapper<Dimension> &solution,

		const InterpolantWrapper<Dimension> &stochastic_control,
		const InterpolantWrapper<Dimension> &impulse_control,

		double execution_time_seconds
	) noexcept :
		solution(solution),

		stochastic_control(stochastic_control),
		impulse_control(impulse_control),

		execution_time_seconds(execution_time_seconds)
	{}

};

struct ControlMethod {
	static constexpr char FULLY_IMPLICIT = 0;
	static constexpr char SEMI_LAGRANGIAN  = 1;
	static constexpr char EXPLICIT_IMPULSE = 2;
	static constexpr char FULLY_EXPLICIT = 1 | 2;
};

template <Index Dimension>
struct HJBQVI {

	const RectilinearGrid<Dimension> spatial_grid;

	// TODO: Extend this to arbitrary dimensions
	const RectilinearGrid1 stochastic_control_grid;
	const RectilinearGrid1 impulse_control_grid;

	const Real expiry;
	const Real discount;

	const ArrayFunction<Dimension+1, Dimension> volatility;
	const ArrayFunction<Dimension+2, Dimension> controlled_drift;
	const ArrayFunction<Dimension+1, Dimension> uncontrolled_drift;
	const Function<Dimension+2> controlled_continuous_flow;
	const Function<Dimension+1> uncontrolled_continuous_flow;
	const ArrayFunction<Dimension+2, Dimension> transition;
	const Function<Dimension+2> impulse_flow;
	const Function<Dimension+1> exit_function;

	const char handling;
	const int timesteps;

	const bool coefficients_time_independent;

	inline bool fully_implicit() const
	{ return handling == 0; }

	inline bool semi_lagrangian() const
	{ return handling & ControlMethod::SEMI_LAGRANGIAN; }

	inline bool explicit_impulse() const
	{ return handling & ControlMethod::EXPLICIT_IMPULSE; }

	inline bool fully_explicit() const
	{ return handling == ControlMethod::FULLY_EXPLICIT; }

	HJBQVI(
		const RectilinearGrid<Dimension> spatial_grid,

		// TODO: Extend this to arbitrary dimensions
		const RectilinearGrid1 stochastic_control_grid,
		const RectilinearGrid1 impulse_control_grid,

		Real expiry,
		Real discount,

		const ArrayFunction<Dimension+1, Dimension> volatility,
		const ArrayFunction<Dimension+2, Dimension> controlled_drift,
		const ArrayFunction<Dimension+1, Dimension> uncontrolled_drift,
		const Function<Dimension+2> controlled_continuous_flow,
		const Function<Dimension+1> uncontrolled_continuous_flow,
		const ArrayFunction<Dimension+2, Dimension> transition,
		const Function<Dimension+2> impulse_flow,
		const Function<Dimension+1> exit_function,

		int handling,
		int timesteps,

		bool coefficients_time_independent = false
	) noexcept :
		spatial_grid(spatial_grid),

		stochastic_control_grid(stochastic_control_grid),
		impulse_control_grid(impulse_control_grid),

		expiry(expiry),
		discount(discount),

		volatility(volatility),
		controlled_drift(controlled_drift),
		uncontrolled_drift(uncontrolled_drift),
		controlled_continuous_flow(controlled_continuous_flow),
		uncontrolled_continuous_flow(uncontrolled_continuous_flow),
		transition(transition),
		impulse_flow(impulse_flow),
		exit_function(exit_function),

		handling(handling),
		timesteps(timesteps),

		coefficients_time_independent(coefficients_time_independent)
	{
		if(
			expiry >= numeric_limits<Real>::infinity()
			&& !fully_implicit()
		) {
			throw "error: only an implicit method can be used for"
					" infinite-horizon problems";
		}
	}

};

template <int Dimension>
struct ControlledOperator final : public
		RawControlledLinearSystem<Dimension, 1> {

	const HJBQVI<Dimension> hjbqvi;
	RectilinearGrid<Dimension> refined_spatial_grid;
	int offsets[Dimension];

	template <typename H, typename R>
	ControlledOperator(
		H &hjbqvi,
		R &refined_spatial_grid
	) noexcept :
		hjbqvi(hjbqvi),
		refined_spatial_grid(refined_spatial_grid)
	{
		// Space between ticks
		offsets[0] = 1;
		for(int d = 1; d < Dimension; ++d) {
			offsets[d] = offsets[d-1]
					* refined_spatial_grid[d].size();
		}
	}

	virtual Matrix A(Real time) {
		Matrix A = refined_spatial_grid.matrix();
		A.reserve(
			IntegerVector::Constant(
				refined_spatial_grid.size(),
				1 + 2 * Dimension
			)
		); // Reserve nonzero entries per row

		// Control as a vector
		const Vector &q = hjbqvi.semi_lagrangian() ?
				refined_spatial_grid.zero() : this->control(0);

		// Iterate through points on grid
		int i[Dimension];
		for(int row = 0; row < refined_spatial_grid.size(); ++row) {
			Real total = 0.;

			// Get coordinates of point
			Real args[Dimension+2];
			args[0] = time; // Time
			for(int d = 0; d < Dimension; ++d) {
				i[d] = (row / offsets[d]) %
						refined_spatial_grid[d].size();
				args[1+d] = refined_spatial_grid[d][i[d]];
			}
			args[1+Dimension] = q(row); // Control

			// Get volatility
			auto volatility = packAndCall<Dimension+1>(
				hjbqvi.volatility,
				args
			);

			// Get drifts
			auto controlled = packAndCall<Dimension+2>(
				hjbqvi.controlled_drift,
				args
			);
			auto uncontrolled = packAndCall<Dimension+1>(
				hjbqvi.uncontrolled_drift,
				args
			);

			for(int d = 0; d < Dimension; ++d) {
				// Skip boundary points
				if(i[d] == 0 || i[d] ==
						refined_spatial_grid[d].size()
						- 1) { continue; }

				const Axis &x = refined_spatial_grid[d];

				const Real
					dxb = x[ i[d]     ] - x[ i[d] - 1 ],
					dxc = x[ i[d] + 1 ] - x[ i[d] - 1 ],
					dxf = x[ i[d] + 1 ] - x[ i[d]     ]
				;

				// Local volatility
				const Real v = volatility[d];

				// Local drift
				const Real mu = (hjbqvi.semi_lagrangian() ? 0.
						: controlled[d])
						+ uncontrolled[d];

				const Real vv = v * v;

				const Real alpha_common = vv / dxb / dxc;
				const Real  beta_common = vv / dxf / dxc;

				// Central
				Real alpha = alpha_common - mu / dxc;
				Real  beta =  beta_common + mu / dxc;
				if(alpha < 0.) {
					alpha = alpha_common;
					 beta =  beta_common + mu / dxf;
				} else if(beta < 0.) {
					alpha = alpha_common - mu / dxb;
					 beta =  beta_common;
				}

				A.insert(row, row - offsets[d]) = -alpha;
				A.insert(row, row + offsets[d]) = - beta;

				total += alpha + beta;
			}

			A.insert(row, row) = total + hjbqvi.discount;
		}

		A.makeCompressed();
		return A;
	}

	virtual Vector b(Real time) {
		Vector b = refined_spatial_grid.vector();

		// Control as a vector
		const Vector &q = hjbqvi.semi_lagrangian() ?
				refined_spatial_grid.zero() : this->control(0);

		// Iterate through points on grid
		int i[Dimension];
		for(int row = 0; row < refined_spatial_grid.size(); ++row) {

			// Get coordinates of point
			Real args[Dimension+2];
			args[0] = time; // Time
			for(int d = 0; d < Dimension; ++d) {
				i[d] = (row / offsets[d]) %
						refined_spatial_grid[d].size();
				args[1+d] = refined_spatial_grid[d][i[d]];
			}
			args[1+Dimension] = q(row); // Control

			// Get flows
			Real controlled = packAndCall<Dimension+2>(
				hjbqvi.controlled_continuous_flow,
				args
			);
			Real uncontrolled = packAndCall<Dimension+1>(
				hjbqvi.uncontrolled_continuous_flow,
				args
			);

			b(row) = (hjbqvi.semi_lagrangian() ? 0 : controlled)
					+ uncontrolled;
		}

		return b;
	}

	virtual bool isATheSame() const {
		return hjbqvi.coefficients_time_independent;
	}

};

// 1D
template <Index Dimension>
unique_ptr<ControlledLinearSystemBase> make_impulse(
	const HJBQVI<Dimension> &hjbqvi,
	const RectilinearGrid<Dimension> &refined_spatial_grid,
	int refinement
) {

	static_assert(Dimension == 1, "Dimension mismatch");

	return unique_ptr<ControlledLinearSystemBase>(
		(ControlledLinearSystemBase*)
		new Impulse1_1(
			refined_spatial_grid,
			[&hjbqvi] (Real t, Real x, Real zeta) {
				Real args[3];
				args[0] = hjbqvi.expiry;
				args[1] = x;
				args[2] = zeta;
				return packAndCall<1+2>(
					hjbqvi.impulse_flow,
					args
				);
			},
			[&hjbqvi] (Real t, Real x, Real zeta) {
				Real args[3];
				args[0] = t;
				args[1] = x;
				args[2] = zeta;
				auto res = packAndCall<1+2>(
					hjbqvi.transition,
					args
				);
				return res[0];
			}
		)
	);

}

// 2D
template <>
unique_ptr<ControlledLinearSystemBase> make_impulse<2>(
	const HJBQVI<2> &hjbqvi,
	const RectilinearGrid2 &refined_spatial_grid,
	int refinement
) {

	return unique_ptr<ControlledLinearSystemBase>(
		(ControlledLinearSystemBase*)
		new Impulse2_1(
			refined_spatial_grid,
			[&hjbqvi] (Real t, Real x, Real y, Real zeta) {
				Real args[4];
				args[0] = hjbqvi.expiry;
				args[1] = x;
				args[2] = y;
				args[3] = zeta;
				return packAndCall<2+2>(
					hjbqvi.impulse_flow,
					args
				);
			},
			[&hjbqvi] (Real t, Real x, Real y, Real zeta) {
				Real args[4];
				args[0] = t;
				args[1] = x;
				args[2] = y;
				args[3] = zeta;
				auto res = packAndCall<2+2>(
					hjbqvi.transition,
					args
				);
				return res[0];
			},
			[&hjbqvi] (Real t, Real x, Real y, Real zeta) {
				Real args[4];
				args[0] = t;
				args[1] = x;
				args[2] = y;
				args[3] = zeta;
				auto res = packAndCall<2+2>(
					hjbqvi.transition,
					args
				);
				return res[1];
			}
		)
	);

}

// 1D
template <Index Dimension>
InterpolantWrapper<Dimension> make_solution(
	const HJBQVI<Dimension> &hjbqvi,
	const RectilinearGrid<Dimension> &refined_spatial_grid,
	Iteration *iteration,
	IterationNode *root
) {

	static_assert(Dimension == 1, "Dimension mismatch");

	BiCGSTABSolver solver;

	return iteration->solve(
		refined_spatial_grid,
		[=] (Real x) {
			Real args[2];
			args[0] = hjbqvi.expiry;
			args[1] = x;
			return packAndCall<1+1>(
				hjbqvi.exit_function,
				args
			);
		},
		*root,
		solver
	);

}

// 2D
template <>
InterpolantWrapper<2> make_solution<2>(
	const HJBQVI<2> &hjbqvi,
	const RectilinearGrid2 &refined_spatial_grid,
	Iteration *iteration,
	IterationNode *root) {

	BiCGSTABSolver solver;

	return iteration->solve(
		refined_spatial_grid,
		[=] (Real x, Real y) {
			Real args[2];
			args[0] = hjbqvi.expiry;
			args[1] = x;
			args[2] = x;
			return packAndCall<2+1>(
				hjbqvi.exit_function,
				args
			);
		},
		*root,
		solver
	);

}

// 1D
template <Index Dimension>
Transform<Dimension> make_event(
	const HJBQVI<Dimension> &hjbqvi,
	const RectilinearGrid1 &refined_stochastic_control_grid,
	const RectilinearGrid1 &refined_impulse_control_grid,
	Real t, Real dt,
	Vector &stochastic_control_vector,
	Vector &impulse_control_vector,
	std::vector<bool> &mask
) {

	static_assert(Dimension == 1, "Dimension mismatch");

	const Axis &Q = refined_stochastic_control_grid[0];

	auto semi_lagrangian = [=, &hjbqvi, &Q, &stochastic_control_vector]
			(const Interpolant1 &u, Real x) {
		static int k;
		if(x == hjbqvi.spatial_grid[0][0]) { k = 0; }

		Real best = -numeric_limits<Real>::infinity();

		for(int i = 0; i < Q.size(); ++i) {
			const Real q = Q[i];
			const Real new_value =
				u(x + hjbqvi.controlled_drift(t, x, q)[0] * dt)
				+ hjbqvi.controlled_continuous_flow(t, x, q)*dt
			;
			if(new_value > best) {
				best = new_value;
				stochastic_control_vector(k) = q;
			}
		}

		++k;

		return best;
	};

	const Axis &Z = refined_impulse_control_grid[0];

	auto explicit_impulse = [=, &hjbqvi, &Z, &impulse_control_vector]
			(const Interpolant1 &u, Real x) {
		static int k;
		if(x == hjbqvi.spatial_grid[0][0]) { k = 0; }

		Real best = -numeric_limits<Real>::infinity();

		for(int i = 0; i < Z.size(); ++i) {
			const Real zeta = Z[i];
			const Real new_value =
				u(hjbqvi.transition(t, x, zeta)[0])
				+ hjbqvi.impulse_flow(t, x, zeta)
			;
			if(new_value > best) {
				best = new_value;
				impulse_control_vector(k) = zeta;
			}
		}

		++k;

		return best;
	};


	if(hjbqvi.fully_explicit()) {
		return [=, &hjbqvi, &mask] (const Interpolant1 &u, Real x) {
			// Clear mask
			if(x == hjbqvi.spatial_grid[0][0]) {
				mask.clear();
			}

			Real a = semi_lagrangian(u, x);
			Real b = explicit_impulse(u, x);
			if(a >= b) {
				mask.push_back(false);
				return a;
			}
			mask.push_back(true);
			return b;
		};
	}

	if(hjbqvi.semi_lagrangian()) {
		return semi_lagrangian;
	}

	if(hjbqvi.explicit_impulse()) {
		return [=, &hjbqvi, &mask] (const Interpolant1 &u, Real x) {
			// Clear mask
			if(x == hjbqvi.spatial_grid[0][0]) {
				mask.clear();
			}

			Real a = u(x);
			Real b = explicit_impulse(u, x);
			if(a >= b) {
				mask.push_back(false);
				return a;
			}
			mask.push_back(true);
			return b;
		};
	}

	throw "error: unexpected error";

}

template <Index Dimension>
Result<Dimension> solve(const HJBQVI<Dimension> &hjbqvi, int refinement = 0) {

	static_assert(Dimension <= 2 && Dimension >= 1,
			"Only dimensions 1 and 2 are currently supported.");

	const bool finite_horizon =
			hjbqvi.expiry < numeric_limits<Real>::infinity();

	// Refine grid
	auto refined_spatial_grid = hjbqvi.spatial_grid.refined(refinement);
	auto refined_stochastic_control_grid =
			hjbqvi.stochastic_control_grid.refined(refinement);
	auto refined_impulse_control_grid =
			hjbqvi.impulse_control_grid.refined(refinement);

	Vector stochastic_control_vector = refined_spatial_grid.vector();
	Vector impulse_control_vector = refined_spatial_grid.vector();

	std::vector<bool> mask;
	mask.reserve(refined_spatial_grid.size());

	// Refine timesteps
	int timesteps = hjbqvi.timesteps;
	for(int i = 0; i < refinement; ++i) { timesteps *= 2; }

	ToleranceIteration tolerance_iteration;

	ControlledOperator<Dimension> controlled_operator(
		hjbqvi,
		refined_spatial_grid
	);

	MinPolicyIteration<Dimension, 1> stochastic_policy(
		refined_spatial_grid,
		refined_stochastic_control_grid,
		controlled_operator
	);

	unique_ptr<ControlledLinearSystemBase> impulse =
			make_impulse<Dimension>(hjbqvi, refined_spatial_grid,
			refinement);

	MinPolicyIteration<Dimension, 1> impulse_policy(
		refined_spatial_grid,
		refined_impulse_control_grid,
		*impulse
	);

	stochastic_policy.setIteration(tolerance_iteration);
	impulse_policy.setIteration(tolerance_iteration);

	unique_ptr<ReverseTimeIteration> stepper;
	if(finite_horizon) {
		stepper = unique_ptr<ReverseTimeIteration>(
			(ReverseTimeIteration*)
			new ReverseConstantStepper(
				0.,
				hjbqvi.expiry,
				hjbqvi.expiry / timesteps
			)
		);

		if(!hjbqvi.fully_explicit()) {
			stepper->setInnerIteration(tolerance_iteration);
		}
	}

	LinearSystem *discretize;
	if(hjbqvi.semi_lagrangian()) {
		// Semi-Lagrangian
		discretize = &controlled_operator;
	} else {
		discretize = &stochastic_policy;
	}

	ReverseBDFOne<Dimension> discretization(
		refined_spatial_grid,
		*discretize
	);

	// Finite or infinite horizon
	Iteration *iteration;
	IterationNode *penalized;
	if(finite_horizon) {
		discretization.setIteration(*stepper);

		iteration = stepper.get();
		penalized = &discretization;
	} else {
		iteration = &tolerance_iteration;
		penalized = &stochastic_policy;
	}

	// Penalty method
	PenaltyMethod penalty(
		refined_spatial_grid,
		*penalized,
		impulse_policy
	);
	penalty.setIteration(tolerance_iteration);

	// Pick root
	IterationNode *root;
	if(hjbqvi.explicit_impulse()) {
		// Explicit impulse
		root = penalized;
	} else {
		// Implicit impulse
		root = &penalty;
	}

	// Add events
	if(!hjbqvi.fully_implicit()) {
		const Real dt = hjbqvi.expiry / timesteps;
		for(int e = 0; e < timesteps; ++e) {
			const Real t = e * dt;
			stepper->add(
				t,
				make_event(
					hjbqvi,
					refined_stochastic_control_grid,
					refined_impulse_control_grid,
					t, dt,
					stochastic_control_vector,
					impulse_control_vector,
					mask
				),
				refined_spatial_grid
			);
		}
	}

	// Timing
	auto start = chrono::steady_clock::now();
	Real seconds;

	auto u = make_solution(hjbqvi, refined_spatial_grid, iteration, root);

	// Timing
	auto end = chrono::steady_clock::now();
	auto diff = end - start;
	seconds = chrono::duration<Real>(diff).count();

	// Implicit stochastic control
	if(!hjbqvi.semi_lagrangian()) {
		stochastic_control_vector = controlled_operator.control(0);
	}

	// Implicit impulse control
	if(!hjbqvi.explicit_impulse()) {
		if(Dimension == 1) {
			impulse_control_vector = ((Impulse1_1 *)impulse.get())
					->control(0);
		} // TODO: Other dimensions
		mask = penalty.constraintMask();
	}

	// Apply mask
	for(int i = 0; i < refined_spatial_grid.size(); ++i) {
		if(mask[i]) {
			stochastic_control_vector(i) = nan("");
		} else {
			// Impulse control is not active here
			impulse_control_vector(i) = nan("");
		}
	}

	auto factory = refined_spatial_grid.defaultInterpolantFactory();

	// Return
	return Result<Dimension>(
		u,
		factory.make(stochastic_control_vector),
		factory.make(impulse_control_vector),
		seconds
	);

}

template <Index Dimension>
int hjbqvi_main(const HJBQVI<Dimension> &hjbqvi) {

	try {

		auto result = solve<Dimension>(hjbqvi, 0);

		cout
			<< accessor(
				hjbqvi.spatial_grid,
				result.stochastic_control
			)
			<< endl
			<< accessor(
				hjbqvi.spatial_grid,
				result.impulse_control
			)
			<< endl
			<< accessor(
				hjbqvi.spatial_grid,
				result.solution
			)
			<< endl
		;

	} catch(...) {
		// TODO: Error reporting
		return 1;
	}

	return 0;

}

