////////////////////////////////////////////////////////////////////////////////
//
// Code to solve HJBQVI problems.
//
// Author: Parsiad Azimzadeh
//
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>

////////////////////////////////////////////////////////////////////////////////

#include <array>            // std::array
#include <chrono>           // std::chrono
#include <cmath>            // std::nan
#include <functional>       // std::function
#include <iostream>         // std::ostream, std::cout
#include <iomanip>          // std::setw
#include <initializer_list> // std::initializer_list
#include <limits>           // std::numeric_limits
#include <memory>           // std::unique_ptr
#include <numeric>          // std::accumulate
#include <string>           // std::string
#include <vector>           // std::vector

////////////////////////////////////////////////////////////////////////////////

namespace QuantPDE {
namespace HJBQVI {

template <unsigned N, unsigned M>
using ArrayFunction = std::array< Function<N>, M >;

struct ControlMethod {
	static constexpr char FULLY_IMPLICIT = 0;
	static constexpr char SEMI_LAGRANGIAN  = 1;
	static constexpr char EXPLICIT_IMPULSE = 2;
	static constexpr char FULLY_EXPLICIT = 1 | 2;
};

template <Index Dimension>
struct Result {

	const InterpolantWrapper<Dimension> solution; // u(t,x)

	const InterpolantWrapper<Dimension> stochastic_control;
	const InterpolantWrapper<Dimension> impulse_control;

	const Real mean_iterations;

	const Real execution_time_seconds;

	Result(
		const InterpolantWrapper<Dimension> &solution,

		const InterpolantWrapper<Dimension> &stochastic_control,
		const InterpolantWrapper<Dimension> &impulse_control,

		Real mean_iterations,

		Real execution_time_seconds
	) noexcept :
		solution(solution),

		stochastic_control(stochastic_control),
		impulse_control(impulse_control),

		mean_iterations(mean_iterations),

		execution_time_seconds(execution_time_seconds)
	{}

};

template <Index Dimension>
struct Options {

	const std::array<Real, Dimension> test_point;
	const int max_refinement;
	const int print_refinement;

	Options(
		std::array<Real, Dimension> test_point,
		int max_refinement = 0,
		int print_refinement = 0
	) noexcept:
		test_point(test_point),
		max_refinement(max_refinement),
		print_refinement(print_refinement)
	{}

};

template <Index Dimension>
struct Description {

	const RectilinearGrid<Dimension> spatial_grid;

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

	const bool refine_stochastic_control_grid;
	const bool refine_impulse_control_grid;

	const bool coefficients_time_independent;

	inline bool fully_implicit() const
	{ return handling == 0; }

	inline bool semi_lagrangian() const
	{ return handling & ControlMethod::SEMI_LAGRANGIAN; }

	inline bool explicit_impulse() const
	{ return handling & ControlMethod::EXPLICIT_IMPULSE; }

	inline bool fully_explicit() const
	{ return handling == ControlMethod::FULLY_EXPLICIT; }

	Description(
		const std::array<Axis, Dimension> spatial_axes,

		const std::array<Axis, 1> stochastic_control_axes,
		const std::array<Axis, 1> impulse_control_axes,

		Real expiry,
		Real discount,

		const ArrayFunction<Dimension+1, Dimension> volatility,
		const ArrayFunction<Dimension+2, Dimension> controlled_drift,
		const ArrayFunction<Dimension+1, Dimension> uncontrolled_drift,
		const Function<Dimension+2> controlled_continuous_flow,
		const Function<Dimension+1> uncontrolled_continuous_flow,
		const std::array<Function<Dimension+2>, Dimension> transition,
		const Function<Dimension+2> impulse_flow,
		const Function<Dimension+1> exit_function,

		int handling,
		int timesteps,

		bool refine_stochastic_control_grid = true,
		bool refine_impulse_control_grid = true,

		bool coefficients_time_independent = false
	) noexcept :
		spatial_grid(spatial_axes),

		stochastic_control_grid(stochastic_control_axes),
		impulse_control_grid(impulse_control_axes),

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

		refine_stochastic_control_grid(refine_stochastic_control_grid),
		refine_impulse_control_grid(refine_impulse_control_grid),

		coefficients_time_independent(coefficients_time_independent)
	{
		if(
			expiry >= std::numeric_limits<Real>::infinity()
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

	const Description<Dimension> hjbqvi;
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
		Real args[Dimension+2];
		for(int row = 0; row < refined_spatial_grid.size(); ++row) {
			Real total = 0.;

			// Get coordinates of point
			args[0] = time; // Time
			for(int d = 0; d < Dimension; ++d) {
				i[d] = (row / offsets[d]) %
						refined_spatial_grid[d].size();
				args[1+d] = refined_spatial_grid[d][i[d]];
			}
			args[1+Dimension] = q(row); // Control

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

				// Get volatility
				const Real v = packAndCall<Dimension+1>(
					hjbqvi.volatility[d],
					args
				);

				// Get drifts
				Real mu = packAndCall<Dimension+1>(
					hjbqvi.uncontrolled_drift[d],
					args
				);
				if(!hjbqvi.semi_lagrangian()) {
					mu += packAndCall<Dimension+2>(
						hjbqvi.controlled_drift[d],
						args
					);
				}

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
		Real args[Dimension+2];
		for(int row = 0; row < refined_spatial_grid.size(); ++row) {

			// Get coordinates of point
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

template <Index Dimension>
class ExplicitEvent : public EventBase {

	const Description<Dimension> &hjbqvi;
	const RectilinearGrid<Dimension> &refined_spatial_grid;
	const RectilinearGrid1 &refined_stochastic_control_grid;
	const RectilinearGrid1 &refined_impulse_control_grid;
	Vector &stochastic_control_vector;
	Vector &impulse_control_vector;
	Real time;
	const ReverseTimeIteration &stepper;
	std::vector<bool> &mask;

	int offsets[Dimension];

	template <typename V>
	Vector _doEvent(V &&vector) const {

		mask.clear();

		Vector best = refined_spatial_grid.vector();

		const Real dt = stepper.timestep();

		auto factory = refined_spatial_grid.defaultInterpolantFactory();
		auto u = factory.make(vector);

		Real args[Dimension+2];
		int i[Dimension];
		for(int row = 0; row < refined_spatial_grid.size(); ++row) {

		////////////////////////////////////////////////////////////////
		// begin row loop
		////////////////////////////////////////////////////////////////

		// Get coordinates of point
		args[0] = time; // Time
		for(int d = 0; d < Dimension; ++d) {
			i[d] = (row / offsets[d]) %
					refined_spatial_grid[d].size();
			args[1+d] = refined_spatial_grid[d][i[d]];
		}

		Real a = -std::numeric_limits<Real>::infinity();
		if(hjbqvi.semi_lagrangian()) {

			// Find optimal control
			for(
				int k = 0;
				k < refined_stochastic_control_grid.size();
				++k
			) {

				args[1+Dimension] =
					refined_stochastic_control_grid
					[0][k];

				std::array<Real, Dimension> new_state;
				for(int d = 0; d < Dimension; ++d) {
					const Real m = packAndCall<Dimension+2>(
						hjbqvi.controlled_drift[d],
						args
					);
					new_state[d] = args[1+d] + m * dt;
				}

				const Real flow = packAndCall<Dimension+2>(
					hjbqvi.controlled_continuous_flow,
					args
				);

				const Real new_value =
					u.interpolate(new_state)
					+ flow * dt
				;

				if(new_value > a) {
					a = new_value;
					stochastic_control_vector(row) =
							args[1+Dimension];
				}
			}

		} else {
			a = vector(row);
		}


		Real b = -std::numeric_limits<Real>::infinity();
		if(hjbqvi.explicit_impulse()) {

			// Find optimal control
			for(
				int k = 0;
				k < refined_impulse_control_grid.size();
				++k
			) {

				args[1+Dimension] =
						refined_impulse_control_grid
						[0][k];

				std::array<Real, Dimension> new_state;
				for(int d = 0; d < Dimension; ++d) {
					new_state[d] = packAndCall<Dimension+2>(
						hjbqvi.transition[d],
						args
					);
				}

				const Real flow = packAndCall<Dimension+2>(
					hjbqvi.impulse_flow,
					args
				);

				const Real new_value =
					u.interpolate(new_state)
					+ flow
				;

				if(new_value > b) {
					b = new_value;
					impulse_control_vector(row) =
							args[1+Dimension];
				}
			}

		}

		if(a >= b) {
			mask.push_back(false);
			best(row) = a;
		} else {
			mask.push_back(true);
			best(row) = b;
		}

		////////////////////////////////////////////////////////////////
		// end row loop
		////////////////////////////////////////////////////////////////

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

	template <typename H, typename R, typename S, typename I, typename T>
	ExplicitEvent(
		H &hjbqvi,
		R &refined_spatial_grid,
		S &refined_stochastic_control_grid,
		I &refined_impulse_control_grid,
		Vector &stochastic_control_vector,
		Vector &impulse_control_vector,
		Real time,
		T &stepper,
		std::vector<bool> &mask
	) noexcept :
		hjbqvi(hjbqvi),
		refined_spatial_grid(refined_spatial_grid),
		refined_stochastic_control_grid(
				refined_stochastic_control_grid),
		refined_impulse_control_grid(refined_impulse_control_grid),
		stochastic_control_vector(stochastic_control_vector),
		impulse_control_vector(impulse_control_vector),
		time(time),
		stepper(stepper),
		mask(mask)
	{
		// Space between ticks
		offsets[0] = 1;
		for(int d = 1; d < Dimension; ++d) {
			offsets[d] = offsets[d-1]
					* refined_spatial_grid[d].size();
		}
	}

};

template <Index Dimension>
Result<Dimension> solve(const Description<Dimension> &hjbqvi,
		int refinement = 0) {

	static_assert(Dimension <= 2 && Dimension >= 1,
			"Only dimensions 1 and 2 are currently supported.");

	const bool finite_horizon =
			hjbqvi.expiry < std::numeric_limits<Real>::infinity();

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

	Impulse<Dimension, 1> impulse(
		refined_spatial_grid,
		hjbqvi.impulse_flow,
		hjbqvi.transition
	);

	MinPolicyIteration<Dimension, 1> impulse_policy(
		refined_spatial_grid,
		refined_impulse_control_grid,
		impulse
	);

	stochastic_policy.setIteration(tolerance_iteration);
	impulse_policy.setIteration(tolerance_iteration);

	std::unique_ptr<ReverseTimeIteration> stepper;
	if(finite_horizon) {
		stepper = std::unique_ptr<ReverseTimeIteration>(
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

	ReverseBDFOne discretization(
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
	const Real dt = hjbqvi.expiry / timesteps;
	for(int e = 0; e < timesteps; ++e) {
		const Real time = e * dt;

		stepper->add(
			time,
			std::unique_ptr<EventBase>(
				new ExplicitEvent<Dimension>(
					hjbqvi,
					refined_spatial_grid,
					refined_stochastic_control_grid,
					refined_impulse_control_grid,
					stochastic_control_vector,
					impulse_control_vector,
					time,
					*stepper,
					mask
				)
			)
		);
	}

	// Linear system solver
	BiCGSTABSolver solver;

	// Bind exit function to expiry time to get payoff
	auto cauchy_data = curry<Dimension+1>(
		hjbqvi.exit_function,
		hjbqvi.expiry
	);

	// Timing
	Real seconds;
	auto start = std::chrono::steady_clock::now();

	// Solve
	auto u = iteration->solve(
		refined_spatial_grid,
		cauchy_data,
		*root,
		solver
	);

	// Timing
	auto end = std::chrono::steady_clock::now();
	auto diff = end - start;
	seconds = std::chrono::duration<Real>(diff).count();

	// Implicit stochastic control
	if(!hjbqvi.semi_lagrangian()) {
		stochastic_control_vector = controlled_operator.control(0);
	}

	// Implicit impulse control
	if(!hjbqvi.explicit_impulse()) {
		impulse_control_vector = impulse.control(0);
		mask = penalty.constraintMask();
	}

	// Mean iterations
	Real mean_iterations = std::nan("");
	if(!hjbqvi.fully_explicit()) {
		auto its = tolerance_iteration.iterations();
		mean_iterations = std::accumulate(its.begin(), its.end(), 0.)
				/ its.size();
	}

	// Apply mask
	for(int i = 0; i < refined_spatial_grid.size(); ++i) {
		if(mask[i]) {
			stochastic_control_vector(i) = std::nan("");
		} else {
			// Impulse control is not active here
			impulse_control_vector(i) = std::nan("");
		}
	}

	auto factory = refined_spatial_grid.defaultInterpolantFactory();

	// Return
	return Result<Dimension>(
		u,
		factory.make(stochastic_control_vector),
		factory.make(impulse_control_vector),
		mean_iterations,
		seconds
	);

}

template <Index Dimension>
void run(
	const Description<Dimension> &hjbqvi,
	const Options<Dimension> &opts,
	std::ostream &out = std::cout
) {

	out.precision(12);

	const int spacing = 23;
	auto space = [=] () { return std::setw(spacing); };

	// Headers
	out
		<< space() << "Spatial Nodes"
		<< space() << "Stochastic Ctrl Nodes"
		<< space() << "Impulse Ctrl Nodes"
		<< space() << "Timesteps"
		<< space() << "Mean Iterations"
		<< space() << "Value"
		<< space() << "Change"
		<< space() << "Ratio"
		<< space() << "Execution Time (sec)"
		<< std::endl
	;

	Real
		previousValue = nan(""),
		previousChange = nan(""),
		value,
		change,
		ratio
	;

	for(
		int
			refinement = 0,
			timesteps = hjbqvi.timesteps,
			spatial_nodes = hjbqvi.spatial_grid.size(),
			q_nodes = hjbqvi.stochastic_control_grid.size(),
			zeta_nodes = hjbqvi.impulse_control_grid.size()
		;
			refinement <= opts.max_refinement
		;
			++refinement,
			timesteps *= 2,
			spatial_nodes = spatial_nodes * 2 - 1,
			q_nodes = q_nodes * 2 - 1,
			zeta_nodes = zeta_nodes * 2 - 1
	) {
		auto result = solve<Dimension>(hjbqvi, refinement);

		// Get value of function at the test point
		value = result.solution.interpolate(opts.test_point);
		change = value - previousValue;
		ratio = previousChange / change;
		previousValue = value;
		previousChange = change;

		// Print
		out
			<< space() << spatial_nodes
			<< space() << q_nodes
			<< space() << zeta_nodes
			<< space() << timesteps
			<< space() << result.mean_iterations
			<< space() << value
			<< space() << change
			<< space() << ratio
			<< space() << result.execution_time_seconds
			<< std::endl
		;

		if(refinement == opts.max_refinement) {
			// Print header
			out << std::endl; // Extra spacing
			for(int d = 0; d < Dimension; ++d) {
				out << space() << ("x_" + std::to_string(d+1));
			}
			out
				<< space() << "Value u(t=0, x)"
				<< space() << "Stochastic Control"
				<< space() << "Impulse Control"
				<< std::endl
			;

			auto print_grid = hjbqvi.spatial_grid.refined(
					opts.print_refinement);
			for(auto node : print_grid) {
				for(int d = 0; d < Dimension; ++d) {
					out << space() << node[d];
				}
				out
					<< space() << result.solution
							.interpolate(node)
					<< space() << result.stochastic_control
							.interpolate(node)
					<< space() << result.impulse_control
							.interpolate(node)
					<< std::endl
				;
			}
		}
	}
}

} // namespace HJBQVI
} // namespace QuantPDE

// TODO: Make Axis.hpp safe

