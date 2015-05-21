////////////////////////////////////////////////////////////////////////////////
//
// Author: Parsiad Azimzadeh
//
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>

using namespace QuantPDE;

////////////////////////////////////////////////////////////////////////////////

#include <array>      // array
#include <chrono>     // duration
#include <cmath>      // fabs
#include <functional> // std::function
#include <iostream>   // cerr
#include <limits>     // numeric_limits
#include <memory>     // unique_ptr
#include <string>     // string

using namespace std;

////////////////////////////////////////////////////////////////////////////////

template <unsigned N, unsigned M>
using ArrayFunction = std::function<
	NaryFunctionSignature<
		std::array<Real, M>,
		N,
		Real
	>
>;

enum class ControlHandling {
	IMPLICIT,
	EXPLICIT,
	MIXED
};

template <Index Dimension>
struct HJBQVI {

	const RectilinearGrid<Dimension> spatial_grid;

	// TODO: Extend this to arbitrary dimensions
	const RectilinearGrid1 stochastic_control_grid;
	const RectilinearGrid1 impulse_control_grid;

	const double expiry;
	const double discount;

	const ArrayFunction<Dimension+1, Dimension> volatility;
	const ArrayFunction<Dimension+2, Dimension> controlled_drift;
	const ArrayFunction<Dimension+1, Dimension> uncontrolled_drift;
	const Function<Dimension+2> controlled_continuous_flow;
	const Function<Dimension+1> uncontrolled_continuous_flow;
	const ArrayFunction<Dimension+2, Dimension> transition;
	const Function<Dimension+2> impulse_flow;
	const Function<Dimension+1> exit_function;

	const ControlHandling handling;
	const int timesteps;

	HJBQVI(
		const RectilinearGrid<Dimension> spatial_grid,

		// TODO: Extend this to arbitrary dimensions
		const RectilinearGrid1 stochastic_control_grid,
		const RectilinearGrid1 impulse_control_grid,

		const double expiry,
		const double discount,

		const ArrayFunction<Dimension+1, Dimension> volatility,
		const ArrayFunction<Dimension+2, Dimension> controlled_drift,
		const ArrayFunction<Dimension+1, Dimension> uncontrolled_drift,
		const Function<Dimension+2> controlled_continuous_flow,
		const Function<Dimension+1> uncontrolled_continuous_flow,
		const ArrayFunction<Dimension+2, Dimension> transition,
		const Function<Dimension+2> impulse_flow,
		const Function<Dimension+1> exit_function,

		const ControlHandling handling,
		const int timesteps
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
		timesteps(timesteps)
	{}

};

template <int Dimension>
struct ControlledOperator final : public
		RawControlledLinearSystem<Dimension, 1> {

	const HJBQVI<Dimension> hjbqvi;
	int offsets[Dimension];

	template <typename H>
	ControlledOperator(H &hjbqvi) noexcept : hjbqvi(hjbqvi) {
		// Space between ticks
		offsets[0] = 1;
		for(int d = 1; d < Dimension; ++d) {
			offsets[d] = offsets[d-1]
					* hjbqvi.spatial_grid[d].size();
		}
	}

	virtual Matrix A(double time) {
		Matrix A = hjbqvi.spatial_grid.matrix();
		A.reserve(
			IntegerVector::Constant(
				hjbqvi.spatial_grid.size(),
				1 + 2 * Dimension
			)
		); // Reserve nonzero entries per row

		// Control as a vector
		const Vector &q = hjbqvi.handling == ControlHandling::IMPLICIT ?
			this->control(0)
			: hjbqvi.spatial_grid.zero()
		;

		// Iterate through points on grid
		int i[Dimension];
		for(int row = 0; row < hjbqvi.spatial_grid.size(); ++row) {
			double total = 0.;

			// Get coordinates of point
			Real args[Dimension+2];
			args[0] = time; // Time
			for(int d = 0; d < Dimension; ++d) {
				i[d] = (row / offsets[d])
						% hjbqvi.spatial_grid[d].size();
				args[1+d] = hjbqvi.spatial_grid[d][i[d]];
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
						hjbqvi.spatial_grid[d].size()
						- 1) { continue; }

				const Axis &x = hjbqvi.spatial_grid[d];

				const double
					dxb = x[ i[d]     ] - x[ i[d] - 1 ],
					dxc = x[ i[d] + 1 ] - x[ i[d] - 1 ],
					dxf = x[ i[d] + 1 ] - x[ i[d]     ]
				;

				// Local volatility
				const double v = volatility[d];

				// Local drift
				const double mu = (hjbqvi.handling ==
					ControlHandling::IMPLICIT ?
						controlled[d]
						: 0.
				) + uncontrolled[d];

				const double vv = v * v;

				const double alpha_common = vv / dxb / dxc;
				const double  beta_common = vv / dxf / dxc;

				// Central
				double alpha = alpha_common - mu / dxc;
				double  beta =  beta_common + mu / dxc;
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

	virtual Vector b(double time) {
		Vector b = hjbqvi.spatial_grid.vector();

		// Control as a vector
		const Vector &q = hjbqvi.handling == ControlHandling::IMPLICIT ?
			this->control(0)
			: hjbqvi.spatial_grid.zero()
		;


		// Iterate through points on grid
		int i[Dimension];
		for(int row = 0; row < hjbqvi.spatial_grid.size(); ++row) {

			// Get coordinates of point
			Real args[Dimension+2];
			args[0] = time; // Time
			for(int d = 0; d < Dimension; ++d) {
				i[d] = (row / offsets[d])
						% hjbqvi.spatial_grid[d].size();
				args[1+d] = hjbqvi.spatial_grid[d][i[d]];
			}
			args[1+Dimension] = q(row); // Control

			// Get flows
			double controlled = packAndCall<Dimension+2>(
				hjbqvi.controlled_continuous_flow,
				args
			);
			double uncontrolled = packAndCall<Dimension+1>(
				hjbqvi.uncontrolled_continuous_flow,
				args
			);

			b(row) = (hjbqvi.handling == ControlHandling::IMPLICIT ?
				controlled
				: 0.
			) + uncontrolled;

		}

		return b;
	}

	virtual bool isATheSame() const {
		// TODO: Check
		return false;
	}

};

// 1D
template <Index Dimension>
std::unique_ptr<ControlledLinearSystemBase> make_impulse(
		const HJBQVI<Dimension> &hjbqvi, int refinement) {

	static_assert(Dimension == 1, "Dimension mismatch");

	return std::unique_ptr<ControlledLinearSystemBase>(
		(ControlledLinearSystemBase*)
		new Impulse1_1(
			hjbqvi.spatial_grid,
			[&hjbqvi] (double t, double x, double zeta) {
				Real args[3];
				args[0] = hjbqvi.expiry;
				args[1] = x;
				args[2] = zeta;
				return packAndCall<1+2>(
					hjbqvi.impulse_flow,
					args
				);
			},
			[&hjbqvi] (double t, double x, double zeta) {
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
std::unique_ptr<ControlledLinearSystemBase> make_impulse<2>(
		const HJBQVI<2> &hjbqvi, int refinement) {

	return std::unique_ptr<ControlledLinearSystemBase>(
		(ControlledLinearSystemBase*)
		new Impulse2_1(
			hjbqvi.spatial_grid,
			[&hjbqvi] (double t, double x, double y,
					double zeta) {
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
			[&hjbqvi] (double t, double x, double y,
					double zeta) {
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
			[&hjbqvi] (double t, double x, double y,
					double zeta) {
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

template <Index Dimension>
InterpolantWrapper<Dimension> make_solution(const HJBQVI<Dimension> &hjbqvi,
		Iteration *iteration, IterationNode *root) {

	static_assert(Dimension == 1, "Dimension mismatch");

	BiCGSTABSolver solver;

	return iteration->solve(
		hjbqvi.spatial_grid,
		[=] (double x) {
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

template <>
InterpolantWrapper<2> make_solution(const HJBQVI<2> &hjbqvi,
		Iteration *iteration, IterationNode *root) {

	BiCGSTABSolver solver;

	return iteration->solve(
		hjbqvi.spatial_grid,
		[=] (double x, double y) {
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

template <Index Dimension>
void solve(const HJBQVI<Dimension> &hjbqvi, int refinement) {

	static_assert(Dimension <= 2 && Dimension >= 1,
			"Only dimensions 1 and 2 are currently supported.");

	const bool finite_horizon =
			hjbqvi.expiry < std::numeric_limits<double>::infinity();

	ToleranceIteration tolerance_iteration;

	ControlledOperator<Dimension> controlled_operator(hjbqvi);

	MinPolicyIteration<Dimension, 1> stochastic_policy(
		hjbqvi.spatial_grid,
		hjbqvi.stochastic_control_grid,
		controlled_operator
	);

	std::unique_ptr<ControlledLinearSystemBase> impulse =
			make_impulse<Dimension>(hjbqvi, refinement);

	MinPolicyIteration<Dimension, 1> impulse_policy(
		hjbqvi.spatial_grid,
		hjbqvi.impulse_control_grid,
		*impulse
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
				hjbqvi.expiry / hjbqvi.timesteps
			)
		);

		if(hjbqvi.handling != ControlHandling::EXPLICIT) {
			stepper->setInnerIteration(tolerance_iteration);
		}
	}

	LinearSystem *discretize;
	if(hjbqvi.handling == ControlHandling::IMPLICIT) {
		discretize = &stochastic_policy;
	} else {
		// Semi-Lagrangian
		discretize = &controlled_operator;
	}

	ReverseBDFOne<Dimension> discretization(
		hjbqvi.spatial_grid,
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
		hjbqvi.spatial_grid,
		*penalized,
		impulse_policy
	);
	penalty.setIteration(tolerance_iteration);

	// Pick root
	IterationNode *root;
	if(hjbqvi.handling == ControlHandling::EXPLICIT) {
		// Explicit impulse
		root = penalized;
	} else {
		// Implicit impulse
		root = &penalty;
	}

	// TODO: Events

	// Timing
	auto start = std::chrono::steady_clock::now();
	double seconds;

	auto u = make_solution(hjbqvi, iteration, root);

	// Timing
	auto end = chrono::steady_clock::now();
	auto diff = end - start;
	seconds = chrono::duration<double>(diff).count();

	// Print solution
	std::cout << accessor( hjbqvi.spatial_grid, u ) << std::endl;

}

int main(int argc, char **argv) {

	const RectilinearGrid1 spatial_grid = RectilinearGrid1(
		Axis::cluster(
			-6., // Left-hand boundary
			0.,  // Feature to cluster around
			+6., // Right-hand boundary
			32,  // Number of points
			10.  // Clustering intensity
		)
	);

	constexpr int dimension = 1;
	HJBQVI<dimension> hjbqvi(
		// Spatial grid
		spatial_grid,

		// Stochastic control grid
		RectilinearGrid1(
			Axis::uniform(
				0.,   // Left-hand boundary
				0.05, // Right-hand boundary
				32    // Number of points
			)
		),

		// Impulse control grid
		spatial_grid,

		// Expiry
		10., //numeric_limits<double>::infinity(),

		// Discount
		0.02,

		// Volatility
		[] (Real t, Real x) {
			const std::array<Real, 1> res { 0.3 };
			return res;
		},

		// Controlled drift
		[] (Real t, Real x, Real q) {
			const Real a = 0.25;
			const std::array<Real, 1> res { -a * q };
			return res;
		},

		// Uncontrolled drift
		[] (Real t, Real x) {
			const std::array<Real, 1> res { 0. };
			return res;
		},

		// Controlled continuous flow
		[] (Real t, Real x, Real q) {
			const Real b = 3.;
			return -q * q * b;
		},

		// Uncontrolled continuous flow
		[] (Real t, Real x) {
			const Real m = 0.;
			const Real tmp = max(x - m, 0.);
			return -tmp * tmp;
		},

		// Transition
		[] (Real t, Real x, Real x_new) {
			const std::array<Real, 1> res { x_new };
			return res;
		},

		// Impulse flow
		[] (Real t, Real x, Real x_new) {
			const Real lambda = 1., c = 0.;
			const Real zeta = x_new - x;
			return - lambda * fabs(zeta) - c;
		},

		// Exit function
		[] (Real t, Real x) { return 0.; },

		// Implicit method?
		ControlHandling::IMPLICIT,

		// Number of timesteps
		32
	);

	solve<1>(hjbqvi, 0);

	return 0;

}

