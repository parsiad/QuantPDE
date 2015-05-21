#include <QuantPDE/Core>

////////////////////////////////////////////////////////////////////////////////

#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/toplev.h> // for do_octave_atexit

////////////////////////////////////////////////////////////////////////////////

#include <chrono>   // duration
#include <iostream> // cerr
#include <limits>   // numeric_limits
#include <memory>   // unique_ptr
#include <string>   // string

////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<QuantPDE::DomainBase> parse_grid(int d, int refinement,
		const std::string &s) {
	octave_value_list list = get_top_level_value( s, false );
	Cell cell = list(0).cell_value();

	if(cell.length() != d) {
		throw 1; // FIXME
	}

	if(d == 1) {
		Matrix a0 = cell(0).matrix_value();

		QuantPDE::RectilinearGrid1 initial_grid(
			QuantPDE::Axis(
				a0.data(),
				a0.length()
			)
		);

		return std::unique_ptr<QuantPDE::DomainBase>(
			(QuantPDE::DomainBase*)
			new QuantPDE::RectilinearGrid1(
				initial_grid.refined( refinement )
			)
		);
	} // TODO: Other dimensions

	throw 1; // FIXME
}

template <int Dimension, int ControlDimension = 1>
struct ControlledOperator final : public QuantPDE
		::RawControlledLinearSystem<Dimension, ControlDimension> {

	const QuantPDE::RectilinearGrid<Dimension> &spatial_grid;
	const bool is_controlled;
	int offsets[Dimension];

	template <typename G>
	ControlledOperator(
		G &spatial_grid,
		bool is_controlled
	) noexcept :
		spatial_grid( spatial_grid ),
		is_controlled( is_controlled )
	{
		// Space between tickets
		offsets[0] = 1;
		for(int d = 1; d < Dimension; ++d) {
			offsets[d] = offsets[d-1] * spatial_grid[d].size();
		}
	}

	virtual QuantPDE::Matrix A(double time) {
		QuantPDE::Matrix A = spatial_grid.matrix();
		A.reserve(
			QuantPDE::IntegerVector::Constant(
				spatial_grid.size(),
				1 + 2 * Dimension
			)
		); // Reserve nonzero entries per row

		// Control as a vector
		const QuantPDE::Vector &q = is_controlled ? this->control(0)
				: spatial_grid.zero();

		// Iterate through points on grid
		int i[Dimension];
		for(int row = 0; row < spatial_grid.size(); ++row) {
			double total = 0.;

			// Get coordinates of point
			octave_value_list args;
			args(0) = time; // Time
			for(int d = 0; d < Dimension; ++d) {
				i[d] = (row / offsets[d])
						% spatial_grid[d].size();
				args(1+d) = spatial_grid[d][i[d]];
			}
			args(1+Dimension) = q(row); // Control

			// Get volatility
			Matrix volatility = (
				feval("volatility", args, 1)
			)(0).matrix_value();

			// TODO: Throw error on size mismatch
			// TODO: Throw error on cross-derivative terms

			// Get drifts
			Matrix controlled, uncontrolled;

			{
				auto tmp = feval("drift", args, 1);
				controlled = tmp(0).matrix_value();
				uncontrolled = tmp(1).matrix_value();
			}

			// TODO: Throw error on size mismatch

			for(int d = 0; d < Dimension; ++d) {
				// Skip boundary points
				if(
					   i[d] == 0
					|| i[d] == spatial_grid[d].size() - 1
				) { continue; }

				const QuantPDE::Axis &x = spatial_grid[d];

				const double
					dxb = x[ i[d]     ] - x[ i[d] - 1 ],
					dxc = x[ i[d] + 1 ] - x[ i[d] - 1 ],
					dxf = x[ i[d] + 1 ] - x[ i[d]     ]
				;

				// Local volatility
				const double v = volatility(d, d);

				// Local drift
				const double mu =
					(is_controlled ? controlled(d) : 0.)
					+ uncontrolled(d)
				;

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

			const double rho = (
				feval("discount", args, 1)
			)(0).scalar_value();

			A.insert(row, row) = total + rho;
		}

		A.makeCompressed();
		return A;
	}

	virtual QuantPDE::Vector b(double time) {
		QuantPDE::Vector b = spatial_grid.vector();

		// Control as a vector
		const QuantPDE::Vector &q = is_controlled ? this->control(0)
				: spatial_grid.zero();

		// Iterate through points on grid
		int i[Dimension];
		for(int row = 0; row < spatial_grid.size(); ++row) {
			// Get coordinates of point
			octave_value_list args;
			args(0) = time; // Time
			for(int d = 0; d < Dimension; ++d) {
				i[d] = (row / offsets[d])
						% spatial_grid[d].size();
				args(1+d) = spatial_grid[d][i[d]];
			}
			args(1+Dimension) = q(row); // Control

			// Get flows
			double controlled, uncontrolled;
			{
				auto tmp = feval("continuous_flow", args, 1);
				controlled = tmp(0).scalar_value();
				uncontrolled = tmp(1).scalar_value();
			}

			b(row) =
				(is_controlled ? controlled : 0.)
				+ uncontrolled
			;
		}

		return b;
	}

	virtual bool isATheSame() const {
		// TODO: Check
		return false;
	}

};

void solve(int refinement) {

	////////////////////////////////////////////////////////////////////////
	// Parse constant parameters
	////////////////////////////////////////////////////////////////////////

	// Read spatial dimension
	const int d = get_top_level_value("spatial_dimension").int_value();
	// Error is thrown if dimension is not supported

	// Read grids
	auto spatial_grid = parse_grid(
		d,
		refinement,
		"spatial_grid"
	);
	auto stochastic_control_grid = parse_grid(
		1,
		refinement,
		"stochastic_control_grid"
	);
	auto impulse_control_grid = parse_grid(
		1,
		refinement,
		"impulse_control_grid"
	);

	const bool implicit = get_top_level_value("implicit").int_value();
	const double expiry = get_top_level_value("expiry").scalar_value();

	const int timesteps = get_top_level_value("timesteps").int_value()
			* refinement;

	////////////////////////////////////////////////////////////////////////

	const bool finite_horizon =
			expiry < std::numeric_limits<double>::infinity();

	QuantPDE::ToleranceIteration tolerance_iteration;

	std::unique_ptr<QuantPDE::ControlledLinearSystemBase>
		controlled_operator,
		impulse
	;

	std::unique_ptr<QuantPDE::IterationNode>
		stochastic_policy,
		impulse_policy,
		discretization
	;

	if(d == 1) {
		// Operator
		controlled_operator = std::unique_ptr<
				QuantPDE::ControlledLinearSystemBase>(
			(QuantPDE::ControlledLinearSystemBase*)
			new ControlledOperator<1>(
				*(
					(QuantPDE::RectilinearGrid1*)
					spatial_grid.get()
				),
				implicit // FIXME: method != SEMILAGRANGIAN
			)
		);

		// Policy iteration for stochastic control
		stochastic_policy = std::unique_ptr<QuantPDE::IterationNode>(
			(QuantPDE::IterationNode*)
			new QuantPDE::MinPolicyIteration1_1(
				*(
					(QuantPDE::RectilinearGrid1*)
					spatial_grid.get()
				),
				*(
					(QuantPDE::RectilinearGrid1*)
					stochastic_control_grid.get()
				),
				*(
					(QuantPDE::RawControlledLinearSystem1_1*)
					controlled_operator.get()
				)
			)
		);

		// Policy iteration for stochastic control
		impulse = std::unique_ptr<QuantPDE::ControlledLinearSystemBase>(
			(QuantPDE::ControlledLinearSystemBase*)
			new QuantPDE::Impulse1_1(
				*(
					(QuantPDE::RectilinearGrid1*)
					spatial_grid.get()
				),
				[] (double t, double x, double zeta) {
					octave_value_list args;
					args(0) = t;
					args(1) = x;
					args(2) = zeta;
					const double res = (
						feval("impulse_flow", args, 1)
					)(0).scalar_value();
					return res;
				},
				[] (double t, double x, double zeta) {
					octave_value_list args;
					args(0) = t;
					args(1) = x;
					args(2) = zeta;
					Matrix transition = (
						feval("transition", args, 1)
					)(0).matrix_value();
					return transition(0);
				}
			)
		);
		impulse_policy = std::unique_ptr<QuantPDE::IterationNode>(
			(QuantPDE::IterationNode*)
			new QuantPDE::MinPolicyIteration1_1(
				*(
					(QuantPDE::RectilinearGrid1*)
					spatial_grid.get()
				),
				*(
					(QuantPDE::RectilinearGrid1*)
					impulse_control_grid.get()
				),
				*(
					(QuantPDE::RawControlledLinearSystem1_1*)
					impulse.get()
				)
			)
		);
	} // TODO: Other dimensions

	stochastic_policy->setIteration(tolerance_iteration);
	impulse_policy->setIteration(tolerance_iteration);

	std::unique_ptr<QuantPDE::ReverseTimeIteration> stepper;
	if(finite_horizon) {
		stepper = std::unique_ptr<QuantPDE::ReverseTimeIteration>(
			(QuantPDE::ReverseTimeIteration*)
			new QuantPDE::ReverseConstantStepper(
				0.,                // Initial time
				expiry,            // Expiry time
				expiry / timesteps // Timestep size
			)
		);

		if(implicit) { // FIXME: method != EXPLICIT
			stepper->setInnerIteration(tolerance_iteration);
		}
	}

	QuantPDE::LinearSystem *discretize;
	if(implicit) { // FIXME: method != SEMILAGRANGIAN
		discretize = stochastic_policy.get();
	} else {
		discretize = controlled_operator.get();
	}

	// Finite or infinite horizon
	QuantPDE::Iteration *iteration;
	QuantPDE::IterationNode *penalized;
	if(finite_horizon) {
		if(d == 1) {
			discretization = std::unique_ptr<QuantPDE::IterationNode>(
				(QuantPDE::IterationNode*)
				new QuantPDE::ReverseBDFOne1(
					*(
						(QuantPDE::RectilinearGrid1*)
						spatial_grid.get()
					),
					*discretize
				)
			);
		} // TODO: Other dimensions
		discretization->setIteration(*stepper);

		iteration = stepper.get();
		penalized = discretization.get();
	} else {
		iteration = &tolerance_iteration;
		penalized = stochastic_policy.get();
	}

	// Penalty method
	QuantPDE::PenaltyMethod penalty(
		*(spatial_grid.get()),
		*penalized,
		*(impulse_policy.get())
	);
	penalty.setIteration(tolerance_iteration);

	// Pick root
	QuantPDE::IterationNode *root;
	if(implicit) { // method == EXPLICIT_IMPULSE
		root = &penalty;
	} else {
		root = penalized;
	}

	// TODO: Events

	// Solve
	auto start = std::chrono::steady_clock::now();
	QuantPDE::BiCGSTABSolver solver;
	if(d == 1) {
		auto u = iteration->solve(
			*( (QuantPDE::RectilinearGrid1*) spatial_grid.get() ),
			[=] (double x) {
				octave_value_list args;
				args(0) = expiry; // Payoff
				args(1) = x;
				return (
					feval("exit_flow", args, 1)
				)(0).scalar_value();
			},
			*root,
			solver
		);

		// Print solution
		/*
		std::cout << accessor(
			*(
				(QuantPDE::RectilinearGrid1*)
				spatial_grid.get()
			),
			u
		) << std::endl;
		*/
	}
	auto end = std::chrono::steady_clock::now();
	auto diff = end - start;

	/*
	std::cout << "nodes: " << spatial_grid->size() << std::endl;
	std::cout << "q nodes: " << stochastic_control_grid->size() << std::endl;
	std::cout << "zeta nodes: " << impulse_control_grid->size() << std::endl;

	std::cout << "iterations: " << tolerance_iteration.iterations().back() << std::endl;
	*/

	std::cout
		<< "seconds: "
		<< std::chrono::duration <double> (diff).count()
		<< std::endl
	;
}

int main(int argc, char **argv) {
	// Check number of arguments
	if(argc != 2) {
		std::cerr << "usage: hjbqvi OCTAVE_SOURCE_FILE" << std::endl;
		return 1;
	}

	// Start octave
	const char *argvv[] = { "", "--silent" };
	octave_main(2, (char **) argvv, true /* embedded */);

	// Load source file
	source_file(argv[1]);
	// TODO: Check status and die gracefully if unable to load

	// TODO: Change this to allow for different program options
	solve(0);
	solve(1);
	solve(2);
	solve(3);
	solve(4);

	// Exit octave
	clean_up_and_exit(0);

	return 0;
}

