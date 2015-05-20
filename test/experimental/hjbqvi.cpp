//#include <QuantPDE/Core>

////////////////////////////////////////////////////////////////////////////////

#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/toplev.h> // for do_octave_atexit

////////////////////////////////////////////////////////////////////////////////

#include <stdlib.h> // mkstemp
#include <string.h> // strdup
#include <iostream> // cerr

////////////////////////////////////////////////////////////////////////////////

#if 0
/**
 * Specifies the HJBQVI problem for u(t,x) completely.
 * @tparam Dimension The dimension of the independent spatial variable x.
 * @tparam StochasticControlDimension The dimension of the stochastic control.
 * @tparam ImpulseControlDimension The dimension of the impulse control.
 */
template <
	Index Dimension,
	Index StochasticControlDimension,
	Index ImpulseControlDimension
>
struct ProblemSpecification {

	/**
	 * If true, the stochastic control grid is refined along with the
	 * spatial grid.
	 */
	const bool refineStochasticControl;

	/**
	 * If true, the impulse control grid is refined along with the spatial
	 * grid.
	 */
	const bool refineImpulseControl;

	/**
	 * The spatial grid.
	 */
	const RectilinearGrid<Dimension> grid;

	/**
	 * The stochastic control grid.
	 */
	const RectilinearGrid<Dimension> stochasticControlGrid;

	/**
	 * The impulse control grid.
	 */
	const RectilinearGrid<Dimension> impulseControlGrid;

	/**
	 * The local volatility.
	 */
	const Function<1+Dimension> volatility;

	/**
	 * The controlled component of the local drift.
	 */
	const Function<1+Dimension+StochasticControlDimension> controlledDrift;

	/**
	 * The uncontrolled component of the local drift.
	 */
	const Function<1+Dimension> uncontrolledDrift;

	/**
	 * The controlled component of the local discount.
	 */
	const Function<1+Dimension+StochasticControlDimension>
			controlledDiscount;

	/**
	 * The uncontrolled component of the local discount.
	 */
	const Function<1+Dimension> uncontrolledDiscount;

	/**
	 * The intervention state transition function.
	 */
	const Function<1+Dimension+ImpulseControlDimension> stateTransition;

	/**
	 * The intervention flow (e.g. cash flow) function.
	 */
	const Function<1+Dimension+ImpulseControlDimension> interventionFlow;

}
#endif

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

	// Get dimension
	// octave_value_list a = get_top_level_value("dimension", false);

	octave_value_list functionArguments;
	functionArguments(0) = 0.;
	functionArguments(1) = 10.;
	functionArguments(2) = 20.;
	const octave_value_list result = feval("transition", functionArguments, 1);


	// Exit octave
	clean_up_and_exit(0);

	return 0;
}
