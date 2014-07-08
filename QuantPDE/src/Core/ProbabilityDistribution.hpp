#ifndef QUANT_PDE_CORE_PROBABILITY_DISTRIBUTION_HPP
#define QUANT_PDE_CORE_PROBABILITY_DISTRIBUTION_HPP

#include <array>

namespace QuantPDE {

/**
 * A probability distribution.
 */
template <Index Dimension>
class ProbabilityDistribution {

	/*
	Function<Dimension + 1> _cdf() const {
		Function<Dimension + 1> pdf = pdf();
		return [ QUANT_PDE_MOVE_CAPTURE(Function<Dimension + 1>, pdf) ]
				(Real t, Ts ...coordinates) {

		}
	}
	*/

public:

	/**
	 * Returns the support of the associated probability density function.
	 *
	 * Only probability densities defined on cells are supported (i.e. the
	 * support must be a Cartesian product of intervals).
	 *
	 * The first Dimension indices are the lower bounds of integration. The
	 * next Dimension indices are the upper bounds of integration.
	 *
	 * @return The support of the associated probability density function.
	 */
	virtual std::array<Real, Dimension * 2> support() const = 0;

	/**
	 * Returns the associated probability density function.
	 * @return A lambda function whose first argument is time.
	 */
	virtual Function<Dimension + 1> pdf() const = 0;

	/**
	 * Returns the associated cumulative distribution function.
	 * @return A lambda function whose first argument is time.
	 */
	virtual Function<Dimension + 1> cdf() const {
		// TODO
	}


	virtual Real mean() const {
		auto sup = support();

		// TODO: Unroll
		for(Index i = 0; i < Dimension; i++) {

		}
	}

	virtual Real variance() const {
		// TODO: Integrate
	}

	virtual Real median() const {
		// TODO: Bisection to find where CDF = 1/2
	}



};

} // QuantPDE

#endif

