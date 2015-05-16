#ifndef QUANT_PDE_CORE_PROBABILITY_DISTRIBUTION_HPP
#define QUANT_PDE_CORE_PROBABILITY_DISTRIBUTION_HPP

#include <array>  // std::array

namespace QuantPDE {

/**
 * A probability distribution.
 */
template <Index Dimension>
class ProbabilityDistribution {

	/*
	template <typename F, int ...Indices>
	static inline Real integrate(F &&function, const Real *array,
			Sequence<Indices...>) const {
		typedef TrapezoidalRule<Dimension, 1> T;
		typedef AdaptiveQuadrature<Dimension, T> Integral;

		return Integral(function, array[Indices]...)(
				array[Indices + Dimension]...);
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
	virtual Function<Dimension> pdf() const = 0;

	// TODO: Default behaviour for cdf. mean, variance, and median using
	//       support() and pdf()

	/**
	 * Returns the associated cumulative distribution function.
	 * @return A lambda function whose first argument is time.
	 */
	//virtual Function<Dimension + 1> cdf() const = 0;

	/**
	 * @return The mean.
	 */
	virtual Real mean() const = 0;

	/**
	 * @return The variance.
	 */
	//virtual Real variance() const = 0;

	/**
	 * @return The median.
	 */
	//virtual Real median() const = 0;

};

// TODO: Probability distribution that changes with time

class Lognormal1 : public ProbabilityDistribution<1> {

	Real mu, sigma;

public:

	Lognormal1(Real mu, Real sigma) noexcept : mu(mu), sigma(sigma) {
		assert(sigma > 0.);
	}

	virtual std::array<Real, 2> support() const {
		return {{ 0., std::numeric_limits<Real>::infinity() }};
	}

	virtual Function1 pdf() const {
		return [this] (Real x) -> Real {
			assert(x > 0.);

			if(x == 0.) {
				return 0;
			}

			const Real y = std::log(x) - mu;

			return 1. / (x * sigma * std::sqrt(2. * M_PI))
					* std::exp( -y * y )
					/ (2. * sigma * sigma);
		};
	}

	virtual Real mean() const {
		return std::exp( mu * sigma * sigma / 2. );
	}

};

} // QuantPDE

#endif

