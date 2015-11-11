#ifndef QUANT_PDE_MODULES_LAMBDAS_DENSITIES_HPP
#define QUANT_PDE_MODULES_LAMBDAS_DENSITIES_HPP

#include <cmath> // std::exp, std::log, std::sqrt, std::pow

namespace QuantPDE {

namespace Modules {

/**
 * Returns a lambda function for the probability density
 * \f$
 * 	f(J) =
 * 		\frac{1}{J\sigma\sqrt{2\pi}}
 * 		\exp\left(\frac{\left(-\ln J-\mu\right)^2}{2\sigma^2}\right)
 * \f$
 * @param mu
 * @param sigma
 * @return A lambda function.
 */
inline Function1 lognormal(Real mu = 0., Real sigma = 1.) {
	assert(sigma > 0.);

	return [mu, sigma] (Real x) -> Real {
		assert(x > 0.);

		if(x == 0.) {
			return 0;
		}

		return 1. / (x * sigma * std::sqrt(2. * M_PI))
				* std::exp( -((std::log(x)-mu)*(std::log(x)-mu))
				/ (2. * sigma * sigma) );
	};
}

/**
 * Returns a lambda function for the probability density
 * \f$
 * 	f(J) =
 * 	       p  \eta_1 J^{-\eta_1 - 1} 1_{\{ J \geq 1 \} }
 * 	+ (1 - p) \eta_2 J^{ \eta_2 - 1} 1_{\{ J <    1 \} }
 * $\f
 * @param p Probability of upward jump.
 * @param eta_1 The upward jump random variable has mean 1/eta_1.
 * @param eta_2 The downward jump random variable has mean 1/eta_2.
 * @return A lambda function.
 */
inline Function1 doubleExponential(Real p, Real eta_1, Real eta_2) {
	assert(p >= 0.);
	assert(eta_1 > 1.); // Ensures finite expectation
	assert(eta_2 > 0.);

	return [p, eta_1, eta_2] (Real x) -> Real {
		assert(x > 0.);

		if(x == 0.) {
			return 0.;
		}

		return (x >= 1)
			? (   p  * eta_1 * std::pow(x, -eta_1 - 1))
			: ((1-p) * eta_2 * std::pow(x,  eta_2 - 1));
	};
}

} // Modules

} // QuantPDE

#endif

