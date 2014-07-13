#ifndef QUANT_PDE_MODULES_LAMBDAS_DENSITIES_HPP
#define QUANT_PDE_MODULES_LAMBDAS_DENSITIES_HPP

#include <cmath> // std::exp, std::log, std::sqrt

namespace QuantPDE {

namespace Modules {

/**
 * @param mu The mean.
 * @param sigma The standard deviation.
 * @return Lambda function for a lognormal probability density.
 */
Function1 lognormal(Real mu = 0., Real sigma = 1.) {
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

} // Modules

} // QuantPDE

#endif

