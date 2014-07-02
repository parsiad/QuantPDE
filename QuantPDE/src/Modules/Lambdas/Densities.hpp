#ifndef QUANT_PDE_MODULES_LAMBDAS_DENSITIES
#define QUANT_PDE_MODULES_LAMBDAS_DENSITIES

#include <cmath> // std::exp, std::log, std::sqrt

namespace QuantPDE {

namespace Modules {

/**
 * @param mu The mean.
 * @param sigma The standard deviation.
 * @return Lambda function for a lognormal density.
 */
Function1 lognormalDensity(Real mu = 0., Real sigma = 1.) {
	return [mu, sigma] (Real x) {
		return 1. / (x * sigma * std::sqrt(2. * M_PI))
				* std::exp( -((std::log(x)-mu)*(std::log(x)-mu))
				/ (2. * sigma * sigma) );
	};
}

} // Modules

} // QuantPDE

#endif

