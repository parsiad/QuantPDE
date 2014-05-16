#ifndef QUANT_PDE_MODULES_PAYOFFS
#define QUANT_PDE_MODULES_PAYOFFS

#include <functional> // std::bind, std::function

namespace QuantPDE {

namespace Modules {

namespace Payoffs {

/**
 * The payoff for a vanilla call option, \f$\max\left(S - K, 0\right)\f$.
 */
const std::function<Real (Real, Real)> call = [] (Real S, Real K) {
	return S > K ? S - K : 0.;
};

/**
 * The payoff for a vanilla put option, \f$\max\left(S - K, 0\right)\f$.
 */
const std::function<Real (Real, Real)> put = [] (Real S, Real K) {
	return K < S ? K - S : 0.;
};

#define QUANT_PDE_MODULES_PAYOFFS_CALL_FIXED_STRIKE(strike) std::bind( \
		QuantPDE::Modules::Payoffs::call, std::placeholders::_1, strike)

#define QUANT_PDE_MODULES_PAYOFFS_PUT_FIXED_STRIKE(strike) std::bind( \
		QuantPDE::Modules::Payoffs::put, std::placeholders::_1, strike)

} // Payoffs

} // Modules

}

#endif
