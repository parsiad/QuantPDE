#ifndef QUANT_PDE_MODULES_PAYOFFS
#define QUANT_PDE_MODULES_PAYOFFS

#include <functional> // std::bind, std::function

namespace QuantPDE {

namespace Modules {

/**
 * @param stock The stock price.
 * @param strike The strike price.
 * @return Payoff for a vanilla call option, \f$\max\left(S - K, 0\right)\f$.
 */
Real call(Real stock, Real strike) {
	return stock > strike ? stock - strike : 0.;
}

/**
 * @param stock The stock price.
 * @param strike The strike price.
 * @return Payoff for a vanilla put option, \f$\max\left(S - K, 0\right)\f$.
 */
Real put(Real stock, Real strike) {
	return call(strike, stock);
}

#define QUANT_PDE_MODULES_PAYOFFS_CALL_FIXED_STRIKE(STRIKE) std::bind( \
		QuantPDE::Modules::Payoffs::call, std::placeholders::_1, STRIKE)

#define QUANT_PDE_MODULES_PAYOFFS_PUT_FIXED_STRIKE(STRIKE) std::bind( \
		QuantPDE::Modules::Payoffs::put, std::placeholders::_1, STRIKE)

} // Modules

}

#endif
