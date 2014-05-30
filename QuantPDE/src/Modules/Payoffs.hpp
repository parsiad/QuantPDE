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
Real callPayoff(Real stock, Real strike) {
	return stock > strike ? stock - strike : 0.;
}

/**
 * @param stock The stock price.
 * @param strike The strike price.
 * @return Payoff for a vanilla put option, \f$\max\left(S - K, 0\right)\f$.
 */
Real putPayoff(Real stock, Real strike) {
	return callPayoff(strike, stock);
}

#define QUANT_PDE_MODULES_CALL_PAYOFF_FIXED_STRIKE(STRIKE) std::bind( \
		QuantPDE::Modules::Payoffs::call, std::placeholders::_1, STRIKE)

#define QUANT_PDE_MODULES_PUT_PAYOFF_FIXED_STRIKE(STRIKE) std::bind( \
		QuantPDE::Modules::Payoffs::put, std::placeholders::_1, STRIKE)

} // Modules

}

#endif
