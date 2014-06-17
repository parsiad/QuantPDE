#ifndef QUANT_PDE_MODULES_PAYOFFS_VANILLA
#define QUANT_PDE_MODULES_PAYOFFS_VANILLA

#include <functional> // std::bind, std::function

namespace QuantPDE {

namespace Modules {

namespace Payoffs {

/**
 * @param S The stock price.
 * @param K The strike price.
 * @return Payoff for a vanilla call option, \f$\max\left(S - K, 0\right)\f$.
 */
Real call(Real S, Real K) {
	return S < K ? 0. : S - K;
}

#define QUANT_PDE_MODULES_PAYOFFS_CALL_FIXED_STRIKE(STRIKE) std::bind( \
		QuantPDE::Modules::Payoffs::call, std::placeholders::_1, STRIKE)

/**
 * @param S The stock price.
 * @param K The strike price.
 * @return Payoff for a vanilla put option, \f$\max\left(S - K, 0\right)\f$.
 */
Real put(Real S, Real K) {
	return S < K ? K - S : 0.;
}

#define QUANT_PDE_MODULES_PAYOFFS_PUT_FIXED_STRIKE(STRIKE) std::bind( \
		QuantPDE::Modules::Payoffs::put, std::placeholders::_1, STRIKE)

/**
 * @param stock The stock price.
 * @param strike The strike price.
 * @return Payoff for a straddle; the sum of the payoffs of a call and a put.
 */
Real straddle(Real S, Real K) {
	return S < K ? K - S : S - K;
}

#define QUANT_PDE_MODULES_PAYOFFS_STRADDLE_FIXED_STRIKE(STRIKE) std::bind(   \
		QuantPDE::Modules::Payoffs::straddle, std::placeholders::_1, \
		STRIKE)

} // Payoffs

} // Modules

}

#endif
