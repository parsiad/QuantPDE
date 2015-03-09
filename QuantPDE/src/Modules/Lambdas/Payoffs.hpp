#ifndef QUANT_PDE_MODULES_LAMBDAS_PAYOFFS_HPP
#define QUANT_PDE_MODULES_LAMBDAS_PAYOFFS_HPP

namespace QuantPDE {

namespace Modules {

/**
 * @param K strike price.
 * @return Lambda function payoff for a digital call option, \f$1_{S \geq K}\f$.
 */
inline Function1 digitalCallPayoff(Real K) {
	return [K] (Real S) { return S < K ? 0. : 1.; };
}

/**
 * @param K strike price.
 * @return Lambda function payoff for a digital put option, \f$1_{S \leq K}\f$.
 */
inline Function1 digitalPutPayoff(Real K) {
	return [K] (Real S) { return S > K ? 0. : 1.; };
}

/**
 * @param K Strike price.
 * @return Lambda function payoff for a vanilla call option,
 *         \f$\left(S\right)\equiv max\left(S - K, 0\right)\f$.
 */
inline Function1 callPayoff(Real K) {
	return [K] (Real S) { return S < K ? 0. : S - K; };
}

/**
 * @param K Strike price.
 * @return Lambda function payoff for a vanilla put option,
 *         \f$\left(S\right)\equiv max\left(S - K, 0\right)\f$.
 */
inline Function1 putPayoff(Real K) {
	return [K] (Real S) { return S < K ? K - S : 0.; };
}

/**
 * @param K Strike price.
 * @return Lambda function payoff for a straddle; the sum of the payoffs of a
 *         call and a put.
 */
inline Function1 straddlePayoff(Real K) {
	return [K] (Real S) { return S < K ? K - S : S - K; };
}

} // Modules

} // QuantPDE

#endif
