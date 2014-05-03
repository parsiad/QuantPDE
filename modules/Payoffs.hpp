#ifndef QUANT_PDE_MODULES_PAYOFFS
#define QUANT_PDE_MODULES_PAYOFFS

namespace QuantPDE {

namespace Modules {

namespace Payoffs {

/**
 * The payoff for a vanilla call option, \f$\max\left(x - K, 0\right)\f$.
 */
class Call : public Function1d {

	double strike;

	virtual double get(double x) const {
		double p = x - strike;
		return p > 0. ? p : 0.;
	}

public:

	/**
	 * Constructor.
	 * @param strike The strike price \f$K\f$.
	 */
	Call(double strike) : strike(strike) {
	}

	/**
	 * Copy constructor.
	 */
	Call(const Call &that) : strike(that.strike) {
	}

};

/**
 * The payoff for a vanilla put option, \f$\max\left(K - x, 0\right)\f$.
 */
class Put : public Function1d {

	double strike;

	virtual double get(double x) const {
		double p = strike - x;
		return p > 0. ? p : 0.;
	}

public:

	/**
	 * Constructor.
	 * @param strike The strike price \f$K\f$.
	 */
	Put(double strike) : strike(strike) {
	}

	/**
	 * Copy constructor.
	 */
	Put(const Put &that) : strike(that.strike) {
	}

};

}

}

}

#endif
