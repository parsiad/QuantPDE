#ifndef QUANT_PDE_MODULES_CONSTRAINTS_BLACK_SCHOLES_EQUATION
#define QUANT_PDE_MODULES_CONSTRAINTS_BLACK_SCHOLES_EQUATION

namespace QuantPDE {

namespace Modules {

class BlackScholesEquation final : public Constraint {

public:

	const Function2 interest, volatility, dividends;

	template <typename F1, typename F2, typename F3>
	BlackScholesEquation(F1 &&interest, F2 &&volatility, F3 &&dividends)
			noexcept : interest( std::forward<F1>(interest) ),
			volatility( std::forward<F2>(volatility) ),
			dividends( std::forward<F3>(dividends) ) {
	}

	virtual std::string identifier() const {
		return QUANT_PDE_IDENTIFIER(BlackScholesEquation);
	}

};

} // Modules

} // QuantPDE

#endif

