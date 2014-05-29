#ifndef QUANT_PDE_MODULES_EUROPEAN_OPTION
#define QUANT_PDE_MODULES_EUROPEAN_OPTION

#include <memory> // std::unique_ptr

namespace QuantPDE {

namespace Modules {

class EuropeanOption final : public Problem1 {

public:

	template <typename F1, typename F2, typename F3, typename F4>
	EuropeanOption(F1 &&payoff, F2 &&interest, F3 &&volatility,
			F4 &&dividends, const Real &start, const Real &end)
			noexcept : Problem1(std::forward<F1>(payoff)) {

		// dV/dt + d^2V/dS^2 sigma^2 S^2 / 2 + dV/dS rS - rV = 0
		this->add(
			std::unique_ptr<Constraint>(
				new BlackScholesEquation(
					std::forward<F2>(interest),
					std::forward<F3>(volatility),
					std::forward<F4>(dividends)
				)
			),
			start,
			end
		);

	}

};

} // Modules

} // QuantPDE

#endif

