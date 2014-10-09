#ifndef QUANT_PDE_CORE_IMPULSE
#define QUANT_PDE_CORE_IMPULSE

namespace QuantPDE {

template <Index Dimension, Index ControlDimension>
class Impulse final : public ControlledLinearSystem<Dimension> {

	const RectilinearGrid2 &grid;

	// t, x1, ..., xn, q1, ..., qn
	Function<1 + Dimension + ControlDimension> transition, cashflow;

	// TODO: Change name of Coefficient
	std::array<Coefficient<Dimension>, ControlDimension> controls;

public:

	template <typename G, typename F1, typename F2, typename ...C>
	Impulse(G &grid, F1 &&transition, F2 &&cashflow, C &&...controls)
			noexcept : grid(grid), transition(transition),
			cashflow(cashflow), controls(controls...) {
	}

	virtual Matrix A(Real t) {
		Matrix A;

		// TODO

		A.makeCompressed();
		return A;
	}

	virtual Vector b(Real t) {
		Vector b;

		// TODO

		return b;
	}

};

}

#endif

