#ifndef QUANT_PDE_TEST_GMWB_HPP
#define QUANT_PDE_TEST_GMWB_HPP

#include <algorithm> // max
#include <tuple>     // get

namespace QuantPDE {

class ImpulseWithdrawal final : public ControlledLinearSystem2,
		public IterationNode {

	RectilinearGrid2 &grid;
	Noncontrollable2 kappa;

	Controllable2 control;

public:

	template <typename G, typename F1>
	ImpulseWithdrawal(G &grid, F1 &&kappa) noexcept
			: grid(grid), kappa(kappa), control( Control2(grid) ) {
		registerControl( control );
	}

	virtual Matrix A(Real t) {
		Matrix M = grid.matrix();
		M.reserve(IntegerVector::Constant(grid.size(), 4));

		Index k = 0;
		for(auto node : grid) {
			const Real S = node[0]; // Investment
			const Real W = node[1]; // Withdrawal

			// Amount withdrawn pre-penalty
			const Real gamma = control(t, S, W) * W;

			const Real Splus = std::max(S - gamma, 0.);
			const Real Wplus = W - gamma;

			// Interpolation data
			auto data = interpolationData<2>(grid, {{Splus,Wplus}});

			const Index i0 = std::get<0>( data[0] );
			const Index i1 = std::get<0>( data[1] );
			const Real  w0 = std::get<1>( data[0] );
			const Real  w1 = std::get<1>( data[1] );

			assert( (grid[0][i0+1] - Splus)
					/ (grid[0][i0+1] - grid[0][i0]) == w0 );
			assert( (grid[1][i1+1] - Wplus)
					/ (grid[1][i1+1] - grid[1][i1]) == w1 );

			const Index j = grid.index(i0, i1);

			M.insert(k, j                     ) =    w0  *    w1 ;
			M.insert(k, j     + grid[0].size()) =    w0  * (1-w1);
			M.insert(k, j + 1                 ) = (1-w0) *    w1 ;
			M.insert(k, j + 1 + grid[0].size()) = (1-w0) * (1-w1);

			++k;
		}

		M.makeCompressed();
		return grid.identity() - M;
	}

	virtual Vector b(Real t) {
		Vector b = grid.vector();

		for(auto node : accessor(grid, b)) {
			const Real S = (&node)[0]; // Investment
			const Real W = (&node)[1]; // Withdrawal

			// Amount withdrawn, pre-penalty
			const Real gamma = control(t, S, W) * W;

			// Cashflow minus adjustment
			*node = (1 - kappa(t, S, W)) * gamma - epsilon;
		}

		return b;
	}

};

}

#endif
