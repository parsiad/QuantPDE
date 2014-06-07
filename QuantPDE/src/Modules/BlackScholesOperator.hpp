#ifndef QUANT_PDE_MODULES_BLACK_SCHOLES_OPERATOR
#define QUANT_PDE_MODULES_BLACK_SCHOLES_OPERATOR

namespace QuantPDE {

namespace Modules {

class BlackScholesOperator : public LinearOperator {

	const RectilinearGrid1 *G;
	Function2 r, v, q;

public:

	template <typename G, typename F1, typename F2, typename F3>
	BlackScholesOperator(G &grid, F1 &&interest, F2 &&volatility,
			F3 &&dividends) noexcept : G(&grid),
			r( std::forward<F1>(interest) ),
			v( std::forward<F2>(volatility) ),
			q( std::forward<F3>(dividends) ) {
	}

	virtual Matrix discretize(Real time) const {
		auto M_G = G->builder( IntegerVector::Constant(G->size(), 3) );

		const Axis &S = (*G)[0];

		// Interior points
		// alpha_i dt V_{i-1}^{n+1} + (1 + (alpha_i + beta_i + r) dt)
		// 		V_i^{n+1} + beta_i dt V_{i+1}^{n+1} = V_i^n

		for(Index i = 1; i < S.size() - 1; i++) {

			const double r_i = r( time, S[i] );
			const double v_i = v( time, S[i] );
			const double q_i = q( time, S[i] );

			const double
				dSb = S[i]     - S[i - 1],
				dSc = S[i + 1] - S[i - 1],
				dSf = S[i + 1] - S[i]
			;

			const double alpha_common = v_i * v_i * S[i] * S[i]
					/ dSb / dSc;
			const double  beta_common = v_i * v_i * S[i] * S[i]
					/ dSf / dSc;

			// Central
			double alpha_i = alpha_common - (r_i - q_i) * S[i]
					/ dSc;
			double beta_i  =  beta_common + (r_i - q_i) * S[i]
					/ dSc;
			if(alpha_i < 0) {
				// Forward
				alpha_i = alpha_common;
				beta_i  =  beta_common + (r_i - q_i) * S[i]
						/ dSf;
			} else if(beta_i < 0) {
				// Backward
				alpha_i = alpha_common - (r_i - q_i) * S[i]
						/ dSb;
				beta_i  =  beta_common;
			}

			M_G(i, i - 1) = -alpha_i;
			M_G(i, i)     =  alpha_i + beta_i + r_i;
			M_G(i, i + 1) = -beta_i;

		}

		// Boundaries
		// Left:  (1 + r dt) V_i^{n+1} = V_i^n
		// Right:            V_i^{n+1} = V_i^n (from linearity)

		M_G(0, 0) = r( time, S[0] );

		return M_G.matrix();
	}

};

} // Modules

} // QuantPDE

#endif

