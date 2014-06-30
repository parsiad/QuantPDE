#ifndef QUANT_PDE_MODULES_BLACK_SCHOLES
#define QUANT_PDE_MODULES_BLACK_SCHOLES

namespace QuantPDE {

namespace Modules {

/**
Represents the operator $L$ in
\f$V_{t}+LV\equiv V_{t}+\frac{1}{2}\sigma^{2}S^{2}V_{SS}+\left(r-q\right)SV_{S}-rV\f$
\f$r\f$, \f$\sigma\f$, and \f$q\f$ are the usual interest rate, volatility
and (continuous) rate of dividends.

The boundary conditions imposed are \f$V_t - rV = 0\f$ at \f$S=0\f$ and
\f$V_t - qV = 0\f$ at \f$S\rightarrow \infty\f$. The latter is derived by
assuming that the option is linear in the asset for large enough values of the
asset.
**/
class BlackScholes : public ControlledLinearSystem1 {

protected:

	const RectilinearGrid1 &G;
	MultiFunction1 r, v, q, l;

	/**
	 * Constructor for jump-diffusion process. Jumps occur according to a
	 * Poisson process.
	 * @param grid The underlying spatial grid.
	 * @param interest The risk-free interest rate.
	 * @param volatility The volatility of the underlying asset.
	 * @param dividends The continuous dividend rate.
	 * @param meanArrivalTime The mean arrival time of the Poisson process.
	 */
	template <typename G1, typename F1, typename F2, typename F3,
			typename F4>
	BlackScholes(
		G1 &grid,
		F1 &&interest,
		F2 &&volatility,
		F3 &&dividends,
		F4 &&meanArrivalTime
	) noexcept :
		G( grid ),
		r( std::forward<F1>(interest) ),
		v( std::forward<F2>(volatility) ),
		q( std::forward<F3>(dividends) ),
		l( std::forward<F4>(meanArrivalTime) ) {

		registerControl(r);
		registerControl(v);
		registerControl(q);
		registerControl(l);

	}

public:

	/**
	 * Constructor.
	 * @param grid The underlying spatial grid.
	 * @param interest The risk-free interest rate.
	 * @param volatility The volatility of the underlying asset.
	 * @param dividends The continuous dividend rate.
	 */
	template <typename G1, typename F1, typename F2, typename F3>
	BlackScholes(
		G1 &grid,
		F1 &&interest,
		F2 &&volatility,
		F3 &&dividends
	) noexcept :
		G( grid ),
		r( std::forward<F1>(interest) ),
		v( std::forward<F2>(volatility) ),
		q( std::forward<F3>(dividends) ),
		l( 0. ) {

		registerControl(r);
		registerControl(v);
		registerControl(q);

	}

	virtual Matrix A(Real t) {
		auto M_G = G.builder( IntegerVector::Constant(G.size(), 3) );

		const Axis &S = G[0];
		const Index n = S.size();

		// Interior points
		// alpha_i dt V_{i-1}^{n+1} + (1 + (alpha_i + beta_i + r) dt)
		// 		V_i^{n+1} + beta_i dt V_{i+1}^{n+1} = V_i^n

		for(Index i = 1; i < n - 1; i++) {

			const Real r_i = r( t, S[i] );
			const Real v_i = v( t, S[i] );
			const Real q_i = q( t, S[i] );
			const Real l_i = l( t, S[i] );

			const Real kappa = 0.;

			const Real
				dSb = S[i]     - S[i - 1],
				dSc = S[i + 1] - S[i - 1],
				dSf = S[i + 1] - S[i]
			;

			const Real alpha_common = v_i * v_i * S[i] * S[i]
					/ dSb / dSc;
			const Real  beta_common = v_i * v_i * S[i] * S[i]
					/ dSf / dSc;

			// Central
			Real alpha_i = alpha_common - (r_i - q_i - l_i * kappa)
					* S[i] / dSc;
			Real beta_i  =  beta_common + (r_i - q_i - l_i * kappa)
					* S[i] / dSc;
			if(alpha_i < 0) {
				// Forward
				alpha_i = alpha_common;
				beta_i  =  beta_common + (r_i - q_i
						- l_i * kappa) * S[i] / dSf;
			} else if(beta_i < 0) {
				// Backward
				alpha_i = alpha_common - (r_i - q_i
						- l_i * kappa) * S[i] / dSb;
				beta_i  =  beta_common;
			}

			M_G(i, i - 1) = -alpha_i;
			M_G(i, i)     = alpha_i + beta_i + r_i + l_i;
			M_G(i, i + 1) = -beta_i;

		}

		// Boundaries
		// Left:  (1 + r dt) V_i^{n+1} = V_i^n
		// Right:              V(t, S) = g(t) S (linearity assumption)

		M_G(0, 0) = r(t, S[0]);
		M_G(n-1, n-1) = q(t, S[n-1]);

		return M_G.matrix();
	}

	virtual Vector b(Real) {
		return G.zero();
	}

	virtual bool isATheSame() const {
		return r.isConstantInTime() && v.isConstantInTime()
				&& q.isConstantInTime() && l.isConstantInTime();
	}

};

/**
Represents the operator $L$ in
\f$V_{t}+LV\equiv V_{t}+\frac{1}{2}\sigma^{2}S^{2}V_{SS}+\left(r-q-\lambda\kappa\right)SV_{S}-\left(r+\lambda\right)V+\lambda\int_{0}^{\infty}V\left(t,S\eta\right)g\left(\eta\right)d\eta.\f$
\f$r\f$, \f$\sigma\f$, and \f$q\f$ are the usual interest rate, volatility
and (continuous) rate of dividends. \f$\lambda\f$ is the mean arrival
time of the Poisson process responsible for generating the jumps.
Assuming a jump has occured, \f$g\left(\eta\right)\f$ is the probability
density of a jump of amplitude \f$\eta\f$, with
\f$\kappa\equiv E\left[\eta\right]-1\f$ describing the expected relative change
in the stock.

The boundary conditions imposed are \f$V_t - rV = 0\f$ at \f$S=0\f$ and
\f$V_t - qV = 0\f$ at \f$S\rightarrow \infty\f$. The latter is derived by
assuming that the option is linear in the asset for large enough values of the
asset.
**/
class BlackScholesJumpDiffusion final : public BlackScholes {

	Function1 g;

public:

	/**
	 * Constructor for jump-diffusion process. Jumps occur according to a
	 * Poisson process.
	 * @param grid The underlying spatial grid.
	 * @param interest The risk-free interest rate.
	 * @param volatility The volatility of the underlying asset.
	 * @param dividends The continuous dividend rate.
	 * @param meanArrivalTime The mean arrival time of the Poisson process.
	 * @param density The probability density of the jump amplitude.
	 */
	template <typename G, typename F1, typename F2, typename F3,
			typename F4, typename F5>
	BlackScholesJumpDiffusion(
		G &grid,
		F1 &&interest,
		F2 &&volatility,
		F3 &&dividends,
		F4 &&meanArrivalTime,
		F5 &&density
	) noexcept : BlackScholes(
		grid,
		std::forward<F1>(interest),
		std::forward<F2>(volatility),
		std::forward<F3>(dividends),
		std::forward<F4>(meanArrivalTime)
	), g(std::forward<F5>(density)) {
	}

	virtual Vector b(Real) {
		//const Vector &v = this->iterand(0);

		// TODO: Explicit jump-diffusion
		return G.zero();
	}

};

} // Modules

} // QuantPDE

#endif

