#ifndef QUANT_PDE_MODULES_BLACK_SCHOLES
#define QUANT_PDE_MODULES_BLACK_SCHOLES

#include <complex> // std::complex, std::conj
#include <vector>  // std::vector

namespace QuantPDE {

namespace Modules {

/**
 * Represents the operator $L$ in
 * \f$V_{t}+LV\equiv V_{t}+\frac{1}{2}\sigma^{2}S^{2}V_{SS}+\left(r-q\right)SV_{S}-rV\f$
 * \f$r\f$, \f$\sigma\f$, and \f$q\f$ are the usual interest rate, volatility
 * and (continuous) rate of dividends.
 *
 * The boundary conditions imposed are \f$V_t - rV = 0\f$ at \f$S=0\f$ and
 * \f$V_t - qV = 0\f$ at \f$S\rightarrow \infty\f$. The latter is derived by
 * assuming that the option is linear in the asset for large enough values of
 * the asset.
**/
class BlackScholes : public ControlledLinearSystem1 {

	Coefficient1 r, v, q;
	Real kappa;

	void (BlackScholes::*_computeKappa)(Real);

protected:

	// Integration rule
	typedef AdaptiveQuadrature1<TrapezoidalRule1<>> Integral;

	const RectilinearGrid1 &G;
	Coefficient1 l;
	MultiFunction1 g;

	void pass(Real) {
	}

	inline void computeKappa(Real t) {
		// Computes (E[y]-1) where y is an r.v. with probability density
		// g : [0, Infinity) -> [0, Infinity)
		kappa = Integral(
			[=] (Real y) { return exp(2 * y) * g(t, exp(y)); },
			-std::numeric_limits<Real>::infinity()
		)( std::numeric_limits<Real>::infinity() ) - 1.;
	}

	/**
	 * Constructor for jump-diffusion process. Jumps occur according to a
	 * Poisson process.
	 * @param grid The underlying spatial grid.
	 * @param interest The risk-free interest rate.
	 * @param volatility The volatility of the underlying asset.
	 * @param dividends The continuous dividend rate.
	 * @param meanArrivalTime The mean arrival time of the Poisson process.
	 * @param jumpAmplitudeDensity The jump amplitude probability density.
	 */
	template <typename G1, typename F1, typename F2, typename F3,
			typename F4, typename F5>
	BlackScholes(
		G1 &grid,
		F1 &&interest,
		F2 &&volatility,
		F3 &&dividends,
		F4 &&meanArrivalTime,
		F5 &&jumpAmplitudeDensity
	) noexcept :
		r( std::forward<F1>(interest) ),
		v( std::forward<F2>(volatility) ),
		q( std::forward<F3>(dividends) ),
		G( grid ),
		l( std::forward<F4>(meanArrivalTime) ),
		g( std::forward<F5>(jumpAmplitudeDensity) )
	{
		registerControl(r);
		registerControl(v);
		registerControl(q);
		registerControl(l);
		// g is not controllable

		if(g.isConstantInTime()) {
			computeKappa(-1.);
			_computeKappa = &BlackScholes::pass;
		} else {
			_computeKappa = &BlackScholes::computeKappa;
		}
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
		r( std::forward<F1>(interest) ),
		v( std::forward<F2>(volatility) ),
		q( std::forward<F3>(dividends) ),
		G( grid ),
		l( 0. ),
		g( [] (Real S) { return 0.; } )
	{
		registerControl(r);
		registerControl(v);
		registerControl(q);
	}

	virtual Matrix A(Real t) {
		auto M_G = G.builder( IntegerVector::Constant(G.size(), 3) );

		const Axis &S = G[0];
		const Index n = S.size();

		(this->*_computeKappa)(t);

		// Interior points
		// alpha_i dt V_{i-1}^{n+1} + (1 + (alpha_i + beta_i + r) dt)
		// 		V_i^{n+1} + beta_i dt V_{i+1}^{n+1} = V_i^n

		for(Index i = 1; i < n - 1; i++) {

			const Real r_i = r( t, S[i] );
			const Real v_i = v( t, S[i] );
			const Real q_i = q( t, S[i] );
			const Real l_i = l( t, S[i] );

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
 * Represents the operator $L$ in
 * \f$V_{t}+LV\equiv V_{t}+\frac{1}{2}\sigma^{2}S^{2}V_{SS}+\left(r-q-\lambda\kappa\right)SV_{S}-\left(r+\lambda\right)V+\lambda\int_{0}^{\infty}V\left(t,S\eta\right)g\left(\eta\right)d\eta.\f$
 * \f$r\f$, \f$\sigma\f$, and \f$q\f$ are the usual interest rate, volatility
 * and (continuous) rate of dividends. \f$\lambda\f$ is the mean arrival
 * time of the Poisson process responsible for generating the jumps.
 * Assuming a jump has occured, \f$g\left(\eta\right)\f$ is the probability
 * density of a jump of amplitude \f$\eta\f$, with
 * \f$\kappa\equiv E\left[\eta\right]-1\f$ describing the expected relative
 * change in the stock.
 *
 * The boundary conditions imposed are \f$V_t - rV = 0\f$ at \f$S=0\f$ and
 * \f$V_t - qV = 0\f$ at \f$S\rightarrow \infty\f$. The latter is derived by
 * assuming that the option is linear in the asset for large enough values of
 * the asset.
 *
 * The integral introduced by the jump term is handled using the FFT correlation
 * integral method described in [1].
 *
 * [1] d'Halluin, Yann, Peter A. Forsyth, and Kenneth R. Vetzal. "Robust
 * numerical methods for contingent claims under jump diffusion processes." IMA
 * Journal of Numerical Analysis 25.1 (2005): 87-112.
**/
class BlackScholesJumpDiffusion final : public IterationNode,
		public BlackScholes {

	inline RectilinearGrid1 initializeGrid() {
		// Spatial axis
		const Axis &S = G[0];
		const Index n = S.size();

		// Take the number of points in the frequency domain to be the
		// smallest power of 2 larger than n
		N = 2;
		while(N < n) {
			N *= 2;
		}

		// Find leftmost point
		int left = 0;
		for(Index i = 0; i < n; i++) {
			if( S[i] > 0. ) {
				left = i;
				break;
			}
		}
		assert( left < n - 1 );

		// Min:Step:Max
		x0 = std::log( S[left] );
		const Real xf = std::log( S[n-1] );
		dx = (xf - x0) / (N - 1);

		// Create and return grid
		return RectilinearGrid1( Axis::uniform(x0, xf, N) );
	}

	void computeDensityFFT(Real t) {
		// Transformed density
		auto fbar = [&] (Real x) {
			return g(t, std::exp(x)) * std::exp(x);
		};

		// Compute FFT of transformed density
		std::vector<Real> fprime;
		fprime.reserve(N);
		for(Index i = 0; i <= N/2; i++) {
			// Integrate around x_i
			const Real a = x0 + dx * (-.5 + i);
			const Real b = x0 + dx * ( .5 + i);
			fprime.push_back( Integral(fbar, a)(b) );
		}
		for(Index i = N/2+1; i < N; i++) {
			// Integrate around x_{i - N}
			const Real a = x0 + dx * (-.5 + i - N);
			const Real b = x0 + dx * ( .5 + i - N);
			fprime.push_back( Integral(fbar, a)(b) );
		}
		fprimeFFT.reserve(N);
		fft.fwd(fprimeFFT, fprime);
	}

	Eigen::FFT<Real> fft;

	Index N;
	Real x0, dx;
	RectilinearGrid1 F;

	std::vector<std::complex<Real>> fprimeFFT;

	void (BlackScholesJumpDiffusion::*_computeDensityFFT)(Real);

public:

	/**
	 * Constructor for jump-diffusion process. Jumps occur according to a
	 * Poisson process.
	 * @param grid The underlying spatial grid.
	 * @param interest The risk-free interest rate.
	 * @param volatility The volatility of the underlying asset.
	 * @param dividends The continuous dividend rate.
	 * @param meanArrivalTime The mean arrival time of the Poisson process.
	 * @param jumpAmplitudeDensity The jump amplitude probability density.
	 */
	template <typename G, typename F1, typename F2, typename F3,
			typename F4, typename F5>
	BlackScholesJumpDiffusion(
		G &grid,
		F1 &&interest,
		F2 &&volatility,
		F3 &&dividends,
		F4 &&meanArrivalTime,
		F5 &&jumpAmplitudeDensity
	) noexcept :
		BlackScholes(
			grid,
			std::forward<F1>(interest),
			std::forward<F2>(volatility),
			std::forward<F3>(dividends),
			std::forward<F4>(meanArrivalTime),
			std::forward<F5>(jumpAmplitudeDensity)
		),
		F( initializeGrid() )
	{
		if(g.isConstantInTime()) {
			// Precompute once
			computeDensityFFT(-1.);
			_computeDensityFFT = &BlackScholesJumpDiffusion::pass;
		} else {
			// Compute every time b() is computed
			_computeDensityFFT = &BlackScholesJumpDiffusion
					::computeDensityFFT;
		}
	}

	virtual Matrix A(Real t) {
		return BlackScholes::A(t);
	}

	virtual Vector b(Real) {
		// Discretize the jump term explicitly using [1]

		// Explicit time
		Real t0 = this->time(0);

		// Compute FFT of density
		(this->*_computeDensityFFT)(t0);

		// Transformed solution
		PiecewiseLinear1 V(G, this->iterand(0));
		auto Vbar = [&] (Real x) { return V( std::exp(x) ); };

		// Build FFT vectors
		std::vector<Real> buffer;
		buffer.reserve(N);
		for(Index i = 0; i < N; i++) {
			buffer.push_back( Vbar(x0 + i * dx) );
		}
		std::vector<std::complex<Real>> bufferFFT;
		bufferFFT.reserve(N);
		fft.fwd(bufferFFT, buffer);

		// Multiplication (overwrite bufferFFT to save space)
		for(Index i = 0; i < N; i++) {
			bufferFFT[i] = bufferFFT[i]*std::conj(fprimeFFT[i]);
		}

		// Inverse fft (overwrite buffer to save space)
		fft.inv(buffer, bufferFFT);

		// Copy results to vector so that we can use existing
		// interpolation methods
		// Vector correlation = F.vector();
		//for(Index i = 0; i < N; i++) {
		//	correlation(i) = buffer[i];
		//}
		//return G.image(PiecewiseLinear1(F, correlation));

		// No need to do the above any longer since PiecewiseLinear
		// is templated to accept any structure indexable by operator[]
		// (slightly more efficient)

		PiecewiseLinear< 1, std::vector<Real> > h(
			F,                // Frequency grid
			std::move(buffer) // Data points
		);

		return G.image([&] (Real S) {
			return l(t0, S) * h(std::log(S));
		}) * dx;
	}

};

} // Modules

} // QuantPDE

#endif

