#ifndef QUANT_PDE_MODULES_BLACK_SCHOLES_HPP
#define QUANT_PDE_MODULES_BLACK_SCHOLES_HPP

#include <complex> // std::complex, std::conj
#include <utility> // std::forward
#include <vector>  // std::vector

namespace QuantPDE {

namespace Modules {

/**
 * Represents the (discretized) operator \f$\mathcal{L}\f$ defined by
 * \f$
 *     \mathcal{L} V \equiv
 *     - \frac{1}{2} \sigma^2 S^2 V_{SS} - \left( r - q \right) S V_S + r V
 * \f$
 * where \f$r\f$, \f$\sigma\f$, and \f$q\f$ are the usual interest rate,
 * volatility, and (continuous) dividend rate.
 *
 * The boundary conditions imposed are \f$ V_t - rV = 0 \f$ at \f$ S = 0 \f$ and
 * \f$ V_t - qV = 0 \f$ as \f$ S\rightarrow \infty \f$. The latter is derived by
 * assuming that the option is linear in the asset for large enough values of
 * the asset. See @cite windcliff2004analysis for discussion on the linear
 * boundary condition.
 *
 * @tparam SIndex The index of the risky asset.
**/
template <Index Dimension, Index SIndex>
class BlackScholes : public ControlledLinearSystem<Dimension> {

	static_assert(Dimension > 0, "Dimension must be positive");
	static_assert(SIndex >=0 && SIndex < Dimension,
			"The asset index must be between 0 (inclusive) and "
			"Dimension (exclusive)");

	Controllable<Dimension> r, v, q;
	Real kappa;

	void (BlackScholes::*_computeKappa)(Real);

protected:

	const RectilinearGrid<Dimension> &G;

	// TODO: Increase dimensions
	Controllable<Dimension> l;
	Noncontrollable<Dimension> g;

	void pass(Real) {
	}

	inline void computeKappa(Real t) {
		typedef AdaptiveQuadrature1<TrapezoidalRule1<>> Integral;

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
		this->registerControl(r);
		this->registerControl(v);
		this->registerControl(q);
		this->registerControl(l);
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
		kappa( 0. ),
		G( grid ),
		l( 0. ),
		g( 0. )
	{
		this->registerControl(r);
		this->registerControl(v);
		this->registerControl(q);

		_computeKappa = &BlackScholes::pass;
	}

	virtual Matrix A(Real t) {
		// 3 nonzeros per row
		Matrix M(G.size(), G.size());
		M.reserve( IntegerVector::Constant(G.size(), 3) );

		// S axis
		const Axis &S = G[SIndex];
		const Index n = S.size();

		(this->*_computeKappa)(t);

		// Take the images of curried coefficient functions
		auto rvec = G.image( curry<Dimension+1>(r, t) );
		auto vvec = G.image( curry<Dimension+1>(v, t) );
		auto qvec = G.image( curry<Dimension+1>(q, t) );
		auto lvec = G.image( curry<Dimension+1>(l, t) );

		// Interior points
		// alpha_i dt V_{i-1}^{n+1} + (1 + (alpha_i + beta_i + r) dt)
		// 		V_i^{n+1} + beta_i dt V_{i+1}^{n+1} = V_i^n

		// Boundaries
		// Left:  (1 + r dt) V_i^{n+1} = V_i^n
		// Right:              V(t, S) = g(t) S (linearity assumption)

		// Space between S ticks
		Index offset = 1;
		for(Index d = 0; d < SIndex; ++d) {
			offset *= G[d].size();
		}

		// Iterate through nodes on the grid
		for(Index idx = 0; idx < G.size(); ++idx) {
			// Retrieve index of S tick
			Index i = (idx / offset) % G[SIndex].size();

			// TODO: Remove branching
			if(i == 0) {
				// Left boundary
				// M(0, 0) = r(t, S[0]);
				M.insert(idx, idx) = rvec[idx];
			} else if(i == n-1) {
				// Right boundary
				// M(n-1, n-1) = q(t, S[n-1]);
				M.insert(idx, idx) = qvec[idx];
			} else {
				// Interior point
				const Real r_i = rvec[idx];
				const Real v_i = vvec[idx];
				const Real q_i = qvec[idx];
				const Real l_i = lvec[idx];

				const Real
					dSb = S[i]     - S[i - 1],
					dSc = S[i + 1] - S[i - 1],
					dSf = S[i + 1] - S[i]
				;

				const Real tmp1 = v_i * v_i * S[i] * S[i];
				const Real tmp2 = (r_i - q_i - l_i*kappa)*S[i];

				const Real alpha_common = tmp1 / dSb / dSc;
				const Real  beta_common = tmp1 / dSf / dSc;

				// Central
				Real alpha_i = alpha_common - tmp2 / dSc;
				Real beta_i  =  beta_common + tmp2 / dSc;
				if(alpha_i < 0.) {
					// Forward
					alpha_i = alpha_common;
					beta_i  =  beta_common + tmp2 / dSf;
				} else if(beta_i < 0.) {
					// Backward
					alpha_i = alpha_common - tmp2 / dSb;
					beta_i  =  beta_common;
				}

				// M(i, i - 1) = -alpha_i;
				// M(i, i)     = alpha_i + beta_i + r_i + l_i;
				// M(i, i + 1) = -beta_i;

				M.insert(idx, idx - offset) = -alpha_i;
				M.insert(idx, idx)          =  alpha_i + beta_i
				                               + r_i + l_i;
				M.insert(idx, idx + offset) = -beta_i;
			}
		}

		M.makeCompressed();
		return M;
	}

	virtual Vector b(Real) {
		return G.zero();
	}

	virtual bool isATheSame() const {
		return r.isConstantInTime() && v.isConstantInTime()
				&& q.isConstantInTime() && l.isConstantInTime();
	}

};

typedef BlackScholes<1, 0> BlackScholes1;

/**
 * Represents the (discretized) operator \f$ \mathcal{L} \f$ defined by
 * \f$
 *     \mathcal{L} V \equiv
 *     -\frac{1}{2} \sigma^2 S^2 V_{SS}
 *     - \left( r - q - \lambda \kappa \right) S V_S
 *     + \left( r + \lambda \right) V
 *     - \lambda \int_0^{\infty}
 *         V \left( t, S \eta \right)
 *         g \left( \eta \right)
 *     d\eta.
 * \f$
 * where \f$ r \f$, \f$ \sigma \f$, and \f$ q \f$ are the usual interest rate,
 * volatility and (continuous) rate of dividends.
 * \f$ \lambda \f$ is the mean arrival time of the Poisson process responsible
 * for generating the jumps.
 * Assuming a jump has occured, \f$ g \left( \eta \right) \f$ is the probability
 * density of a jump of amplitude \f$ \eta \f$, with
 * \f$ \kappa \equiv E \left[ \eta \right] - 1 \f$ describing the expected
 * relative change in the stock.
 *
 * The boundary conditions imposed are \f$ V_t - rV = 0 \f$ at \f$ S = 0 \f$ and
 * \f$ V_t - qV = 0 \f$ as \f$ S\rightarrow \infty \f$. The latter is derived by
 * assuming that the option is linear in the asset for large enough values of
 * the asset. See @cite windcliff2004analysis for discussion on the linear
 * boundary condition.
 *
 * The integral introduced by the jump term is handled using the FFT correlation
 * integral method described in @cite d2005robust .
**/
class BlackScholesJumpDiffusion final : public IterationNode,
		public BlackScholes1 {
		// TODO: Extend this to "multidimensional"

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

		assert(S[0] >= 0.); // Can we take its log?

		// Min:Step:Max
		const int i0 = S[0] < epsilon ? 1 : 0;
		x0 = std::log( S[i0] );
		const Real xf = std::log( S[n-2] );
		dx = (xf - x0) / (N - 1);

		// Create and return grid
		return RectilinearGrid1( Axis::uniform(x0, xf, N) );
	}

	void computeDensityFFT(Real t) {
		// Tested 2014-07-05

		//typedef TrapezoidalRule1<> Integral;
		typedef AdaptiveQuadrature1<TrapezoidalRule1<>> Integral;

		// Transformed density
		auto fbar = [&] (Real x) {
			return g(t, std::exp(x)) * std::exp(x);
		};

		// Integrate density around grid points
		std::vector<Real> fprime;
		fprime.reserve(N);
		for(Index i = 0; i <= N/2; ++i) {
			// Integrate around x_i
			const Real a = dx * (-.5 + i);
			const Real b = dx * ( .5 + i);
			fprime.push_back( Integral(fbar, a)(b) );
		}
		for(Index i = N/2+1; i < N; ++i) {
			// Integrate around x_{i - N}
			const Real a = dx * (-.5 + i - N);
			const Real b = dx * ( .5 + i - N);
			fprime.push_back( Integral(fbar, a)(b) );
		}

		// Compute FFT of transformed density
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

	// TODO: Deprecate this in lieu of something more clever
	virtual Vector b(Real) {
		// Spatial axis
		const Axis &S = G[0];
		const Index n = S.size();

		// Discretize the jump term explicitly using [1]

		// Explicit time
		Real t0 = this->time(0);

		// Compute FFT of density
		(this->*_computeDensityFFT)(t0);

		// Transformed solution
		Vector v = this->iterand(0);
		PiecewiseLinear1 V(G, v);
		auto Vbar = [&] (Real x) { return V( std::exp(x) ); };

		// Buffers
		std::vector<Real> buffer;
		std::vector<std::complex<Real>> bufferFFT;

		// Evaluate Vbar at the grid points
		buffer.reserve(N);
		for(Index i = 0; i < N; ++i) {
			buffer.push_back( Vbar(x0 + i * dx) );
		}

		// Forward transform
		bufferFFT.reserve(N);
		fft.fwd(bufferFFT, buffer);

		// Multiplication (overwrite bufferFFT to save space)
		for(Index i = 0; i < N; ++i) {
			bufferFFT[i] *= std::conj(fprimeFFT[i]);
		}

		// Inverse fft (overwrite buffer to save space)
		fft.inv(buffer, bufferFFT);

		// Copy results to vector so that we can use existing
		// interpolation methods
		/*Vector correlation = F.vector();
		for(Index i = 0; i < N; ++i) {
			correlation(i) = buffer[i];
		}
		PiecewiseLinear1 h(F, std::move(correlation));*/

		// No need to do the above any longer since PiecewiseLinear
		// is templated to accept any structure indexable by operator[]
		// (slightly more efficient)
		PiecewiseLinear< 1, std::vector<Real> > h(
			F,                // Frequency grid
			std::move(buffer) // Data points
		);

		Vector b = G.vector();

		// Left
		b(0) = 0.;

		// Interior points
		for(Index i = 1; i < n - 1; ++i) {
			b(i) = l(t0, S[i]) * h(std::log(S[i]));
		}

		// Right
		b(n - 1) = 0.;

		return b;
	}

};

typedef BlackScholesJumpDiffusion BlackScholesJumpDiffusion1;

} // Modules

} // QuantPDE

#endif

