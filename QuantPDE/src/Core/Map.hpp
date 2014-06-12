#include <cmath>       // exp, sqrt
#include <complex>     // std::complex
#include <utility>     // std::forward, std::move
#include <vector>      // std::vector

#ifndef QUANT_PDE_CORE_MAP
#define QUANT_PDE_CORE_MAP

namespace QuantPDE {

/**
 * Represents a mapping from a function on a domain (of particular dimension)
 * to a vector.
 * @see QuantPDE::Domain
 */
template <Index Dimension>
class Map {

public:

	virtual ~Map() {
	}

	/**
	 * Maps a function to a vector.
	 * @param function The function.
	 */
	virtual Vector operator()(const Function<Dimension> &function) const
			= 0;

	/**
	 * Maps a function to a vector.
	 * @param function The function (rvalue reference).
	 */
	virtual Vector operator()(Function<Dimension> &&function) const = 0;

};

typedef Map<1> Map1;
typedef Map<2> Map2;
typedef Map<3> Map3;

/**
 * Maps a function to a domain pointwise.
 */
template <Index Dimension>
class PointwiseMap final : public Map<Dimension> {

	const Domain<Dimension> *domain;

public:

	/**
	 * Constructor.
	 */
	template <typename D>
	PointwiseMap(D &domain) noexcept : domain(&domain) {
	}

	/**
	 * Copy constructor.
	 */
	PointwiseMap(const PointwiseMap &that) noexcept : domain(that.domain) {
	}

	/**
	 * Assignment operator.
	 */
	PointwiseMap &operator=(const PointwiseMap &that) & noexcept {
		domain = that.domain;
	}

	virtual Vector operator()(const Function<Dimension> &function) const {
		return domain->image(function);
	}

	virtual Vector operator()(Function<Dimension> &&function) const {
		return domain->image( std::move(function) );
	}

};

typedef PointwiseMap<1> PointwiseMap1;
typedef PointwiseMap<2> PointwiseMap2;
typedef PointwiseMap<3> PointwiseMap3;

/**
 * Performs a convolution with \f$\varphi\left(x/epsilon\right)/\epsilon\f$
 * and maps the result to the grid.
 */
class MollifierConvolution1 : public Map1 {

	const RectilinearGrid1 *grid;
	Function2 mollifier;
	Real epsilon;

	template <typename F>
	Vector map(F &&function) const {
		// 1. Create temporary grid
		// 2. Perform convolution on temporary grid
		// 3. Use linear interpolation to transfer the function to the
		//    grid

		// TODO: Fix this!

		////////////////////////////////////////////////////////////////
		// Creating the temporary Grid
		////////////////////////////////////////////////////////////////

		// Left and right endpoints
		const Axis &X = (*grid)[0];
		Real
			left  = X[0]           , // - 4 * a,
			right = X[X.size() - 1]  // + 4 * a
		;
		size_t n = (right - left) / epsilon + 1;

		// Build it
		Real *ticks = new Real[n];
		for(size_t i = 0; i < n; i++) {
			ticks[i] = left + i * epsilon;
		}
		RectilinearGrid1 tmp( Axis(ticks, n) );
		delete [] ticks;

		////////////////////////////////////////////////////////////////
		// Convolution and linear interpolation
		////////////////////////////////////////////////////////////////

		Vector u = tmp.image( std::forward<F>(function) );

		Vector v = tmp.image([&] (Real x) {
			return mollifier(epsilon, x - (left + right) / 2.);
		}) * epsilon;

		////////////////////////////////////////////////////////////////
		// FFT
		////////////////////////////////////////////////////////////////

		Eigen::FFT<Real> fft;

		std::vector<Real> u0, v0;
		for(size_t i = 0; i < n; i++) {
			u0.push_back( (u.data())[i] );
			v0.push_back( (v.data())[i] );
		}

		// Pad with zeros
		for(size_t i = 0; i < n - 1; i++) {
			u0.push_back( 0. );
			v0.push_back( 0. );
		}

		std::vector<std::complex<Real>> uk, vk;

		// Forward transform
		fft.fwd(uk, u0);
		fft.fwd(vk, v0);

		// Component-wise multiplication
		for(size_t i = 0; i < uk.size(); i++) {
			uk[i] = uk[i] * vk[i];
		}

		// Inverse transform
		fft.inv(u0, uk);

		// Transfer data to a vector
		Vector w = tmp.vector();
		for(size_t i = 0; i < n; i++) {
			w(i) = u0[i + n/2];
		}

		////////////////////////////////////////////////////////////////

		// Transfer the function to the original grid via piecewise
		// linear interpolation and return it
		return grid->image( PiecewiseLinear1(tmp, w) );
	}

public:

	/**
	 * Constructor.
	 */
	template <typename G, typename F>
	MollifierConvolution1(G &grid, F &&mollifier, Real epsilon) noexcept
			: grid(&grid), mollifier( std::forward<F>(mollifier) ),
			epsilon(epsilon) {
	}

	/**
	 * Copy constructor.
	 */
	MollifierConvolution1(const MollifierConvolution1 &that) noexcept
			: grid(that.grid), mollifier(that.mollifier),
			epsilon(that.epsilon) {
	}

	/**
	 * Move constructor.
	 */
	MollifierConvolution1(MollifierConvolution1 &&that) noexcept
			: grid(that.grid), mollifier(std::move(that.mollifier)),
			epsilon(that.epsilon) {
	}

	/**
	 * Copy assignment operator.
	 */
	MollifierConvolution1 &operator=(const MollifierConvolution1 &that) &
			noexcept {
		grid = that.grid;
		mollifier = that.mollifier;
		epsilon = that.epsilon;
		return *this;
	}

	/**
	 * Move assignment operator.
	 */
	MollifierConvolution1 &operator=(MollifierConvolution1 &&that) &
			noexcept {
		grid = that.grid;
		mollifier = std::move(that.mollifier);
		epsilon = that.epsilon;
		return *this;
	}

	/**
	 * Sets the value of \f$\epsilon\f$ to be used in the convolution.
	 * @param epsilon The new value of \f$\epsilon\f$.
	 */
	void setEpsilon(Real epsilon) {
		this->epsilon = epsilon;
	}

	virtual Vector operator()(const Function1 &function) const {
		return map(function);
	}

	virtual Vector operator()(Function1 &&function) const {
		return map( std::move(function) );
	}

};

/**
 * Uses the function \f$\varphi\left(x\right)\equiv e^{-x^2}/\sqrt{\pi}\f$ as a
 * mollifier.
 * Note that technically, this is not a mollifier (it does not have compact
 * support).
 * @see QuantPDE::MollifierConvolution1
 */
class DiracConvolution1 final : public MollifierConvolution1 {

public:

	/**
	 * Constructor.
	 */
	template <typename G>
	DiracConvolution1(G &grid, Real epsilon) noexcept
			: MollifierConvolution1(
				grid,
				[] (Real a, Real x) {
					return 1 / (a * std::sqrt(M_PI))
							* std::exp( -(x * x)
							/ (a * a) );
				},
				epsilon
			) {
	}

};


} // QuantPDE

#endif