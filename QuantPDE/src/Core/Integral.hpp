#ifndef QUANT_PDE_CORE_INTEGRAL_HPP
#define QUANT_PDE_CORE_INTEGRAL_HPP

#include <array>   // std::array
#include <cstdint> // std::intmax_t
#include <cstdlib> // std::abs
#include <limits>  // std::numeric_limits
#include <utility> // std::forward, std::move

// TODO: static integrate method

namespace QuantPDE {

/**
 * A function defined by
 * \f$F\left(\mathbf{x}\right)\equiv\int_{\mathbf{a}}^{\mathbf{x}}f\left(\mathbf{y}\right)d\mathbf{y}\f$,
 * where \f$f\colon\mathbb{R}^{n}\rightarrow\mathbb{R}\f$ and \f$\mathbb{a}\f$
 * are fixed.
 * @tparam Dimension \f$n\f$.
 */
template <Index Dimension>
class Integral {

	/**
	 * Performs the integration.
	 * @param a Lower bounds of integration.
	 * @param x Upper bounds of integration.
	 * @return The result.
	 */
	virtual Real compute(const std::array<Real, Dimension> &a,
			const std::array<Real, Dimension> &x) const = 0;

protected:

	const std::array<Real, Dimension> a;
	const Function<Dimension> function;

public:

	/**
	 * Constructor.
	 * @param function Function to integrate.
	 * @param a Lower bounds of integration.
	 */
	template <typename F, typename ...Ts>
	Integral(F &&function, Ts ...a) noexcept
			: a( {{a...}} ), function(std::forward<F>(function)) {
	}

	/**
	 * Copy constructor.
	 */
	Integral(const Integral &that) noexcept : a(that.a),
			function(that.function) {
	}

	/**
	 * Move constructor.
	 */
	Integral(Integral &&that) noexcept : a(std::move(that.a)),
			function(std::move(that.function)) {
	}

	/**
	 * Copy assignment operator.
	 */
	Integral &operator=(const Integral &that) noexcept {
		a = that.a;
		function = that.function;
		return *this;
	}

	/**
	 * Move assignment operator.
	 */
	Integral &operator=(Integral &&that) noexcept {
		a = std::move(that.a);
		function = std::move(that.function);
		return *this;
	}

	/**
	 * Destructor.
	 */
	virtual ~Integral() {
	}

	/**
	 * Performs the integration.
	 * @param x Upper bounds of integration.
	 * @return The integral evaluated at x.
	 */
	template <typename ...Ts>
	Real operator()(Ts ...x) const {
		// TODO: Optimize; inefficient to figure out whether an integral
		// is finite at run-time

		const Real tolerance = 1e-6; // TODO: Make this a parameter

		// Buffers for lower and upper bounds of integration
		std::array<Real, Dimension> aa(a);
		std::array<Real, Dimension> xx {{x...}};

		// neginf[i] == true => i-th lower bound is -infinity
		// posinf[i] == true => i-th upper bound is +infinity
		bool neginf[Dimension], posinf[Dimension];

		// The integral consists of only finite lower and upper bounds
		bool finite = true;

		// Used to expand the region of integration on each axis
		Real m[Dimension];

		for(Index i = 0; i < Dimension; ++i) {
			neginf[i] = a[i]  ==
					-std::numeric_limits<Real>::infinity();
			posinf[i] = xx[i] ==
					std::numeric_limits<Real>::infinity();

			if(neginf[i]) {
				if(posinf[i]) {
					// (-inf, inf)
					aa[i] = -1.;
					xx[i] =  1.;
					m[i]  =  0.;

					finite = false;
				} else {
					// (-inf, x]
					aa[i] = xx[i] - 1.;
					m[i]  = xx[i];

					finite = false;
				}
			} else if(posinf[i]) {
				// [a, inf)
				xx[i] = aa[i] + 1.;
				m[i]  = aa[i];

				finite = false;
			}
		}

		// No need to integrate more than once if the region of
		// integration does not change
		if(finite) {
			return compute(aa, xx);
		}

		// TODO: Currently, we are reintegrating over the whole region
		//       of integration at every iteration. This is redundant!
		//       Can save roughly 1/2^N of the time if we reuse the
		//       integral over the subregion. This is particularly
		//       noticeable when N = 1 (one dimension).

		// Keep expanding the region of integration until convergence is
		// attained
		Real previous = compute(aa, xx);
		Real integral;
		while(true) {
			// Expand the region of integration
			for(Index i = 0; i < Dimension; ++i) {
				if(neginf[i]) {
					aa[i] += aa[i] - m[i];
				}

				if(posinf[i]) {
					xx[i] += xx[i] - m[i];
				}
			}

			// Recompute the integral
			integral = compute(aa, xx);

			// Break when the relative error is low enough
			if(std::abs( (integral - previous) / integral )
					< tolerance) {
				break;
			}

			previous = integral;
		}

		return integral;
	}

};

////////////////////////////////////////////////////////////////////////////////

/**
 * Computes the sum appearing in the uniform trapezoidal rule efficiently for
 * arbitrary dimensionality.
 * @tparam N Dimension.
 * @tparam Is The number of intervals to use per dimension.
 * @see QuantPDE::TrapezoidalRule
 */
template <Index N, int ...Is>
struct TrapezoidalRuleSubroutine {
};

template <Index N, int I, int ...Is>
struct TrapezoidalRuleSubroutine<N, I, Is...> {
	static_assert(I > 0, "The number of intervals must be positive");

	template <typename F>
	static inline Real compute(const F &function, const Real *a,
			const Real *x) {
		Real sum;

		// Left endpoint
		sum = function(*a) * TrapezoidalRuleSubroutine<N-1, Is...>
				::compute(function, a+1, x+1);

		// Interior points
		const Real dx = (*x - *a) / I;
		for(int i = 1; i <= I - 1; ++i) {
			sum += 2 * function(*a + i * dx)
					* TrapezoidalRuleSubroutine<N-1, Is...>
					::compute(function, a+1, x+1);
		}

		// Right endpoint
		sum += function(*x) * TrapezoidalRuleSubroutine<N-1, Is...>
				::compute(function, a+1, x+1);

		return sum;
	}
};

template <>
struct TrapezoidalRuleSubroutine<0> {
	template <typename F>
	static inline Real compute(const F &, const Real *, const Real *) {
		return 1.;
	}
};

/**
 * Performs the trapezoidal rule assuming a uniform grid.
 * In one dimension, the trapezoidal rule (on the uniform grid
 * \f$a\equiv x_{1}<x_{2}<\ldots<x_{N+1}\equiv x\f$) is
 * \f$\int_{a}^{x}f\left(y\right)dy\approx\frac{x-a}{2N}\left(f\left(x_{1}\right)+2f\left(x_{2}\right)+\ldots+2f\left(x_{N}\right)+f\left(x_{N+1}\right)\right)\f$.
 * @tparam Dimension \f$n\f$.
 * @tparam Intervals The number of intervals to use per dimension.
 * @see QuantPDE::Integral
 */
template <Index Dimension, int ...Intervals>
class TrapezoidalRule : public Integral<Dimension> {

	static_assert(Dimension == sizeof...(Intervals),
			"The number of arguments must be consistent with the "
			"dimensions");

	virtual Real compute(const std::array<Real, Dimension> &a,
			const std::array<Real, Dimension> &x) const {
		// 2^n
		typedef IntegerPower<2, Dimension> Power;

		// Product of N_1 * N_2 * ... * N_n
		typedef IntegerProduct<Intervals...> Product;

		double scale = 1. / Power::value;
		for(int i = 0; i < Dimension; ++i) {
			scale *= (x[i] - a[i]) / Product::value;
		}

		return scale*TrapezoidalRuleSubroutine<Dimension, Intervals...>
				::compute(this->function, a.data(), x.data());
	}

public:

	/**
	 * Constructor.
	 */
	template <typename F, typename ...Ts>
	TrapezoidalRule(F &&function, Ts ...a) noexcept
			: Integral<Dimension>(std::forward<F>(function), a...) {
	}

	/**
	 * Copy constructor.
	 */
	TrapezoidalRule(const TrapezoidalRule &that) noexcept
			: Integral<Dimension>(that) {
	}

	/**
	 * Move constructor.
	 */
	TrapezoidalRule(TrapezoidalRule &&that) noexcept
			: Integral<Dimension>(std::move(that)) {
	}

	/**
	 * Copy assignment operator.
	 */
	TrapezoidalRule &operator=(const TrapezoidalRule &that) noexcept {
		Integral<Dimension>::operator=(that);
		return *this;
	}

	/**
	 * Move assignment operator.
	 */
	TrapezoidalRule &operator=(TrapezoidalRule &&that) noexcept {
		Integral<Dimension>::operator=(std::move(that));
		return *this;
	}

};

template <int Intervals1 = 1>
using TrapezoidalRule1 = TrapezoidalRule<1, Intervals1>;

template <int Intervals1 = 1, int Intervals2 = 1>
using TrapezoidalRule2 = TrapezoidalRule<2, Intervals1, Intervals2>;

template <int Intervals1 = 1, int Intervals2 = 1, int Intervals3 = 1>
using TrapezoidalRule3 = TrapezoidalRule<3, Intervals1, Intervals2, Intervals3>;

////////////////////////////////////////////////////////////////////////////////

template <Index Dimension, typename T>
class AdaptiveQuadrature : public Integral<Dimension> {

	static constexpr int maxDepth = 64;

	////////////////////////////////////////////////////////////////////////

	template <int ...Indices>
	inline Real packAndCall(const Real *array, Sequence<Indices...>) const {
		return T(this->function, array[Indices]...)(
				array[Indices + Dimension]...);
	}

	////////////////////////////////////////////////////////////////////////

	Real refine(Real previous, Real *p, int n = maxDepth) const {
		// The buffer Real *p is recycled by each call to adapt

		// 2^Dimension
		typedef IntegerPower<2, Dimension> Power;

		Real pp[Power::value - 1][Dimension * 2];
		Real integrals[Power::value];
		Real sum = 0.;

		// Indices 0 to 2^N - 1
		for(std::intmax_t i = 0; i < Power::value - 1; ++i) {
			// Build array
			for(Index j = 0; j < Dimension; ++j) {
				const Real a = p[j];
				const Real x = p[j+Dimension];

				if(i && (1 << j)) {
					// j-th bit of i is 1
					pp[i][j] = (a + x) / 2.;
					pp[i][j+Dimension] = x;
				} else {
					// j-th bit of i is 0
					pp[i][j] = a;
					pp[i][j+Dimension] = (a + x) / 2.;
				}
			}

			// Integrate
			integrals[i] = packAndCall(
				pp[i],
				GenerateSequence<Dimension>()
			);

			sum += integrals[i];
		}

		// Index Power::value - 1
		{
			for(Index j = 0; j < Dimension; ++j) {
				const Real a = p[j];
				const Real x = p[j+Dimension];

				p[j] = (a + x) / 2.;

				// No need to set this since it will already be
				// the top-right corner
				// p[j+Dimension] = x;
			}

			// Integrate
			integrals[Power::value - 1] = packAndCall(
				p,
				GenerateSequence<Dimension>()
			);

			sum += integrals[Power::value - 1];
		}

		const Real error = std::abs((sum - previous)/sum);

		if( n > 0 && error > tolerance ) {
			sum = 0.;

			// Most of the cells
			for(std::intmax_t i = 0; i < Power::value - 1; ++i) {
				sum += refine(integrals[i], pp[i], n - 1);
			}

			// Top-right cell
			sum += refine(integrals[Power::value - 1], p, n - 1);
		}

		return sum;
	}

	virtual Real compute(const std::array<Real, Dimension> &a,
			const std::array<Real, Dimension> &x) const {
		// Outermost integration
		Real p[Dimension * 2];
		for(Index j = 0; j < Dimension; ++j) {
			p[j] = a[j];
			p[j+Dimension] = x[j];
		}

		const Real integral = packAndCall(p,
				GenerateSequence<Dimension>());

		return refine(integral, p);
	}

	const Real tolerance;

public:

	/**
	 * Constructor.
	 */
	template <typename F, typename ...Ts>
	AdaptiveQuadrature(F &&function, Ts ...a) noexcept
			: Integral<Dimension>(std::forward<F>(function), a...),
			tolerance(1e-6) {
		// TODO: Use SFINAE to pass in tolerance
	}

	/**
	 * Copy constructor.
	 */
	AdaptiveQuadrature(const AdaptiveQuadrature &that) noexcept
			: Integral<Dimension>(that) {
	}

	/**
	 * Move constructor.
	 */
	AdaptiveQuadrature(AdaptiveQuadrature &&that) noexcept :
			Integral<Dimension>(std::move(that)),
			tolerance(that.tolerance) {
	}

	/**
	 * Copy assignment.
	 */
	AdaptiveQuadrature &operator=(const AdaptiveQuadrature &that) noexcept {
		Integral<Dimension>::operator=(that);
		tolerance = that.tolerance;
		return *this;
	}

	/**
	 * Move assignment.
	 */
	AdaptiveQuadrature &operator=(AdaptiveQuadrature &&that) noexcept {
		Integral<Dimension>::operator=(std::move(that));
		tolerance = that.tolerance;
		return *this;
	}

};

template <typename T>
using AdaptiveQuadrature1 = AdaptiveQuadrature<1, T>;

template <typename T>
using AdaptiveQuadrature2 = AdaptiveQuadrature<2, T>;

template <typename T>
using AdaptiveQuadrature3 = AdaptiveQuadrature<3, T>;

} // QuantPDE

#endif

