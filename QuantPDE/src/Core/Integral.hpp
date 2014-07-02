#ifndef QUANT_PDE_CORE_INTEGRAL
#define QUANT_PDE_CORE_INTEGRAL

#include <array>   // std::array
#include <cstdint> // std::intmax_t
#include <cstdlib> // std::abs
#include <utility> // std::forward, std::move

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
	 * @param x The upper coordinates of the integral.
	 * @return The integral evaluated at x.
	 */
	virtual Real integrate(const std::array<Real, Dimension> &x) const = 0;

public:

	/**
	 * Destructor.
	 */
	virtual ~Integral() {
	}

	/**
	 * Performs the integration.
	 * @param x The upper coordinates of the integral.
	 * @return The integral evaluated at x.
	 */
	template <typename ...Ts>
	Real operator()(Ts ...x) const {
		return integrate( {{x...}} );
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
	static inline Real integrate(const F &function, const Real *a,
			const Real *x) {
		Real sum;

		// Left endpoint
		sum = function(*a) * TrapezoidalRuleSubroutine<N-1, Is...>
				::integrate(function, a+1, x+1);

		// Interior points
		const Real dx = (*x - *a) / I;
		for(int i = 1; i <= I - 1; i++) {
			sum += 2 * function(*a + i * dx)
					* TrapezoidalRuleSubroutine<N-1, Is...>
					::integrate(function, a+1, x+1);
		}

		// Right endpoint
		sum += function(*x) * TrapezoidalRuleSubroutine<N-1, Is...>
				::integrate(function, a+1, x+1);

		return sum;
	}
};

template <>
struct TrapezoidalRuleSubroutine<0> {
	template <typename F>
	static inline Real integrate(const F &, const Real *, const Real *) {
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

	virtual Real integrate(const std::array<Real, Dimension> &x) const {
		// 2^n
		typedef IntegerPower<2, Dimension> Power;

		// Product of N_1 * N_2 * ... * N_n
		typedef IntegerProduct<Intervals...> Product;

		double scale = 1. / Power::value;
		for(int i = 0; i < Dimension; i++) {
			scale *= (x[i] - a[i]) / Product::value;
		}

		return scale*TrapezoidalRuleSubroutine<Dimension, Intervals...>
				::integrate(function, a.data(), x.data());
	}

	const Function<Dimension> function;
	const std::array<Real, Dimension> a;

public:

	/**
	 * Constructor.
	 * @param function The function to integrate.
	 * @param a The lower coordinates of the integral.
	 */
	template <typename F, typename ...Ts>
	TrapezoidalRule(F &&function, Ts ...a) noexcept
			: function(std::forward<F>(function)), a( {{a...}} ) {
	}

	/**
	 * Copy constructor.
	 */
	TrapezoidalRule(const TrapezoidalRule &that) noexcept
			: function(that.function), a(that.a) {
	}

	/**
	 * Move constructor.
	 */
	TrapezoidalRule(TrapezoidalRule &&that) noexcept
			: function(std::move(that.function)),
			a(std::move(that.a)) {
	}

	/**
	 * Copy assignment.
	 */
	TrapezoidalRule &operator=(const TrapezoidalRule &that) noexcept {
		function = that.function;
		a = that.a;
	}

	/**
	 * Move assignment.
	 */
	TrapezoidalRule &operator=(TrapezoidalRule &&that) noexcept {
		function = std::move(that.function);
		a = std::move(that.a);
	}

};

template <size_t ...Intervals>
using TrapezoidalRule1 = TrapezoidalRule<1, Intervals...>;

template <size_t ...Intervals>
using TrapezoidalRule2 = TrapezoidalRule<2, Intervals...>;

template <size_t ...Intervals>
using TrapezoidalRule3 = TrapezoidalRule<3, Intervals...>;

////////////////////////////////////////////////////////////////////////////////

template <Index Dimension, typename T>
class AdaptiveQuadrature : public Integral<Dimension> {

	static constexpr int maxDepth = 64;

	////////////////////////////////////////////////////////////////////////

	template <int ...Indices>
	inline Real packAndCall(
		const Real *array,
		Sequence<Indices...>
	) const {
		return T(function, array[Indices]...)(
				array[Indices + Dimension]...);
	}

	////////////////////////////////////////////////////////////////////////

	// The buffer Real *p is recycled by each call to adapt
	Real refine(Real previous, Real *p, int n = maxDepth) const {
		/*
		for(int i = 0; i < maxDepth - n; i++) {
			std::cout << "\t";
		}
		std::cout << "Integrating cell ";
		for(Index i = 0; i < Dimension; i++) {
			std::cout << p[i] << " -> " << p[i+Dimension] << " ";
		}
		std::cout << "(" << previous << ")" << std::endl;
		*/

		// 2^Dimension
		typedef IntegerPower<2, Dimension> Power;

		Real pp[Power::value - 1][Dimension * 2];
		Real integrals[Power::value];
		Real sum = 0.;

		// Indicies 0 to 2^N - 1
		for(std::intmax_t i = 0; i < Power::value - 1; i++) {
			// Build array
			for(Index j = 0; j < Dimension; j++) {
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
			for(Index j = 0; j < Dimension; j++) {
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

		/*
		for(std::intmax_t i = 0; i < Power::value - 1; i++) {
			for(Index j = 0; j < Dimension; j++) {
				for(int k = 0; k < maxDepth - n; k++) {
					std::cout << "\t";
				}
				std::cout << pp[i][j] << " -> "
						<< pp[i][j+Dimension] << " ";
			}
			std::cout << "(" << integrals[i] << ")" << std::endl;
		}
		for(Index j = 0; j < Dimension; j++) {
			for(int k = 0; k < maxDepth - n; k++) {
				std::cout << "\t";
			}
			std::cout << p[j] << " -> " << p[j+Dimension] << " " ;
		}
		std::cout << "(" << integrals[Power::value - 1] << ")"
				<< std::endl;
		for(int k = 0; k < maxDepth - n; k++) {
			std::cout << "\t";
		}
		std::cout << "Sum: " << sum << std::endl;

		for(int k = 0; k < maxDepth - n; k++) {
			std::cout << "\t";
		}
		*/

		const Real error = std::abs((sum - previous)/sum);

		//std::cout << "Error: " << error << " "
		//		<< (error > tolerance) << std::endl;

		if( n > 0 && error > tolerance ) {
			sum = 0.;

			// Most of the cells
			for(std::intmax_t i = 0; i < Power::value - 1; i++) {
				sum += refine(integrals[i], pp[i], n - 1);
			}

			// Top-right cell
			sum += refine(integrals[Power::value - 1], p, n - 1);
		}
		/*
		else {
			for(int k = 0; k < maxDepth - n; k++) {
				std::cout << "\t";
			}
			std::cout << "Match" << std::endl;
		}

		for(int i = 0; i < maxDepth - n; i++) {
			std::cout << "\t";
		}
		std::cout << "Result: " << sum << std::endl << std::endl;
		*/

		return sum;
	}

	virtual Real integrate(const std::array<Real, Dimension> &x) const {
		// Outermost integration
		Real p[Dimension * 2];
		for(Index j = 0; j < Dimension; j++) {
			p[j] = a[j];
			p[j+Dimension] = x[j];
		}

		const Real integral = packAndCall(p,
				GenerateSequence<Dimension>());

		return refine(integral, p);
	}

	const Function<Dimension> function;
	const std::array<Real, Dimension> a;
	const Real tolerance;

public:

	// TODO: Use SFINAE to pass in tolerance

	/**
	 * Constructor.
	 * @param function The function to integrate.
	 * @param a The lower coordinates of the integral.
	 */
	template <typename F, typename ...Ts>
	AdaptiveQuadrature(F &&function, Ts ...a) noexcept :
		function(std::forward<F>(function)),
		a( {{a...}} ),
		tolerance(1e-6)
	{
	}

	/**
	 * Copy constructor.
	 */
	AdaptiveQuadrature(const AdaptiveQuadrature &that) noexcept
			: function(that.function), a(that.a),
			tolerance(that.tolerance) {
	}

	/**
	 * Move constructor.
	 */
	AdaptiveQuadrature(AdaptiveQuadrature &&that) noexcept
			: function(std::move(that.function)),
			a(std::move(that.a)), tolerance(that.tolerance) {
	}

	/**
	 * Copy assignment.
	 */
	AdaptiveQuadrature &operator=(const AdaptiveQuadrature &that) noexcept {
		function = that.function;
		a = that.a;
		tolerance = that.tolerance;
	}

	/**
	 * Move assignment.
	 */
	AdaptiveQuadrature &operator=(AdaptiveQuadrature &&that) noexcept {
		function = std::move(that.function);
		a = std::move(that.a);
		tolerance = that.tolerance;
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

