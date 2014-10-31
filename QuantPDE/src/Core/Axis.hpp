#ifndef QUANT_PDE_CORE_AXIS_HPP
#define QUANT_PDE_CORE_AXIS_HPP

#include <cassert>          // assert
#include <cmath>            // std::asinh, std::floor, std::sinh
#include <cstring>          // std::memcpy
#include <initializer_list> // std::initializer_list
#include <iostream>         // std::ostream

namespace QuantPDE {

/**
 * A set of monotonically increasing values used to represent a partition of an
 * interval (e.g. \f$\left\{x_i\right\}\f$, where
 * \f$a \equiv x_1 < \ldots < x_n \equiv b\f$; the \f$x_i\f$ are referred to
 * as ticks).
 */
class Axis final {

	Index length;
	Real *n;

	Axis() noexcept : length(0), n(nullptr) {
	}

	Axis(Index length) noexcept : length(length), n(new Real[length]) {
	}

	template <typename T>
	inline void initialize(const T &list) {
		length = list.size();
		n = new Real[length];

		assert(length > 0);

		Real *p = n;
		for(Real tick : list) {
			*(p++) = tick;
		}

		#ifndef NDEBUG
		// Make sure this axis is strictly monotonically increasing
		p = n;
		while(p < n + length - 1) {
			assert(*p < *(p+1));
			++p;
		}
		#endif
	}

public:

	/**
	 * Initializer list constructor.
	 * @param ticks The ticks.
	 */
	Axis(const std::initializer_list<Real> &ticks) noexcept {
		initialize(ticks);
	}

	/**
	 * Array constructor.
	 * @param ticks The ticks.
	 * @param length The number of ticks.
	 */
	Axis(const Real *ticks, Index length) : length(length),
			n(new Real[length]) {
		assert(length > 0);
		std::memcpy(n, ticks, sizeof(Real) * length);
	}

	/**
	 * Vector constructor.
	 * @param ticks The ticks.
	 */
	Axis(const Vector &ticks) noexcept : length(ticks.size()),
			n(new Real[length]) {
		assert(length > 0);
		std::memcpy(n, ticks.data(), sizeof(Real) * length);
	}

	Axis(const Axis &that) noexcept : length(that.length),
			n(new Real[length]) {
		std::memcpy(n, that.n, length * sizeof(Real));
	}

	Axis(Axis &&that) noexcept : length(that.length), n(that.n) {
		that.n = nullptr;
	}

	~Axis() {
		delete [] n;
		n = nullptr;
	}

	Axis &operator=(const Axis &that) & noexcept {
		length = that.length;

		Real *p = n;
		n = new Real[length];
		std::memcpy(n, that.n, length * sizeof(Real));

		delete [] p;
		p = nullptr;

		return *this;
	}

	Axis &operator=(Axis &&that) & noexcept {
		length = that.length;
		n = that.n;
		that.n = nullptr;

		return *this;
	}

	/**
	 * Return a (non-const) reference to a node by index.
	 * @param i The index.
	 */
	Real &operator[](Index i) {
		return n[i];
	}

	/**
	 * Return a const reference to a node by index.
	 * @param i The index.
	 */
	const Real &operator[](Index i) const {
		return n[i];
	}

	/**
	 * @return A const pointer to the ticks on this axis.
	 */
	const Real *ticks() const {
		return n;
	}

	/**
	 * @return The total number of ticks on this axis.
	 */
	Index size() const {
		return length;
	}

	/**
	 * Creates an axis with ticks spaced by a step-size.
	 * Similar to MATLAB's colon notation with begin:step:end.
	 * @param begin The first tick (inclusive).
	 * @param step The spacing.
	 * @param end The last tick (included if and only if end - begin is
	 *            divisible by step).
	 * @return An axis.
	 */
	static Axis range(Real begin, Real step, Real end) {
		Axis axis( (Index) std::floor( (end - begin) / step ) + 1 );
		for(Index i = 0; i < axis.length; ++i) {
			axis.n[i] = begin + step * i;
		}
		return axis;
	}

	/**
	 * Creates an axis with uniformly-spaced ticks.
	 * @param begin The first tick (inclusive).
	 * @param end The last tick (inclusive).
	 * @param points The total number of points.
	 */
	static Axis uniform(Real begin, Real end, Index points) {
		const Real dx = (end - begin) / (points - 1);
		Axis axis(points);
		for(Index i = 0; i < points; ++i) {
			axis.n[i] = begin + dx * i;
		}
		return axis;
	}

	/**
	 * Creates an axis with points clustered around a single feature.
	 * @param begin The first tick (inclusive).
	 * @param end The last tick (inclusive).
	 * @param points The total number of points.
	 * @param feature The point to cluster around.
	 * @param intensity Controls the fraction of points that lie around the
	                    feature.
	 */
	static Axis cluster(Real begin, Real end, Index points, Real feature,
			Real intensity = 1.) {
		assert(intensity > 0.);

		const Real xi_0 = asinh( (begin - feature) / intensity);
		const Real dxi = ( asinh( (end - feature) / intensity ) - xi_0 )
				/ (points - 1);

		Axis axis(points);
		for(Index i = 0; i < points; ++i) {
			axis.n[i] = feature + intensity * sinh(xi_0 + i * dxi);
		}

		return axis;
	}

	template <Index> friend class RectilinearGrid;

};

/**
 * Prettifies and prints an axis to an output stream.
 * @param os The output stream.
 * @param axis The axis.
 */
std::ostream &operator<<(std::ostream &os, const Axis &axis) {
	os << '(' << axis[0];
	for(Index i = 1; i < axis.size(); ++i) {
		os << ' ' << axis[i];
	}
	os << ')';
	return os;
}

}

#endif

