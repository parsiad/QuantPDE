#ifndef QUANT_PDE_CORE_AXIS
#define QUANT_PDE_CORE_AXIS

#include <cassert>          // assert
#include <cmath>            // std::floor
#include <cstring>          // std::memcpy
#include <initializer_list> // std::initializer_list
#include <iostream>         // std::ostream

namespace QuantPDE {

/**
 * A set of monotonically increasing values used to represent a partition of an
 * interval (e.g. the set \f$\left\{x_i\right\}\f$, where
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

public:

	/**
	 * Constructor.
	 * @param list The ticks. These should be strictly monotonically
	 *             increasing.
	 */
	Axis(const std::initializer_list<Real> &list) noexcept
			: length(list.size()), n(new Real[length]) {
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
			p++;
		}
		#endif
	}

	/**
	 * Constructor.
	 * @param ticks Ticks (as an array).
	 * @param length The number of ticks.
	 */
	Axis(const Real *ticks, Index length) : length(length),
			n(new Real[length]) {
		assert(length > 0);
		std::memcpy(n, ticks, sizeof(Real) * length);
	}

	/**
	 * Initialize the axis from a vector.
	 * @param vector The vector.
	 */
	Axis(const Vector &vector) noexcept : length(vector.size()),
			n(new Real[length]) {
		assert(length > 0);
		std::memcpy(n, vector.data(), sizeof(Real) * length);
	}

	/**
	 * Copy constructor.
	 */
	Axis(const Axis &that) noexcept : length(that.length),
			n(new Real[length]) {
		std::memcpy(n, that.n, length * sizeof(Real));
	}

	/**
	 * Move constructor.
	 */
	Axis(Axis &&that) noexcept : length(that.length), n(that.n) {
		that.n = nullptr;
	}

	/**
	 * Destructor.
	 */
	~Axis() {
		delete [] n;
	}

	/**
	 * Copy assignment operator.
	 */
	Axis &operator=(const Axis &that) & noexcept {
		length = that.length;

		Real *p = n;
		n = new Real[length];
		std::memcpy(n, that.n, length * sizeof(Real));

		delete [] p;

		return *this;
	}

	/**
	 * Move assignment operator.
	 */
	Axis &operator=(Axis &&that) & noexcept {
		length = that.length;
		n = that.n;
		that.n = nullptr;

		return *this;
	}

	/**
	 * Return a non-const reference to a node by index.
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
	 * @return A pointer to the ticks on this axis.
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
	 * Creates an axis with evenly spaced ticks. Similar to MATLAB's colon
	 * notation with begin:step:end.
	 * @param begin The first tick (included in the axis).
	 * @param step The spacing.
	 * @param end The last tick (included if and only if end - begin is
	 *            divisible by step).
	 * @return An axis.
	 */
	static Axis range(Real begin, Real step, Real end) {
		Axis axis( (Index) std::floor( (end - begin) / step ) + 1 );
		for(Index i = 0; i < axis.length; i++) {
			axis[i] = begin + step * i;
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
	for(Index i = 1; i < axis.size(); i++) {
		os << ' ' << axis[i];
	}
	os << ')';
	return os;
}

}

#endif

