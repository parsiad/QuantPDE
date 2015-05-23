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

	static Axis nonsymmetricCluster(Real begin, Real feature, Real end,
			Index points, Real intensity = 1.) {
		assert(intensity > 0.);

		assert(  begin <= feature);
		assert(feature <= end    );
		assert(begin   <  end    );

		Real K = (feature - begin) / (end - begin);

		bool flipped = false;
		if(K < 0.5) {
			flipped = true;
			K  =  1. - K;
		}

		const Real c = K / intensity;

		const Real xi_0 =   asinh( (0. - K) / c );
		const Real dxi  = ( asinh( (1. - K) / c ) - xi_0 ) / (points-1);

		Axis axis(points);
		axis.n[0] = begin;
		if(flipped) {
			for(Index i = 1; i < points - 1; ++i) {
				axis.n[points-1 - i] = (1. -
					( K + c * sinh(xi_0 + i * dxi) )
				) * (end - begin) + begin;
			}
		} else {
			for(Index i = 1; i < points - 1; ++i) {
				axis.n[           i] = (
					( K + c * sinh(xi_0 + i * dxi) )
				) * (end - begin) + begin;
			}
		}
		axis.n[points - 1] = end;

		return axis;
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
	 * Vector constructor.
	 * @param ticks The ticks.
	 */
	Axis(const std::vector<Real> &ticks) noexcept {
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
	 * @param feature The point to cluster around.
	 * @param end The last tick (inclusive).
	 * @param points The total number of points (if this quantity is even,
	 *               an axis with points+1 points is returned instead).
	 * @param intensity Controls the fraction of points that lie around the
	                    feature.
	 */
	static Axis cluster(Real begin, Real feature, Real end, Index points,
			Real intensity = 1.) {
		Index p = (points + 1 + (points % 2 == 0 ? 1 : 0)) / 2;
		return
			nonsymmetricCluster(
				begin,
				feature,
				feature,
				p,
				intensity
			) +
			nonsymmetricCluster(
				feature,
				feature,
				end,
				p,
				intensity
			)
		;
	}

	/**
	 * A hand-picked axis from 0 to 100 with 33 ticks, clustered close to 1.
	 */
	static Axis special;

	/**
	 * Union of two axes.
	 * @param a One axis.
	 * @param b Another axis.
	 * @return An axis with points from both input axes.
	 */
	friend Axis operator+(const Axis &a, const Axis &b) {
		// Assumes both input axes are strictly increasing

		// First-pass (count unique elements)
		Index i = 0, j = 0, k = 0;
		while(true) {
			if(i >= a.length) {
				// Add all remaining elements in b
				k += b.length - j;
				break;
			} else if(j >= b.length) {
				// Add all remaining elements in a
				k += a.length - i;
				break;
			}

			if     (a.n[i] == b.n[j]) { ++i; ++j; }
			else if(a.n[i]  < b.n[j]) { ++i;      }
			else                      {      ++j; }

			++k;
		}

		// Second-pass (populate axis)
		Axis c(k);
		i = j = k = 0;
		while(true) {
			if(i >= a.length) {
				// Add all remaining elements in b
				while(j < b.length) { c.n[k++] = b.n[j++]; }
				break;
			} else if(j >= b.length) {
				// b is finished
				while(i < a.length) { c.n[k++] = a.n[i++]; }
				break;
			}

			if     (a.n[i] == b.n[j]) { c.n[k++] = a.n[i++]; ++j; }
			else if(a.n[i]  < b.n[j]) { c.n[k++] = a.n[i++];      }
			else                      { c.n[k++] = b.n[j++];      }
		}

		return c;
	}

	/**
	 * Scales the points on this axis by a constant.
	 * @param axis The axis.
	 * @param c The constant.
	 */
	friend Axis &&operator*(Axis &&axis, Real c) {
		assert(c != 0.);
		for(Index i = 0; i < axis.length; ++i) {
			axis.n[i] *= c;
		}
		return std::move(axis);
	}

	/**
	 * Scales the points on this axis by a constant.
	 * @param axis The axis.
	 * @param c The constant.
	 */
	friend Axis operator*(const Axis &axis, Real c) {
		assert(c != 0.);
		Axis newAxis(axis.length);
		for(Index i = 0; i < axis.length; ++i) {
			newAxis.n[i] = axis.n[i] * c;
		}
		return newAxis;
	}

	/**
	* Scales the points on this axis by a constant.
	* @param c The constant.
	* @param axis The axis.
	*/
	friend Axis &&operator*(Real c, Axis &&axis) {
		return std::move(axis) * c;
	}

	/**
	 * Scales the points on this axis by a constant.
	 * @param axis The axis.
	 * @param c The constant.
	 */
	friend Axis operator*(Real c, const Axis &axis) {
		return axis * c;
	}

	/**
	* Translate the points on this axis by a constant.
	* @param axis The axis.
	* @param c The constant.
	*/
	friend Axis &&operator+(Axis &&axis, Real c) {
		for(Index i = 0; i < axis.length; ++i) {
			axis.n[i] += c;
		}
		return std::move(axis);
	}

	/**
	 * Translate the points on this axis by a constant.
	 * @param axis The axis.
	 * @param c The constant.
	 */
	friend Axis operator+(const Axis &axis, Real c) {
		Axis newAxis(axis.length);
		for(Index i = 0; i < axis.length; ++i) {
			newAxis.n[i] = axis.n[i] + c;
		}
		return newAxis;
	}

	/**
	* Translate the points on this axis by a constant.
	* @param c The constant.
	* @param axis The axis.
	*/
	friend Axis &&operator+(Real c, Axis &&axis) {
		return std::move(axis) + c;
	}

	/**
	 * Translate the points on this axis by a constant.
	 * @param axis The axis.
	 * @param c The constant.
	 */
	friend Axis operator+(Real c, const Axis &axis) {
		return axis + c;
	}

	/**
	* Translate the points on this axis by the negative of a constant.
	* @param axis The axis.
	* @param c The constant.
	*/
	friend Axis &&operator-(Axis &&axis, Real c) {
		return std::move(axis) + (-c);
	}

	/**
	* Translate the points on this axis by the negative of a constant.
	* @param axis The axis.
	* @param c The constant.
	*/
	friend Axis operator-(const Axis &axis, Real c) {
		return axis + (-c);
	}

	/**
	* Translate the points on this axis by the negative of a constant.
	* @param c The constant.
	* @param axis The axis.
	*/
	friend Axis &&operator-(Real c, Axis &&axis) {
		return std::move(axis) - c;
	}

	/**
	* Translate the points on this axis by the negative of a constant.
	* @param axis The axis.
	* @param c The constant.
	*/
	friend Axis operator-(Real c, const Axis &axis) {
		return axis + (-c);
	}

	template <Index> friend class RectilinearGrid;

};

Axis Axis::special {
	0., .1, .2, .3, .4, .5, .6, .7,
	/*.75,*/ .80,
	.84, .88, .92,
	.94, .96, .98, 1., 1.02, 1.04, 1.06, 1.08, 1.10,
	1.14, 1.18,
	1.23,
	1.30, 1.40, 1.50,
	1.75,
	2.25,
	3.00,
	7.50,
	20.,
	100.
};

/**
 * Prettifies and prints an axis to an output stream.
 * @param os The output stream.
 * @param axis The axis.
 */
std::ostream &operator<<(std::ostream &os, const Axis &axis) {
	os << axis[0];
	for(Index i = 1; i < axis.size(); ++i) {
		os << ' ' << axis[i];
	}
	return os;
}

}

#endif
