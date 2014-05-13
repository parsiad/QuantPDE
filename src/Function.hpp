#ifndef QUANT_PDE_FUNCTION_HPP
#define QUANT_PDE_FUNCTION_HPP

#include <cstring> // std::memcpy

namespace QuantPDE {

typedef std::function< NRealToReal<1> > Function1;
typedef std::function< NRealToReal<2> > Function2;
typedef std::function< NRealToReal<3> > Function3;
typedef std::function< NRealToReal<4> > Function4;

template <Index N> // TODO: Check if positive
using Function = std::function< NRealToReal<N> >;

/**
 * A function \f$f\f$ defined piecewise whose pieces are affine functions.
 * Suppose \f$f\f$ is defined on a grid composed of \f$n\f$ axes, and that the
 * \f$i\f$-th axis is a partition consisting of \f$m_i\f$ nodes.
 * Querying \f$f\f$ at a point involves a binary search on each axis.
 * This procedure has worst-case complexity
 * \f$\sim \sum_i \lg m_i \leq n \lg m\f$, where
 * \f$m \equiv \max\left\{m_i\right\}\f$.
 */
template <Index dim>
class PiecewiseLinear {

	const Grid<dim> *grid;
	const Vector *vector;

	Real interpolate(Index *indices, const Real *weights, Index shift,
			Index n = 0) const {
		// Base case
		if(n == dim) {
			return grid->accessor(*vector)(indices);
		}

		Index *stride = indices + shift;
		std::memcpy(stride, indices, sizeof(Index) * dim);
		stride[n]++;

		return weights[n] * interpolate(indices, weights, shift / 2,
				n + 1) + (1 - weights[n]) * interpolate(stride,
				weights, shift / 2, n + 1);
	}

	template <typename ...Args>
	Real get(Args... args) const {
		static_assert(dim == sizeof...(Args),
				"The number of arguments must be consistent "
				"with the dimensions");

		typedef IntegerPower<2, dim> dimpow;
		static_assert(!dimpow::overflow, "Overflow detected");

		Real arguments[] {args...};
		Real *coordinates = arguments;

		Real weights[dim];
		Index indices[dimpow::value];

		// For the i-th coordinate, find the ticks on the i-th axis that
		// it lies between along with the distance from the leftmost
		// tick
		for(Index i = 0; i < dim; i++, coordinates++) { // Unroll loop
			const Axis &x = (*grid)[i];
			Index length = x.size();

			if(*coordinates <= x[0]) {
				indices[i] = 0;
				weights[i] = 1.;
				continue;
			}

			if(*coordinates >= x[length - 1]) {
				indices[i] = length - 2;
				weights[i] = 0.;
				continue;
			}

			// Binary search to find tick
			Index lo = 0, hi = length - 2, mid = 0;
			Real weight = 0.;
			while(lo <= hi) {
				mid = (lo + hi) / 2;
				if(*coordinates < x[mid]) {
					hi = mid - 1;
				} else if(*coordinates >= x[mid + 1]) {
					lo = mid + 1;
				} else {
					weight = (x[mid + 1] - *coordinates)
							/ (x[mid + 1] - x[mid]);
					break;
				}
			}

			indices[i] = mid;
			weights[i] = weight;
		}

		////////////////////////////////////////////////////////////////
		// Recursive version
		////////////////////////////////////////////////////////////////

		// return interpolate(indices, weights, dimpow::value / 2);

		////////////////////////////////////////////////////////////////
		// Nonrecursive version
		////////////////////////////////////////////////////////////////

		// At first glance, the bit-shifting method below seems slow,
		// but since 2^dim is a small number (e.g. 8 for a 3 dimensional
		// PDE), the loop is unrolled and this is much faster than using
		// the call to interpolate(...) above.

		// Past a certain dimension-threshold, the call to
		// interpolate(...) will be more efficient as it uses Horner's
		// form.

		auto accessor = grid->accessor(*vector);
		Index idxs[dim];
		Real interpolated = 0.;

		// Unroll loop
		for(Index i = 0; i < dimpow::value; i++) {

			Real factor = 1.;

			// Unroll loop
			for(Index j = 0; j < dim; j++) {
				if(i && (1 << j)) {
					// j-th bit of i is 1
					factor *= weights[j];
					idxs[j] = indices[j];
				} else {
					// j-th bit of i is 0
					factor *= 1 - weights[j];
					idxs[j] = indices[j] + 1;
				}
			}

			interpolated += factor * accessor(idxs);

		}

		return interpolated;
	}

public:

	PiecewiseLinear(const Grid<dim> &grid, const Vector &vector)
			: grid(&grid), vector(&vector) {
	}

	// TODO: Figure out what to do with rvalue arguments

	PiecewiseLinear(const PiecewiseLinear &that) : grid(that.grid),
			vector(that.vector) {
	}

	PiecewiseLinear &operator=(const PiecewiseLinear &that) {
		grid = that.grid;
		vector = that.vector;
		return *this;
	}

	// TODO: The implementation below breaks for dim > N, where N is some
	//       positive integer. Is there a way to generalize this to a return
	//       type of Function<dim>?

	static_assert(dim > 4, "Dimensions >= 5 not implemented");

	operator Function1() {
		return [this] (double x) { return get(x); };
	}

	operator Function2() {
		return [this] (double x, double y) { return get(x, y); };
	}

	operator Function3() {
		return [this] (double x, double y, double z) {
				return get(x, y, z); };
	}

	operator Function4() {
		return [this] (double x1, double x2, double x3, double x4) {
				return get(x1, x2, x3, x4); };
	}

};

template <Index dim> template <bool isConst>
Real Grid<dim>::GridVector<isConst>::operator()(const Real *coordinates) const {
	Function<dim> lin = static_cast< Function<dim> > (
			PiecewiseLinear<dim>(*grid, *vector) );

	return unpackAndCall(lin, coordinates, GenerateSequence<dim>());
}

}

#endif

