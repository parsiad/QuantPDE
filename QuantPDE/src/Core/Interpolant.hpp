#ifndef QUANT_PDE_CORE_INTERPOLANT
#define QUANT_PDE_CORE_INTERPOLANT

#include <array>   // std::array
#include <memory>  // std::unique_ptr
#include <utility> // std::forward, std::move

namespace QuantPDE {

/**
 * A function that interpolates data on domain nodes.
 */
template <Index Dimension>
class Interpolant {

	static_assert(Dimension > 0, "Dimension must be positive");

	/**
	 * Performs interpolation to query the value at the specified
	 * coordinates.
	 * @param coordinates The coordinates.
	 * @return Interpolated value.
	 */
	virtual Real interpolate(const std::array<Real, Dimension> &coordinates)
			const = 0;

	typedef std::unique_ptr<Interpolant> P;

public:

	/**
	 * Destructor.
	 */
	virtual ~Interpolant() {
	}

	/**
	 * Performs interpolation to query the value at the specified
	 * coordinates.
	 * @param coordinates The coordinates.
	 * @return Interpolated value.
	 */
	template <typename ...Ts>
	Real operator()(Ts ...coordinates) const {
		return interpolate( {{coordinates...}} );
	}

	////////////////////////////////////////////////////////////////////////
	// Wrapper
	////////////////////////////////////////////////////////////////////////
	class Wrapper : public Interpolant {
		P p;
		virtual Real interpolate(
				const std::array<Real, Dimension> &) const {
			throw 1; // This should not be visible
		}
	public:
		template <typename ...Ts>
		Real operator()(Ts ...coordinates) const {
			return (*p)(coordinates...);
		}
		QUANT_PDE_CORE_WRAPPER_BODY(P, p)
	};
	////////////////////////////////////////////////////////////////////////

	/**
	 * @return A clone of this interpolant.
	 */
	virtual P clone() const = 0;

};

typedef Interpolant<1> Interpolant1;
typedef Interpolant<2> Interpolant2;
typedef Interpolant<3> Interpolant3;

/**
 * A factory for interpolants.
 */
template <Index Dimension>
class InterpolantFactory {

	typedef std::unique_ptr<Interpolant<Dimension>> I;

public:

	/**
	 * Destructor.
	 */
	virtual ~InterpolantFactory() {
	}

	/**
	 * Creates an interpolant for data points on a domain.
	 * @param vector Data points.
	 * @return An interpolant.
	 */
	virtual I make(const Vector &vector) const = 0;
	virtual I make(Vector &&vector) const = 0;

};

typedef InterpolantFactory<1> InterpolantFactory1;
typedef InterpolantFactory<2> InterpolantFactory2;
typedef InterpolantFactory<3> InterpolantFactory3;

////////////////////////////////////////////////////////////////////////////////

/**
 * A function \f$f\f$ defined piecewise whose pieces are affine functions.
 * Suppose \f$f\f$ is defined on a grid composed of \f$n\f$ axes, and that the
 * \f$i\f$-th axis is a partition consisting of \f$m_i\f$ nodes.
 * Querying \f$f\f$ at a point involves a binary search on each axis.
 * This procedure has worst-case complexity
 * \f$\sim \sum_i \lg m_i \leq n \lg m\f$, where
 * \f$m \equiv \max\left\{m_i\right\}\f$.
 */
template <Index Dimension>
class PiecewiseLinear : public Interpolant<Dimension> {

	typedef std::unique_ptr<Interpolant<Dimension>> I;
	typedef std::unique_ptr<InterpolantFactory<Dimension>> F;

	const RectilinearGrid<Dimension> *grid;
	Vector vector;

	/*
	Real interpolate(Index *indices, const Real *weights, Index shift,
			Index n = 0) const {
		// Base case
		if(n == Dimension) {
			return (grid->indexer(*vector))(indices);
		}

		Index *stride = indices + shift;
		std::memcpy(stride, indices, sizeof(Index) * Dimension);
		stride[n]++;

		return weights[n] * interpolate(indices, weights, shift / 2,
				n + 1) + (1 - weights[n]) * interpolate(stride,
				weights, shift / 2, n + 1);
	}
	*/

	virtual Real interpolate(const std::array<Real, Dimension> &coordinates)
			const {
		typedef IntegerPower<2, Dimension> dimpow;
		static_assert(!dimpow::overflow, "Overflow detected");

		Real weights[Dimension];
		Index indices[dimpow::value];

		// For the i-th coordinate, find the ticks on the i-th axis that
		// it lies between along with the distance from the leftmost
		// tick
		// TODO: Optimize (loop unroll)
		for(Index i = 0; i < Dimension; i++) {

			const Axis &x = (*grid)[i];
			Index length = x.size();

			const Real ci = coordinates[i];

			if(ci <= x[0]) {
				indices[i] = 0;
				weights[i] = 1.;
				continue;
			}

			if(ci >= x[length - 1]) {
				indices[i] = length - 2;
				weights[i] = 0.;
				continue;
			}

			// Binary search to find tick
			Index lo = 0, hi = length - 2, mid = 0;
			Real weight = 0.;
			while(lo <= hi) {
				mid = (lo + hi) / 2;
				if(ci < x[mid]) {
					hi = mid - 1;
				} else if(ci >= x[mid + 1]) {
					lo = mid + 1;
				} else {
					weight = ( x[mid + 1] - ci )
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

		Index idxs[Dimension];
		Real interpolated = 0.;

		// TODO: Optimize (loop unroll)
		for(Index i = 0; i < dimpow::value; i++) {
			Real factor = 1.;

			// Unroll loop
			for(Index j = 0; j < Dimension; j++) {
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

			// 1-Dimensional example
			//typedef Index (RectilinearGrid<Dimension>
			//		::*methodType)(Index) const;
			//Index k = packAndCallMethod<Dimension>(*grid,
			//		(methodType) &RectilinearGrid<Dimension>
			//		::index, idxs);

			// n-Dimensional
			NaryMethodConst<Index, RectilinearGrid<Dimension>,
					Dimension, Index> tmp =
					&RectilinearGrid<Dimension>::index;
			Index k = packAndCall<Dimension>(*grid, tmp, idxs);

			interpolated += factor * vector(k);
		};

		return interpolated;
	}

public:

	class Factory : public InterpolantFactory<Dimension> {

		const RectilinearGrid<Dimension> *grid;

	public:

		/**
		 * Constructor.
		 */
		template <typename G>
		Factory(G &grid) noexcept : grid(&grid) {
		}

		/**
		 * Copy constructor.
		 */
		Factory(const Factory &that) noexcept : grid(that.grid) {
		}

		/**
		 * Assignment operator.
		 */
		Factory &operator=(const Factory &that) noexcept {
			grid = that.grid;
			return *this;
		}

		virtual I make(const Vector &vector) const {
			return I(new PiecewiseLinear(*grid, vector));
		}

		virtual I make(Vector &&vector) const {
			return I(new PiecewiseLinear(*grid, std::move(vector)));
		}

	};

	/**
	 * Constructor.
	 */
	template <typename G, typename V>
	PiecewiseLinear(G &grid, V &&vector) noexcept : grid(&grid),
			vector( std::forward<V>(vector) ) {
	}

	/**
	 * Copy constructor.
	 */
	PiecewiseLinear(const PiecewiseLinear &that) noexcept : grid(that.grid),
			vector(that.vector) {
	}

	/**
	 * Move constructor.
	 */
	PiecewiseLinear(PiecewiseLinear &&that) noexcept : grid(that.grid),
			vector( std::move(that.vector) ) {
	}

	/**
	 * Copy assignment operator.
	 */
	PiecewiseLinear &operator=(const PiecewiseLinear &that) & noexcept {
		grid = that.grid;
		vector = that.vector;
		return *this;
	}

	/**
	 * Move assignment operator.
	 */
	PiecewiseLinear &operator=(PiecewiseLinear &&that) & noexcept {
		grid = that.grid;
		vector = std::move(that.vector);
		return *this;
	}

	virtual I clone() const {
		return I(new PiecewiseLinear(*this));
	}

};

template <Index Dimension>
std::unique_ptr<InterpolantFactory<Dimension>>
		RectilinearGrid<Dimension>::defaultInterpolantFactory() const {
	return std::unique_ptr<InterpolantFactory<Dimension>>(
			new typename PiecewiseLinear<Dimension>::Factory(
			*this));
}

typedef PiecewiseLinear<1> PiecewiseLinear1;
typedef PiecewiseLinear<2> PiecewiseLinear2;
typedef PiecewiseLinear<3> PiecewiseLinear3;

} // QuantPDE

#endif

