#ifndef QUANT_PDE_DOMAIN
#define QUANT_PDE_DOMAIN

#include <array>       // std::array
#include <memory>      // std::unique_ptr
#include <type_traits> // std::conditional
#include <utility>     // std::forward, std::move

namespace QuantPDE {

template <Index Dimension> class Refiner;

/**
 * A pure virtual class representing a set of points in \f$\mathbb{R}^n\f$
 * together with an order (by which these points can be indexed and iterated
 * over).
 * @tparam Dimension \f$n\f$.
 */
template <Index Dimension>
class Domain {

	static_assert(Dimension > 0, "Dimension must be positive");

	template <bool IsConst>
	class VectorPosition final {

		typedef typename std::conditional<IsConst, const Vector, Vector>
				::type V;
		typedef typename std::conditional<IsConst, Real, Real &>::type R
				;

		const Domain *domain;
		V *vector;
		Index index;

	public:

		VectorPosition(const Domain &domain, V &vector, Index index)
				noexcept : domain(&domain), vector(&vector),
				index(index) {
			QUANT_PDE_ASSERT_LVALUE(domain);
			QUANT_PDE_ASSERT_LVALUE(vector);
		}

		VectorPosition(const VectorPosition &that) noexcept :
				domain(that.domain), vector(that.vector),
				index(that.index) {
		}

		VectorPosition &operator=(const VectorPosition &that) & noexcept
				{
			domain = that.domain;
			vector = that.vector;
			index = that.index;
			return *this;
		}

		////////////////////////////////////////////////////////////////

		R operator*() const {
			return (*vector)(index);
		}

		std::array<Real, Dimension> operator&() const {
			return domain->coordinates(index);
		}

	};

	template <bool IsConst>
	class VectorIterator final {

		typedef typename std::conditional<IsConst, const Vector, Vector>
				::type V;

		const Domain *domain;
		V *vector;
		Index index;

	public:

		VectorIterator(const Domain &domain, V &vector, Index index = 0)
				noexcept : domain(&domain), vector(&vector),
				index(index) {
			QUANT_PDE_ASSERT_LVALUE(domain);
			QUANT_PDE_ASSERT_LVALUE(vector);
		}

		VectorIterator(const VectorIterator &that) noexcept
				: domain(that.domain), vector(that.vector),
				index(that.index) {
		}

		VectorIterator &operator=(const VectorIterator &that) & noexcept
				{
			domain = that.domain;
			vector = that.vector;
			index = that.index;
			return *this;
		}

		////////////////////////////////////////////////////////////////

		bool operator==(const VectorIterator &that) const {
			return domain == that.domain && vector == that.vector &&
					index == that.index;
		}

		bool operator!=(const VectorIterator &that) const {
			return !(*this == that);
		}

		VectorPosition<IsConst> operator*() const {
			return VectorPosition<IsConst>(*domain, *vector, index);
		}

		VectorIterator &operator++() {
			index++;
			return *this;
		}

		VectorIterator operator++(int) {
			const VectorIterator old(*this);
			++(*this);
			return old;
		}

	};

	template <bool IsConst>
	class VectorAccessor final {

		typedef typename std::conditional<IsConst, const Vector, Vector>
				::type V;

		const Domain *domain;
		V *vector;

	public:

		VectorAccessor(const Domain &domain, V &vector)
				noexcept : domain(&domain), vector(&vector) {
			QUANT_PDE_ASSERT_LVALUE(domain);
			QUANT_PDE_ASSERT_LVALUE(vector);
		}

		VectorAccessor(const VectorAccessor &that) noexcept
				: domain(that.domain), vector(that.vector) {
		}

		VectorAccessor &operator=(const VectorAccessor &that) & noexcept
				{
			domain = that.domain;
			vector = that.vector;
			return *this;
		}

		////////////////////////////////////////////////////////////////

		template <typename ...Ts>
		Real operator()(Ts &&...coordinates) const {
			return domain->value(*vector, {{coordinates...}});
		}

		////////////////////////////////////////////////////////////////

		VectorIterator<IsConst> begin() const {
			return VectorIterator<IsConst>(*domain, *vector);
		}

		VectorIterator<IsConst> end() const {
			return VectorIterator<IsConst>(*domain, *vector,
					domain->size());
		}

		template <bool B>
		friend std::ostream &operator<<(std::ostream &os,
				const VectorAccessor<B> &accessor) {
			std::array<Real, Dimension> coordinates;
			for(auto v_n : accessor) {
				coordinates = &v_n;

				//  TODO: Explicit loop unroll
				for(Index i = 0; i < Dimension; i++) {
					os << coordinates[i] << '\t';
				}

				os << *v_n << std::endl;
			}
			return os;
		}
	};

	typedef VectorAccessor<true> VectorAccessorConst;
	typedef VectorAccessor<false> VectorAccessorNonConst;

	virtual Real value(const Vector &vector,
			const std::array<Real, Dimension> &coordinates) const
			= 0;

public:

	virtual ~Domain() {
	}

	// TODO: Document

	VectorAccessorConst accessor(const Vector &vector) const {
		return VectorAccessorConst(*this, vector);
	}

	VectorAccessorNonConst accessor(Vector &&vector) const = delete;

	VectorAccessorNonConst accessor(Vector &vector) const {
		return VectorAccessorNonConst(*this, vector);
	}

	/**
	 * @return A square matrix of order equal to the number of nodes on this
	 *         grid.
	 */
	Matrix matrix() const {
		return Matrix(size(), size());
	}

	/**
	 * @return A zero vector of length equal to the number of nodes on this
	 *         grid.
	 */
	Vector zero() const {
		return Vector::Zero(size());
	}

	/**
	 * @return A vector of ones of length equal to the number of nodes on
	 *         this grid.
	 */
	Vector ones() const {
		return Vector::Ones(size());
	}

	/**
	 * @return A vector of length equal to the number of nodes on this grid.
	 *         No guarantees are made on the content of this vector.
	 */
	Vector vector() const {
		return Vector(size());
	}

	/**
	 * @param function A function.
	 * @return The image of a function on this grid as a vector.
	 */
	Vector image(const Function<Dimension> &function) const {
		Vector v = vector();
		for(auto node : accessor(v)) {
			*node = unpackAndCall<Dimension>( function,
					(&node).data() );
		}
		return v;
	}

	////////////////////////////////////////////////////////////////////////

	/**
	 * Refines the domain in-place.
	 * @param refiner
	 */
	virtual void refine(const Refiner<Dimension> &refiner) = 0;

	/**
	 * @param index An index.
	 * @return The coordinates associated with this index.
	 */
	virtual std::array<Real, Dimension> coordinates(Index index) const = 0;

	/**
	 * @return The total number of nodes on this grid.
	 */
	virtual Index size() const = 0;

};

////////////////////////////////////////////////////////////////////////////////

/**
 * A pure virtual class representing a routine used to refine a domain.
 */
template <Index Dimension>
class Refiner {

public:

	/**
	 * Constructor
	 */
	Refiner() {
	}

	Refiner(const Refiner &that) = delete;
	Refiner &operator=(const Refiner &that) = delete;

	/**
	 * Destructor.
	 */
	virtual ~Refiner() {
	}

	/**
	 * Creates a refined domain from a given one.
	 * @return A refined domain.
	 */
	virtual std::unique_ptr< Domain<Dimension> > refine(
			const Domain<Dimension> &domain) const = 0;

};

#define QUANT_PDE_DOMAIN_REFINE(DERIVED_CLASS)                             \
		do {                                                       \
			std::unique_ptr< Domain<Dimension> > refined =     \
					refiner.refine( *this );           \
			const DERIVED_CLASS &derived = dynamic_cast<       \
					const DERIVED_CLASS &>( *refined); \
			this->operator=( std::move(derived) );             \
		} while(0)

////////////////////////////////////////////////////////////////////////////////

template <Index Dimension>
class RectilinearGrid : public Domain<Dimension> {

	// vsize = -1 so an error is thrown if used without a prior call to
	// initialize()
	RectilinearGrid() noexcept : vsize(-1) {
	}

	void initialize() {
		vsize = 1;

		// TODO: Explicit loop unroll
		for(Index i = 0; i < Dimension; i++) {
			assert(this->axes[i].size() > 0);
			vsize *= this->axes[i].size();
		}
	}

	////////////////////////////////////////////////////////////////////////

	std::array<Index, Dimension> indices(Index index) const {
		std::array<Index, Dimension> array;

		Index m = 1;

		// TODO: Explicit loop unroll
		for(Index i = 0; i < Dimension - 1; i++) {
			m *= axes[i].size();
		}

		// TODO: Explicit loop unroll
		for(Index i = Dimension - 1; i > 0; i--) {
			array[i] = index / m;
			index -= array[i] * m;
			m /= axes[i].size();
		}
		array[0] = index;

		return array;
	}

	////////////////////////////////////////////////////////////////////////

	virtual Real value(const Vector &vector,
			const std::array<Real, Dimension> &coordinates) const;

	////////////////////////////////////////////////////////////////////////

	Index vsize;
	Axis axes[Dimension]; // Private constructor Axis()

	////////////////////////////////////////////////////////////////////////

public:

	/**
	 * Refines a rectilinear grid by refining each axis as follows: for each
	 * two adjacent (distinct) ticks on the axis, a new tick is inserted
	 * halfway between the two.
	 */
	class NewTickBetweenEachPair : public Refiner<Dimension> {

		virtual std::unique_ptr< Domain<Dimension> > refine(
				const Domain<Dimension> &domain) const {
			std::unique_ptr<Domain<Dimension>> pointer(
					new RectilinearGrid());

			const RectilinearGrid<Dimension> &original
					= dynamic_cast<
					const RectilinearGrid<Dimension> &>(
					domain);

			RectilinearGrid<Dimension> &refined = dynamic_cast<
					RectilinearGrid<Dimension> &>(
					*pointer);

			// TODO: Explicit loop unroll
			for(Index k = 0; k < Dimension; k++) {
				// Refine k-th axis

				const Axis &n = original[k];
				Axis &m = refined.axes[k];

				m = Axis( n.size() * 2 - 1 );

				m[0] = n[0];
				Index i = 1, j = 1;
				while(i < n.size()) {
					m[j++] = ( n[i-1] + n[i] ) / 2.;
					m[j++] = n[i++];
				}
			}

			refined.initialize();

			return pointer;
		}

	};

	/**
	 * Constructor.
	 */
	template <typename ...Ts>
	RectilinearGrid(Ts &&...axes) noexcept
			: axes { std::forward<Ts>(axes)... } {
		static_assert(Dimension == sizeof...(Ts),
				"The number of arguments must be consistent "
				"with the dimensions");
		initialize();
	}

	/**
	 * Copy constructor.
	 */
	RectilinearGrid(const RectilinearGrid &that) noexcept
			: vsize(that.vsize) {
		// TODO: Explicit loop unroll
		for(Index i = 0; i < Dimension; i++) {
			axes[i] = that.axes[i];
		};
	}
	/**
	 * Move constructor.
	 */
	RectilinearGrid(RectilinearGrid &&that) noexcept : vsize(that.vsize) {
		// TODO: Explicit loop unroll
		for(Index i = 0; i < Dimension; i++) {
			axes[i] = std::move(that.axis[i]);
		};
	}

	/**
	 * Copy assignment operator.
	 */
	RectilinearGrid &operator=(const RectilinearGrid &that) & noexcept {
		vsize = that.vsize;

		// TODO: Explicit loop unroll
		for(Index i = 0; i < Dimension; i++) {
			axes[i] = that.axes[i];
		};

		return *this;
	}

	/**
	 * Move assignment operator.
	 */
	RectilinearGrid &operator=(RectilinearGrid &&that) & noexcept {
		vsize = that.vsize;

		// TODO: Explicit loop unroll
		for(Index i = 0; i < Dimension; i++) {
			axes[i] = std::move(that.axes[i]);
		};

		return *this;
	}

	////////////////////////////////////////////////////////////////////////

	// TODO: Matrix builders

	/**
	 * Transforms axes-indices to an index corresponding to the (natural)
	 * order imposed by this grid.
	 * @return The index.
	 */
	template <typename ...Ts>
	Index index(Ts ...indices) const {
		static_assert(Dimension == sizeof...(Ts),
				"The number of arguments must be consistent "
				"with the dimensions");

		Index idxs[] {indices...};

		Index
			index = 0,
			horner = 1
		;

		index = idxs[0];

		// TODO: Explicit loop unroll
		for(Index i = 0; i < Dimension - 1; i++) {
			horner *= axes[i].size();
			index += horner * idxs[i + 1];
		};

		return index;
	}

	/**
	 * @param index An index.
	 * @return The index-th axis.
	 */
	const Axis &operator[](Index index) const {
		return axes[index];
	}

	/**
	 * Prettifies and prints the grid.
	 * @param os The output stream.
	 * @param grid A grid.
	 */
	friend std::ostream &operator<<(std::ostream &os,
			const RectilinearGrid<Dimension> &grid) {
		os << grid[0];
		for(Index i = 1; i < Dimension; i++) {
			os << " x " << grid[i];
		};
		return os;
	}

	////////////////////////////////////////////////////////////////////////

	virtual void refine(const Refiner<Dimension> &refiner) {
		QUANT_PDE_DOMAIN_REFINE(RectilinearGrid);
	}

	virtual std::array<Real, Dimension> coordinates(Index index) const {
		std::array<Real, Dimension> array;

		std::array<Index, Dimension> tmp = indices(index);

		// TODO: Explicit loop unroll
		for(Index i = 0; i < Dimension; i++) {
			array[i] = axes[i][tmp[i]];
		};

		return array;
	}

	virtual Index size() const {
		return vsize;
	}

};

typedef RectilinearGrid<1> RectilinearGrid1;
typedef RectilinearGrid<2> RectilinearGrid2;
typedef RectilinearGrid<3> RectilinearGrid3;

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
class PiecewiseLinear {

	static_assert(Dimension > 0, "Dimension must be positive");

	const RectilinearGrid<Dimension> *grid;
	const Vector *vector;

	/*
	Real interpolate(Index *indices, const Real *weights, Index shift,
			Index n = 0) const {
		// Base case
		if(n == Dimension) {
			return grid->accessor(*vector)(indices);
		}

		Index *stride = indices + shift;
		std::memcpy(stride, indices, sizeof(Index) * Dimension);
		stride[n]++;

		return weights[n] * interpolate(indices, weights, shift / 2,
				n + 1) + (1 - weights[n]) * interpolate(stride,
				weights, shift / 2, n + 1);
	}
	*/

public:

	/**
	 * Constructor.
	 */
	PiecewiseLinear(const RectilinearGrid<Dimension> &grid,
			const Vector &vector) noexcept : grid(&grid),
			vector(&vector) {
		QUANT_PDE_ASSERT_LVALUE(grid);
		QUANT_PDE_ASSERT_LVALUE(vector);
	}

	/**
	 * Copy constructor.
	 */
	PiecewiseLinear(const PiecewiseLinear &that) noexcept : grid(that.grid),
			vector(that.vector) {
	}

	/**
	 * Assignment operator.
	 */
	PiecewiseLinear &operator=(const PiecewiseLinear &that) & noexcept {
		grid = that.grid;
		vector = that.vector;
		return *this;
	}

	////////////////////////////////////////////////////////////////////////

	/**
	 * Performs linear interpolation to query the value at the specified
	 * coordinates.
	 * @param coordinates The coordinates.
	 */
	template <typename ...Args>
	Real operator()(Args ...coordinates) const {
		static_assert(Dimension == sizeof...(Args),
				"The number of arguments must be consistent "
				"with the dimensions");

		Real coords[] {coordinates...};

		typedef IntegerPower<2, Dimension> dimpow;
		static_assert(!dimpow::overflow, "Overflow detected");

		Real weights[Dimension];
		Index indices[dimpow::value];

		// For the i-th coordinate, find the ticks on the i-th axis that
		// it lies between along with the distance from the leftmost
		// tick
		// TODO: Explicit loop unroll
		for(Index i = 0; i < Dimension; i++) {
			const Axis &x = (*grid)[i];
			Index length = x.size();

			if(coords[i] <= x[0]) {
				indices[i] = 0;
				weights[i] = 1.;
				continue;
			}

			if(coords[i] >= x[length - 1]) {
				indices[i] = length - 2;
				weights[i] = 0.;
				continue;
			}

			// Binary search to find tick
			Index lo = 0, hi = length - 2, mid = 0;
			Real weight = 0.;
			while(lo <= hi) {
				mid = (lo + hi) / 2;
				if(coords[i] < x[mid]) {
					hi = mid - 1;
				} else if(coords[i] >= x[mid + 1]) {
					lo = mid + 1;
				} else {
					weight = ( x[mid + 1] - coords[i] )
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

		// TODO: Explicit loop unroll
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
			//Index k = unpackAndCallMethod<Dimension>(*grid,
			//		(methodType) &RectilinearGrid<Dimension>
			//		::index, idxs);

			// n-Dimensional
			NaryMethodConst<Index, RectilinearGrid<Dimension>,
					Dimension, Index> tmp =
					&RectilinearGrid<Dimension>::index;
			Index k = unpackAndCall<Dimension>(*grid, tmp, idxs);

			interpolated += factor * (*vector)(k);
		};

		return interpolated;
	}

};

typedef PiecewiseLinear<1> PiecewiseLinear1;
typedef PiecewiseLinear<2> PiecewiseLinear2;
typedef PiecewiseLinear<3> PiecewiseLinear3;

template <Index Dimension>
Real RectilinearGrid<Dimension>::value(const Vector &vector,
			const std::array<Real, Dimension> &coordinates) const {
	return unpackAndCall<Dimension>( PiecewiseLinear<Dimension>(*this,
			vector), coordinates.data() );
}

}

#endif
