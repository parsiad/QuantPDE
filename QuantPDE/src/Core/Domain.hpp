#ifndef QUANT_PDE_CORE_DOMAIN
#define QUANT_PDE_CORE_DOMAIN

#include <array>       // std::array
#include <cassert>     // assert
#include <cstdlib>     // size_t
#include <iostream>    // std::ostream
#include <memory>      // std::shared_ptr, std::unique_ptr
#include <type_traits> // std::conditional
#include <utility>     // std::forward, std::move

namespace QuantPDE {

template <Index Dimension> class Refiner;

/**
 * A pure virtual class representing a (finite) set of points in some space.
 */
class DomainBase {

public:

	virtual ~DomainBase() {
	}

	/**
	 * @return A square matrix of order equal to the number of nodes on this
	 *         domain.
	 */
	Matrix matrix() const {
		return Matrix(size(), size());
	}

	/**
	 * @return An identity matrix of order equal to the number of nodes on
	 *         this domain.
	 */
	Matrix identity() const {
		Matrix M(size(), size());
		M.setIdentity();
		return M;
	}

	/**
	 * @return A zero vector of length equal to the number of nodes on this
	 *         domain.
	 */
	Vector zero() const {
		return Vector::Zero(size());
	}

	/**
	 * @return A vector of ones of length equal to the number of nodes on
	 *         this domain.
	 */
	Vector ones() const {
		return Vector::Ones(size());
	}

	/**
	 * @return A vector of length equal to the number of nodes on this
	 *         domain. No guarantees are made on the contents of this
	 *         vector.
	 */
	Vector vector() const {
		return Vector(size());
	}

	/**
	 * @return The total number of nodes on this domain.
	 */
	virtual Index size() const = 0;

};

/**
 * A pure virtual class representing a (finite) set of points in
 * \f$R^n\f$ together with an order (by which these points can be
 * indexed and iterated over).
 * @tparam Dimension \f$n\f$.
 */
template <Index Dimension>
class Domain : public DomainBase {

	static_assert(Dimension > 0, "Dimension must be positive");

	////////////////////////////////////////////////////////////////////////

	enum class Ownership { CONST, NON_CONST, SHARED };

	template <Ownership O>
	using V0 = typename std::conditional<
		O == Ownership::CONST,
		const Vector *,                    // CONST
		typename std::conditional<
			O == Ownership::NON_CONST,
			Vector *,                  // NON_CONST
			std::shared_ptr<Vector>    // SHARED
		>::type
	>::type;

	template <Ownership O>
	class VectorPosition final {

		typedef typename std::conditional<
			O == Ownership::CONST,
			Real,
			Real &
		>::type R;

		const Domain *domain;
		V0<O> vector;
		Index index;

	public:

		template <typename D, typename V>
		VectorPosition(D &domain, V &&vector, Index index) noexcept
				: domain(&domain),
				vector( std::forward<V>(vector) ), index(index)
				{
		}

		VectorPosition(const VectorPosition &that) noexcept
				: domain(that.domain), vector(that.vector),
				index(that.index) {
		}

		VectorPosition(VectorPosition &&that) noexcept
				: domain(that.domain),
				vector( std::move(that.vector) ),
				index(that.index) {
		}

		VectorPosition &operator=(const VectorPosition &that) & noexcept
				{
			domain = that.domain;
			vector = that.vector;
			index = that.index;
			return *this;
		}

		VectorPosition &operator=(VectorPosition &&that) & noexcept {
			domain = that.domain;
			vector = std::move(that.vector);
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

	template <Ownership O>
	class VectorIterator final {

		const Domain *domain;
		V0<O> vector;
		Index index;

	public:

		template <typename D, typename V>
		VectorIterator(D &domain, V &&vector, Index index = 0) noexcept
				: domain(&domain),
				vector( std::forward<V>(vector) ), index(index)
				{
		}

		VectorIterator(const VectorIterator &that) noexcept
				: domain(that.domain), vector(that.vector),
				index(that.index) {
		}

		VectorIterator(const VectorIterator &&that) noexcept
				: domain(that.domain),
				vector( std::move(that.vector) ),
				index(that.index) {
		}

		VectorIterator &operator=(const VectorIterator &that) & noexcept
				{
			domain = that.domain;
			vector = that.vector;
			index = that.index;
			return *this;
		}

		VectorIterator &operator=(VectorIterator &&that) & noexcept {
			domain = that.domain;
			vector = std::move(that.vector);
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

		VectorPosition<O> operator*() const {
			return VectorPosition<O>(*domain, vector, index);
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

	template <Ownership O>
	class VectorAccessor final {

		const Domain *domain;
		V0<O> vector;

	public:

		template <typename D, typename V>
		VectorAccessor(D &domain, V &&vector) noexcept
				: domain(&domain),
				vector( std::forward<V>(vector) ) {
		}

		VectorAccessor(const VectorAccessor &that) noexcept
				: domain(that.domain), vector(that.vector) {
		}

		VectorAccessor(VectorAccessor &&that) noexcept
				: domain(that.domain),
				vector( std::move(that.vector) ) {
		}

		VectorAccessor &operator=(const VectorAccessor &that) & noexcept
				{
			domain = that.domain;
			vector = that.vector;
			return *this;
		}

		VectorAccessor &operator=(VectorAccessor &&that) & noexcept {
			domain = that.domain;
			vector = std::move(that.vector);
			return *this;
		}

		////////////////////////////////////////////////////////////////

		template <typename ...Ts>
		Real operator()(Ts &&...coordinates) const {
			return domain->value(*vector, {{coordinates...}});
		}

		////////////////////////////////////////////////////////////////

		VectorIterator<O> begin() const {
			return VectorIterator<O>(*domain, vector);
		}

		VectorIterator<O> end() const {
			return VectorIterator<O>(*domain, vector,
					domain->size());
		}

		////////////////////////////////////////////////////////////////

		template <Ownership O1>
		friend std::ostream &operator<<(std::ostream &os,
				const VectorAccessor<O1> &accessor) {
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

	////////////////////////////////////////////////////////////////////////

	virtual Real value(const Vector &vector,
			const std::array<Real, Dimension> &coordinates) const
			= 0;

public:

	/**
	 * Accessors can be used to iterate over vectors with respect to an
	 * order imposed by a specific domain.
	 *
	 * Consider the following example evaluating the function
	 * \f$f\left(x,y\right)=xy\f$ on two-dimensional domain:
	 * \code{.cpp}
	 * // Initialize a domain (assume My2DDomain is a subtype of Domain2)
	 * My2DDomain D;
	 *
	 * // Create a vector
	 * Vector v = D.vector();
	 *
	 * // For-each loop over the elements in the vector
	 * for(auto v_xy : D.accessor(v)) {
	 * 	// Get the coordinates associated with this node
	 * 	auto xy = &v_xy; // xy is of type std::array<Real, 2>
	 * 	Real x = xy[0], y = xy[1];
	 *
	 * 	// Set the value at this node
	 * 	*v_xy = x * y;
	 * }
	 * \endcode
	 *
	 * This becomes particularly useful if we would like to evaluate an
	 * n-dimensional function on a n-dimensional domain, but we only know
	 * the dimension of the domain at compile-time (as a template argument).
	 *
	 * Consider the following example evaulating the function
	 * \f$f\left(x_1, \ldots, x_n\right)=\prod_{i=1}^n x_i\f$ on an
	 * n-dimensional domain:
	 * \code{.cpp}
	 * // Initialize a domain (assume MyDomain<Dimension> is a subtype of
	 * // Domain<Dimension>)
	 * MyDomain<Dimension> D;
	 *
	 * Vector v = D.vector();
	 * for(auto v_x : D.accessor(v)) {
	 * 	auto x = &v_x; // x is of type std::array<Real, Dimension>
	 * 	Real product = 1.;
	 * 	for(Index i = 0; i < Dimension; i++) {
	 * 		product *= x[i];
	 * 	}
	 * 	*v = product;
	 * }
	 * \endcode
	 * Since the dimension is a template parameter, it is safe to assume
	 * that the innermost loop will be unrolled by any (smart) compiler for
	 * low dimensions, and that autovectorizers will be able to pick up on
	 * potential optimizations.
	 * @return An accessor to a vector.
	 */
	VectorAccessor<Ownership::CONST> accessor(const Vector &vector) const {
		return VectorAccessor<Ownership::CONST>(*this, &vector);
	}

	VectorAccessor<Ownership::SHARED> accessor(Vector &&vector) const {
		return VectorAccessor<Ownership::SHARED>(*this, std::shared_ptr<
				Vector>( new Vector(std::move(vector)) ));
	};

	VectorAccessor<Ownership::NON_CONST> accessor(Vector &vector) const {
		return VectorAccessor<Ownership::NON_CONST>(*this, &vector);
	}

	/**
	 * Used to evaluate the function on the domain.
	 * Consider the following example, evaluating the function
	 * \f$f\left(x,y\right)=xy\f$ on a two-dimensional domain:
	 * \code{.cpp}
	 * // Initialize a domain (assume My2DDomain is a subtype of Domain2)
	 * My2DDomain D;
	 *
	 * // Create a vector with the results
	 * Vector v = D.image( [] (Real x, Real y) { return x * y; } );
	 * \endcode
	 * @param function A function.
	 * @return The image of a function on this domain as a vector.
	 * @see QuantPDE::Domain::accessor
	 */
	template <typename F>
	Vector image(F &&function) const {
		Vector v = vector();
		for(auto node : accessor(v)) {
			*node = packAndCall<Dimension>(
				std::forward<F>(function),
				(&node).data()
			);
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

};

typedef Domain<1> Domain1;
typedef Domain<2> Domain2;
typedef Domain<3> Domain3;

////////////////////////////////////////////////////////////////////////////////

/**
 * A pure virtual class representing a routine used to refine a domain.
 */
template <Index Dimension>
class Refiner {

public:

	/**
	 * Constructor.
	 */
	Refiner() noexcept {
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

#define QUANT_PDE_REFINE_IN_PLACE(DERIVED_CLASS, REFINER)            \
		do {                                                 \
			(*this) = dynamic_cast<DERIVED_CLASS &&>(    \
					*REFINER.refine( *this ) );  \
		} while(0)

////////////////////////////////////////////////////////////////////////////////

/**
 * A (finite) rectilinear grid.
 * @see QuantPDE::Domain
 */
template <Index Dimension>
class RectilinearGrid : public Domain<Dimension> {

	template <bool IsConst>
	class VectorAxesIndexer final {

		typedef typename std::conditional<
			IsConst,
			const Vector,
			Vector
		>::type V;

		typedef typename std::conditional<
			IsConst,
			Real,
			Real &
		>::type R;

		const RectilinearGrid *grid;
		V *vector;

	public:

		template <typename R, typename W>
		VectorAxesIndexer(R &grid, W &vector) noexcept : grid(&grid),
				vector(&vector) {
		}

		VectorAxesIndexer(const VectorAxesIndexer &that) noexcept
				: grid(that.grid), vector(that.vector) {
		}

		VectorAxesIndexer &operator=(const VectorAxesIndexer &that)
				& noexcept {
			grid = that.grid;
			vector = that.vector;
		}

		///////////////////////////////////////////////////////////////

		template <typename ...Ts>
		R operator()(Ts ...indices) const {
			static_assert(Dimension == sizeof...(Ts),
					"The number of arguments must be "
					"consistent with the dimensions");
			return (*vector)(grid->index(indices...));
		}

	};

	typedef VectorAxesIndexer<true> VectorAxesIndexerConst;
	typedef VectorAxesIndexer<false> VectorAxesIndexerNonConst;

	class VectorAxesIndexerValue final {

		const RectilinearGrid *grid;
		Vector vector;

	public:

		template <typename R, typename V>
		VectorAxesIndexerValue(R &grid, V &&vector) noexcept
				: grid(&grid), vector( std::forward<V>(vector) )
				{
		}

		// Only moves are allowed
		VectorAxesIndexerValue(const VectorAxesIndexerValue &)
				= delete;
		VectorAxesIndexerValue &operator=(
				const VectorAxesIndexerValue &) = delete;

		VectorAxesIndexerValue(VectorAxesIndexerValue &&that) noexcept
				: grid(that.grid),
				vector( std::move(that.vector) ) {
		}

		VectorAxesIndexerValue &operator=(VectorAxesIndexerValue &&that)
				& noexcept {
			grid = that.grid;
			vector = std::move(that.vector);
		}

		///////////////////////////////////////////////////////////////

		template <typename ...Ts>
		Real operator()(Ts ...indices) const {
			static_assert(Dimension == sizeof...(Ts),
					"The number of arguments must be "
					"consistent with the dimensions");
			return vector(grid->index(indices...));
		}
	};

	////////////////////////////////////////////////////////////////////////

	class MatrixAxesIndexer final {

		const RectilinearGrid *grid;
		Matrix *matrix;

	public:

		template <typename R>
		MatrixAxesIndexer(R &grid, Matrix &matrix) noexcept
				: grid(&grid), matrix(&matrix) {
			assert(matrix.rows() == grid.size());
			assert(matrix.cols() == grid.size());
		}

		MatrixAxesIndexer(const MatrixAxesIndexer &that) noexcept
				: grid(that.grid), matrix(that.matrix) {
		}

		MatrixAxesIndexer &operator=(const MatrixAxesIndexer &that) &
				noexcept {
			grid = that.grid;
			matrix = that.matrix;
			return *this;
		}

		////////////////////////////////////////////////////////////////

		template <typename ...Ts>
		Real &operator()(Ts ...indices) const {
			static_assert((Dimension * 2) == sizeof...(Ts),
					"The number of arguments must be "
					"consistent with the dimensions");

			// TODO: Currently, we take the argument pack and unpack
			//       it into an array. Each half of the array is
			//       then packed and the member function is called.
			//       Instead, we should split the pack directly.

			Index idxs[] {indices...};

			NaryMethodConst<Index, RectilinearGrid, Dimension,
					Index> pointer
					= &RectilinearGrid::index;

			Index i = packAndCall<Dimension>(*grid, pointer, idxs);
			Index j = packAndCall<Dimension>(*grid, pointer,
					idxs + Dimension);


			return matrix->insert(i, j);
		}

	};

	////////////////////////////////////////////////////////////////////////

	class LazyMatrixBuilder final {

		RectilinearGrid *grid;
		std::vector<Entry> entries;

	public:

		LazyMatrixBuilder(const RectilinearGrid &grid) noexcept
				: grid(&grid) {
		}

		LazyMatrixBuilder(const RectilinearGrid &grid, size_t nonzeros)
				noexcept : grid(&grid) {
			entries.reserve(nonzeros);
		}

		LazyMatrixBuilder(const LazyMatrixBuilder &that) noexcept :
				grid(that.grid), entries(entries) {
		}

		LazyMatrixBuilder(const LazyMatrixBuilder &&that) noexcept :
				grid(that.grid),
				entries( std::move(that.entries) ) {
		}

		LazyMatrixBuilder &operator=(const LazyMatrixBuilder &that) &
				noexcept {
			grid = that.grid;
			entries = that.entries;
		}

		LazyMatrixBuilder &operator=(LazyMatrixBuilder &&that) &
				noexcept {
			grid = that.grid;
			entries = std::move(entries);
		}

		////////////////////////////////////////////////////////////////

		template <typename ...Ts>
		Real &operator()(Ts ...indices) {
			static_assert((Dimension * 2) == sizeof...(Ts),
					"The number of arguments must be "
					"consistent with the dimensions");

			// TODO: Currently, we take the argument pack and unpack
			//       it into an array. Each half of the array is
			//       then packed and the member function is called.
			//       Instead, we should split the pack directly.

			Index idxs[] {indices...};

			NaryMethodConst<Index, RectilinearGrid, Dimension,
					Index> pointer
					= &RectilinearGrid::index;

			Index i = packAndCall<Dimension>(*grid, pointer, idxs);
			Index j = packAndCall<Dimension>(*grid, pointer,
					idxs + Dimension);

			entries.push_back( Entry(i, j, 0.) );
			return entries.back().value();
		}

		Matrix matrix() {
			Matrix M(grid->size(), grid->size());
			M.setFromTriplets(entries.begin(), entries.end());
			return M;
		}

	};

	class FastMatrixBuilder final {

		const RectilinearGrid *grid;
		Matrix M;

	public:

		template <typename S>
		FastMatrixBuilder(const RectilinearGrid &grid,
				S &&reserveSizes) noexcept : grid(&grid),
				M(grid.size(), grid.size()) {
			M.reserve( std::forward<S>(reserveSizes) );
		}

		FastMatrixBuilder(const FastMatrixBuilder &that) noexcept
				: grid(that.grid), M(that.M) {
		}

		FastMatrixBuilder(FastMatrixBuilder &&that) noexcept
				: grid(that.grid), M( std::move(that.M) ) {
		}

		FastMatrixBuilder &operator=(const FastMatrixBuilder &that) &
				noexcept {
			grid = that.grid;
			M = that.M;
		}

		FastMatrixBuilder &operator=(FastMatrixBuilder &&that) &
				noexcept {
			grid = that.grid;
			M = std::move(that.M);
		}

		////////////////////////////////////////////////////////////////

		template <typename ...Ts>
		Real &operator()(Ts ...indices) {
			static_assert((Dimension * 2) == sizeof...(Ts),
					"The number of arguments must be "
					"consistent with the dimensions");

			// TODO: Currently, we take the argument pack and unpack
			//       it into an array. Each half of the array is
			//       then packed and the member function is called.
			//       Instead, we should split the pack directly.

			Index idxs[] {indices...};

			NaryMethodConst<Index, RectilinearGrid, Dimension,
					Index> pointer
					= &RectilinearGrid::index;

			Index i = packAndCall<Dimension>(*grid, pointer, idxs);
			Index j = packAndCall<Dimension>(*grid, pointer,
					idxs + Dimension);

			return M.insert(i, j);
		}

		Matrix &&matrix() {
			M.makeCompressed();

			// This is dangerous! If this method is called twice,
			// we will get a segfault (if we are lucky...)
			return std::move(M);
		}

	};

	////////////////////////////////////////////////////////////////////////

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

	virtual Real value(const Vector &vector,
			const std::array<Real, Dimension> &coordinates) const;

	////////////////////////////////////////////////////////////////////////

	Index vsize;
	Axis axes[Dimension]; // Private constructor Axis()
	                      // (Domain is a friend of Axis)

	////////////////////////////////////////////////////////////////////////

public:

	/**
	 * Refines a rectilinear grid by refining each axis as follows: for each
	 * two adjacent (distinct) ticks on the axis, a new tick is inserted
	 * halfway between the two.
	 */
	class NewTickBetweenEachPair : public Refiner<Dimension> {

	public:

		virtual std::unique_ptr< Domain<Dimension> > refine(
				const Domain<Dimension> &domain) const {
			std::unique_ptr< Domain<Dimension> > pointer(
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

	/**
	 * An indexer can be used to access the elements of a vector on a
	 * rectilinear grid. Indices passed to an index are resolved using the
	 * underlying axes associated with the rectilinear grid.
	 *
	 * Consider the following example evaluating the function
	 * \f$f\left(x,y\right)=xy\f$ on a grid defined by the Cartesian product
	 * \f$\left\{0, 1\right\} \times \left\{-1, 0, 1\right\}\f$:
	 * \code{.cpp}
	 * Axis X { 0., 1. };
	 * Axis Y {-1., 0., 1. };
	 * RectilinearGrid2 G(X, Y);
	 *
	 * // Create a vector
	 * Vector v = G.vector();
	 *
	 * auto indexer = G.indexer(v);
	 * for(Index m = 0; m < X.size(); m++) {
	 * 	for(Index n = 0; n < Y.size(); n++) {
	 * 		indexer(m, n) = X[m] * Y[n];
	 * 	}
	 * }
	 * \endcode
	 *
	 * It should in many situations (such as the one above), an accessor can
	 * be used to achieve the same effect.
	 *
	 * @return An indexer to a vector.
	 * @see QuantPDE::Domain::accessor
	 */
	VectorAxesIndexerNonConst indexer(const Vector &vector) const {
		return VectorAxesIndexerNonConst(*this, vector);
	}

	VectorAxesIndexerValue indexer(Vector &&vector) const {
		return VectorAxesIndexerValue(*this, vector);
	}

	VectorAxesIndexerConst indexer(Vector &vector) const {
		return VectorAxesIndexerConst(*this, vector);
	}

	/**
	 * A matrix indexer can only be used to insert nonzero entries into a
	 * matrix. Indices passed to an index are resolved using the underlying
	 * axes associated with the rectilinear grid.
	 *
	 * Attempting to insert into a location twice will cause undefined
	 * behaviour.
	 *
	 * The example below builds the discrete Laplace operator:
	 * \code{.cpp}
	 * // Create a 2-dimensional grid
	 * RectilinearGrid2 G(X, Y);
	 * Index m = X.size();
	 * Index n = Y.size();
	 *
	 * // Create a matrix with 4 nonzero element reserved per row
	 * Matrix M = G.matrix();
	 * M.reserve(IntegerVector::Constant(G.size(), 4));
	 *
	 * auto indexer = G.indexer(M);
	 *
	 * // Branching inside the loop is expensive and should be eliminated
	 * // in production code
	 * for(Index i = 0; i < X.size(); i++) {
	 * 	for(Index j = 0; j < Y.size(); j++) {
	 * 		if(i > 0)     indexer(i, j, i - 1, j    ) = -1.;
	 * 		if(j > 0)     indexer(i, j, i,     j - 1) = -1.;
	 * 		              indexer(i, j, i,     j    ) =  4.;
	 * 		if(j < n - 1) indexer(i, j, i,     j + 1) = -1.;
	 * 		if(i < m - 1) indexer(i, j, i + 1, j    ) = -1.;
	 * 	}
	 * }
	 *
	 * // The following results in undefined behaviour since we have already
	 * // assigned a value to this position:
	 * // M_G(0, 0, 0, 0) = 4.;
	 *
	 * M.makeCompressed();
	 * \endcode
	 *
	 * @return An indexer to a matrix.
	 */
	MatrixAxesIndexer indexer(Matrix &matrix) const {
		return MatrixAxesIndexer(*this, matrix);
	}

	/**
	 * A matrix builder is used to create a matrix.
	 *
	 * This should be used when the
	 * approximate number of nonzero entries per row (if the underlying
	 * implementation is row-major) is not known ahead of time.
	 *
	 * The example below builds the discrete Laplace operator:
	 * \code{.cpp}
	 * // Create a 1-dimensional grid
	 * Grid2 G(X, Y);
	 * Index m = X.size();
	 * Index n = Y.size();
	 *
	 * auto builder = G.builder();
	 *
	 * // Branching inside the loop is expensive and should be eliminated in
	 * // production code
	 * for(Index i = 0; i < m; i++) {
	 * 	for(Index j = 0; j < n; j++) {
	 * 		if(i > 0)     builder(i, j, i - 1, j    ) = -1.;
	 * 		if(j > 0)     builder(i, j, i,     j - 1) = -1.;
	 * 		              builder(i, j, i,     j    ) =  4.;
	 * 		if(j < n - 1) builder(i, j, i,     j + 1) = -1.;
	 * 		if(i < m - 1) builder(i, j, i + 1, j    ) = -1.;
	 * 	}
	 * }
	 *
	 * Matrix M = builder.matrix();
	 * \endcode
	 *
	 * A call to matrix() should only occur once. The results are undefined
	 * otherwise. Similarly, a particular nonzero entriy should only be
	 * specified once.
	 *
	 * @return A matrix builder.
	 */
	LazyMatrixBuilder builder() const {
		return LazyMatrixBuilder(*this);
	}

	/**
	 * @param nonzeros The approximate total number of nonzero entries in
	 *                 the resulting matrix.
	 */
	LazyMatrixBuilder builder(size_t nonzeros) const {
		return LazyMatrixBuilder(*this, nonzeros);
	}

	/**
	 * When the number of nonzeros per row (if the underlying implementation
	 * is row-major) is known a priori, this method should be used.
	 * @param reserveSizes An integer vector whose i-th element represents
	 *                 the number of elements in the i-th row (if the
	 *                 underlying implementation is row-major) of the
	 *                 matrix.
	 */
	template <typename S>
	FastMatrixBuilder builder(S &&reserveSizes) const {
		return FastMatrixBuilder(*this, std::forward<S>(reserveSizes));
	}

	////////////////////////////////////////////////////////////////////////

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
		QUANT_PDE_REFINE_IN_PLACE(RectilinearGrid, refiner);
	}

	/**
	 * @param index An index corresponding to a node on the grid.
	 * @return An array of indices, each index corresponding to a tick on
	 *         an axis affiliated with the grid.
	 */
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

public:

	/**
	 * Constructor.
	 */
	template <typename R, typename V>
	PiecewiseLinear(R &grid, V &vector) noexcept : grid(&grid),
			vector(&vector) {
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
	template <typename ...Ts>
	Real operator()(Ts ...coordinates) const {
		static_assert(Dimension == sizeof...(Ts),
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
			//Index k = packAndCallMethod<Dimension>(*grid,
			//		(methodType) &RectilinearGrid<Dimension>
			//		::index, idxs);

			// n-Dimensional
			NaryMethodConst<Index, RectilinearGrid<Dimension>,
					Dimension, Index> tmp =
					&RectilinearGrid<Dimension>::index;
			Index k = packAndCall<Dimension>(*grid, tmp, idxs);

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
	return packAndCall<Dimension>( PiecewiseLinear<Dimension>(*this,
			vector), coordinates.data() );
}

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
class PlainMap : public Map<Dimension> {

	const Domain<Dimension> *domain;

public:

	/**
	 * Constructor.
	 */
	template <typename D>
	PlainMap(D &domain) noexcept : domain(&domain) {
	}

	/**
	 * Copy constructor.
	 */
	PlainMap(const PlainMap &that) noexcept : domain(that.domain) {
	}

	/**
	 * Assignment operator.
	 */
	PlainMap &operator=(const PlainMap &that) & noexcept {
		domain = that.domain;
	}

	virtual Vector operator()(const Function<Dimension> &function) const {
		return domain->image(function);
	}

	virtual Vector operator()(Function<Dimension> &&function) const {
		return domain->image( std::move(function) );
	}

};

typedef PlainMap<1> PlainMap1;
typedef PlainMap<2> PlainMap2;
typedef PlainMap<3> PlainMap3;

}

#endif

