#ifndef QUANT_PDE_GRID_HPP
#define QUANT_PDE_GRID_HPP

#include <algorithm>        // std::swap
#include <cstring>          // std::memcpy
#include <cassert>          // assert
#include <initializer_list> // std::initializer_list
#include <iostream>         // std::ostream
#include <utility>          // std::forward

namespace QuantPDE {

/**
 * A rectilinear grid of arbitrary dimension.
 */
template <Index dim>
class Grid final {

	// Check if dimension is positive
	static_assert(dim > 0, "Dimension must be positive");

	template <bool isConst>
	class VectorIndex final {

		typedef typename std::conditional<isConst, const Vector,
				Vector>::type V;
		typedef typename std::conditional<isConst, Real,
				Real &>::type D;

		V *vector;
		Index index;
		Index idxs[dim];

	public:

		VectorIndex(const Grid &grid, V &vector, Index index) noexcept
				: vector(&vector), index(index) {
			grid.roll(index, idxs);
		}

		VectorIndex(const Grid &grid, V &&vector, Index index) = delete;

		VectorIndex(const VectorIndex &that) noexcept
				: vector(that.vector), index(that.index) {
			std::memcpy(idxs, that.idxs, sizeof(Index) * dim);
		}

		VectorIndex &operator=(const VectorIndex &that) & = delete;

		////////////////////////////////////////////////////////////////

		D operator*() const {
			return (*vector)(index);
		}

		const Index *operator&() const {
			return idxs;
		}

	};

	template <bool isConst>
	class VectorIterator final : public std::iterator<
			std::forward_iterator_tag, VectorIndex<isConst>> {

		typedef typename std::conditional<isConst, const Vector,
				Vector>::type V;

		const Grid *grid;
		V *vector;
		Index index;

	public:

		VectorIterator(const Grid &grid, V &vector, Index index = 0)
				noexcept : grid(&grid), vector(&vector),
				index(index) {
		}

		VectorIterator(const Grid &grid, V &&vector, Index index)
				= delete;
		VectorIterator(Grid &&grid, V &vector, Index index) = delete;
		VectorIterator(Grid &&grid, V &&vector, Index index) = delete;

		VectorIterator(const VectorIterator &that) noexcept
				: grid(that.grid), vector(that.vector),
				index(that.index) {
		}

		VectorIterator &operator=(const VectorIterator &that) &
				= delete;

		////////////////////////////////////////////////////////////////

		bool operator==(const VectorIterator &that) const {
			return grid == that.grid && vector == that.vector
					&& index == that.index;
		}

		bool operator!=(const VectorIterator &that) const {
			return !(*this == that);
		}

		VectorIndex<isConst> operator*() const {
			return VectorIndex<isConst>(*grid, *vector, index);
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

	template <bool isConst>
	class GridVector final {

		typedef typename std::conditional<isConst, const Vector,
				Vector>::type V;
		typedef typename std::conditional<isConst, Real,
				Real &>::type D;

		const Grid *grid;
		V *vector;

	public:

		GridVector(const Grid &grid, V &vector) noexcept : grid(&grid),
				vector(&vector) {
			assert(vector.size() == grid.size());
		}

		GridVector(const Grid &grid, V &&vector) = delete;
		GridVector(Grid &&grid, V &&vector) = delete;
		GridVector(Grid &&grid, V &vector) = delete;

		GridVector(const GridVector &that) noexcept : grid(that.grid),
				vector(that.vector) {
		}

		GridVector &operator=(const GridVector &that) & = delete;

		////////////////////////////////////////////////////////////////

		Real operator()(const Real *coordinates) const;

		template <typename ...Args>
		Real operator()(Real c0, Args ...coordinates) const {
			static_assert(dim - 1 == sizeof...(Args),
					"The number of arguments must be "
					"consistent with the dimension");

			Real coords[] {c0, coordinates...};
			return (*this)(coords);
		}

		D operator()(const Index *indices) const {
			return (*vector)(grid->unroll(indices));
		}

		template <typename ...Args>
		D operator()(Index i0, Args ...indices) const {
			static_assert(dim - 1 == sizeof...(Args),
					"The number of arguments must be "
					"consistent with the dimension");

			Index idxs[] {i0, indices...};
			return (*this)(idxs);
		}

		////////////////////////////////////////////////////////////////

		VectorIterator<isConst> begin() const {
			return VectorIterator<isConst>(*grid, *vector);
		}

		VectorIterator<isConst> end() const {
			return VectorIterator<isConst>(*grid, *vector,
					grid->size());
		}

		////////////////////////////////////////////////////////////////

		// TODO: Fix this

		/**
		 * Prettifies and prints the contents of a vector along with the
		 * corresponding grid indices.
		 * @param os The output stream.
		 * @param v_G A A GridVector object.
		 * @see QuantPDE::Grid::accessor
		 */
		/*template <bool B>
		friend std::ostream &operator<<(std::ostream &os,
				const GridVector<B> &v_G) {
			Real coordinates[dim];
			for(auto v_indices : v_G) {
				v_G.grid->coordinates(&v_indices, coordinates);
				for(Index i = 0; i < dim; i++) {
					os << coordinates[i] << '\t';
				}
				os << *v_indices << std::endl;
			}
			return os;
		}*/

	};

	typedef GridVector<true> GridVectorConst;
	typedef GridVector<false> GridVectorNonConst;

	////////////////////////////////////////////////////////////////////////

	class GridMatrix final {

		const Grid *grid;
		Matrix *matrix;

	public:

		GridMatrix(const Grid &grid, Matrix &matrix) : grid(&grid),
				matrix(&matrix) {
			assert(matrix.rows() == grid.size());
			assert(matrix.cols() == grid.size());
		}

		GridMatrix(const GridMatrix &that) = delete;
		GridMatrix &operator=(const GridMatrix &that) & = delete;

		////////////////////////////////////////////////////////////////

		Real &operator()(const Index *indices) const {
			return matrix->insert( grid->unroll(indices),
					grid->unroll(indices + dim) );
		}

		template <typename ...Args>
		Real &operator()(Args ...indices) const {
			static_assert((dim * 2) == sizeof...(Args),
					"The number of arguments must be "
					"consistent with the dimensions");

			Index idxs[] {indices...};
			return matrix->insert( grid->unroll(idxs),
					grid->unroll(idxs + dim) );
		}

	};

	////////////////////////////////////////////////////////////////////////

	public:

	class MatrixBuilder {

		virtual Real &get(const Index *indices) = 0;

	protected:

		const Grid *G;
		Matrix M;

	public:

		MatrixBuilder(const Grid &grid) noexcept : G(&grid) {
		}

		MatrixBuilder(const MatrixBuilder &that) = delete;

		// Do not need virtual constructor
		// virtual ~MatrixBuilder() {
		//}

		MatrixBuilder &operator=(const MatrixBuilder &that) & = delete;

		////////////////////////////////////////////////////////////////

		const Grid &grid() const {
			return *G;
		}

		Real &operator()(const Index *indices) {
			return get(indices);
		}

		template <typename ...Args>
		Real &operator()(Args ...indices) {
			static_assert((dim * 2) == sizeof...(Args),
					"The number of arguments must be "
					"consistent with the dimensions");

			Index idxs[] {indices...};
			return get(idxs);
		}

		virtual const Matrix &matrix() = 0;

	};

	private:

	class LazyMatrixBuilder final : public MatrixBuilder {

		std::vector<Entry> entries;

		virtual Real &get(const Index *indices) {
			entries.push_back(Entry(this->G->unroll(indices),
					this->G->unroll(indices + dim), 0.));
			return entries.back().value();
		}

	public:

		LazyMatrixBuilder(const Grid &grid) noexcept
				: MatrixBuilder(grid) {
		}

		LazyMatrixBuilder(const Grid &grid, size_t nonzeros) noexcept
				: MatrixBuilder(grid) {
			entries.reserve(nonzeros);
		}

		////////////////////////////////////////////////////////////////

		virtual const Matrix &matrix() {
			this->M = Matrix(this->G->size(),
					this->G->size());
			this->M.setFromTriplets(entries.begin(), entries.end());
			return this->M;
		}

	};

	class FastMatrixBuilder final : public MatrixBuilder {

		virtual Real &get(const Index *indices) {
			return this->M.insert(this->G->unroll(indices),
					this->G->unroll(indices + dim));
		}

	public:

		template <typename IV>
		FastMatrixBuilder(const Grid &grid, IV &&profile) noexcept
				: MatrixBuilder(grid) {
			this->M = Matrix(grid.size(), grid.size());
			this->M.reserve(std::forward(profile));
		}

		////////////////////////////////////////////////////////////////

		virtual const Matrix &matrix() {
			this->M.makeCompressed();
			return this->M;
		}

	};

	////////////////////////////////////////////////////////////////////////

	Index vsize;
	Axis axes[dim]; // Grid has access to private empty constructor Axis()

	////////////////////////////////////////////////////////////////////////
	// (un)roll
	////////////////////////////////////////////////////////////////////////

	void roll(Index index, Index *indices) const {
		Index m = 1;

		// Unroll loop
		for(Index i = 0; i < dim - 1; i++) {
			m *= axes[i].size();
		}

		// Unroll loop
		for(Index i = dim - 1; i > 0; i++) {
			indices[i] = index / m;
			index -= indices[i] * m;
			m /= axes[i--].size();
		}
		indices[0] = index;
	}

	Index unroll(const Index *indices) const {
		Index
			index = 0,
			horner = 1
		;

		index = indices[0];

		// Unroll loop
		for(Index i = 1; i < dim; i++) {
			horner *= axes[i - 1].size();
			index += horner * indices[i];
		}

		return index;
	}

	////////////////////////////////////////////////////////////////////////

public:

	/**
	 * Constructor.
	 */
	template <typename ...Args>
	Grid(Args ...args) noexcept {
		static_assert(dim == sizeof...(Args),
				"The number of arguments must be consistent "
				"with the dimensions");

		vsize = 1;
		Axis axes[] = {args...};

		// Unroll loop
		for(Index i = 0; i < dim; i++) {
			assert(axes[i].size() > 0);

			vsize *= axes[i].size();
			this->axes[i] = axes[i];
		}
	}

	/**
	 * Copy constructor.
	 */
	Grid(const Grid &that) noexcept : vsize(that.vsize) {
		// Unroll loop
		for(Index i = 0; i < dim; i++) {
			axes[i] = that.axis[i];
		}
	}
	/**
	 * Move constructor.
	 */
	Grid(Grid &&that) noexcept : vsize(that.vsize) {
		// Unroll loop
		for(Index i = 0; i < dim; i++) {
			axes[i] = std::move(that.axis[i]);
		}
	}

	/**
	 * Copy assignment operator.
	 */
	Grid &operator=(const Grid &that) & noexcept {
		vsize = that.vsize;

		// Unroll loop
		for(Index i = 0; i < dim; i++) {
			axes[i] = that.axes[i];
		}

		return *this;
	}

	/**
	 * Move assignment operator.
	 */
	Grid &operator=(Grid &&that) & noexcept {
		vsize = that.vsize;

		// Unroll loop
		for(Index i = 0; i < dim; i++) {
			axes[i] = std::move(that.axes[i]);
		}

		return *this;
	}

	/**
	 * Creates and returns an accessor to a vector.
	 * Accessors can be used to easily index into vectors. For example,
	 * \code{.cpp}
	 * // Create a 2-dimensional grid
	 * Axis X {0., .25, .5, .75, 1.};
	 * Axis Y {0., .5, 1.};
	 * Grid2 G(X, Y);
	 *
	 * // Create a vector on the grid and an accessor to go with it
	 * Vector v = G.vector();
	 * auto v_G = G.accessor(v);
	 *
	 * // v_{i,j} = X_i * Y_j
	 * for(Index i = 0; i < X.size(); i++) {
	 * 	for(Index j = 0; j < Y.size(); j++) {
	 * 		v_G(i, j) = X(i) * Y(j);
	 * 	}
	 * }
	 *
	 * // Query a value not necessarily on the grid using linear
	 * // interpolation
	 * std::cout << v_g(.6, .75) << std::endl;
	 * \endcode
	 * Alternatively, one can also use accessors as iterable objects to
	 * accomplish the same goal:
	 * \code{.cpp}
	 * // Foreach loop over the elements in the vector
	 * for(auto v_ij : G.accessor(v)) {
	 * 	// Get the indices i and j associated with this node
	 * 	auto ij = &v_ij;
	 * 	Index i = ij[0], j = ij[1];
	 *
	 * 	// Set the value at this node
	 * 	*v_ij = X(i) * Y(j);
	 * }
	 * \endcode
	 * The latter method is more useful when the dimension (i.e. number of
	 * axes) of the grid is not known ahead of time.
	 * @return An accessor to a constant vector.
	 */
	GridVectorConst accessor(const Vector &vector) const {
		return GridVectorConst(*this, vector);
	}

	/**
	 * @return An accessor to a vector.
	 * @see QuantPDE::Grid::accessor
	 */
	GridVectorNonConst accessor(Vector &vector) const {
		return GridVectorNonConst(*this, vector);
	}

	// TODO: handle rvalue refeferences with std::move
	GridVectorNonConst accessor(Vector &&vector) = delete;

	/**
	 * Creates and returns an accessor to a matrix.
	 * A matrix accessor can only be used to insert nonzero entries into
	 * a matrix. Furthermore, attempting to insert into a location twice
	 * will cause undefined behaviour. The example below builds the discrete
	 * Laplace operator.
	 * \code{.cpp}
	 * // Create a 2-dimensional grid
	 * Grid2 G(X, Y);
	 * Index m = X.size();
	 * Index n = Y.size();
	 *
	 * // Create a matrix with 4 nonzero element reserved per row
	 * Matrix M = G.matrix();
	 * M.reserve(IntegerVector::Constant(G.size(), 4));
	 *
	 * auto M_G = G.accessor(M);
	 *
	 * // Branching inside the loop is expensive and should be eliminated
	 * // in production code
	 * for(Index i = 0; i < X.size(); i++) {
	 * 	for(Index j = 0; j < Y.size(); j++) {
	 * 		if(i > 0)     M_G(i, j, i - 1, j    ) = -1.;
	 * 		if(j > 0)     M_G(i, j, i,     j - 1) = -1.;
	 * 		              M_G(i, j, i,     j    ) =  4.;
	 * 		if(j < n - 1) M_G(i, j, i,     j + 1) = -1.;
	 * 		if(i < m - 1) M_G(i, j, i + 1, j    ) = -1.;
	 * 	}
	 * }
	 *
	 * // The following results in undefined behaviour since we have already
	 * // assigned a value to this position:
	 * // M_G(0, 0, 0, 0) = 4.;
	 *
	 * M.makeCompressed();
	 * \endcode
	 * @return An accessor to a matrix.
	 */
	GridMatrix accessor(Matrix &matrix) const {
		return GridMatrix(*this, matrix);
	}

	// No good reason to pass an rvalue to a matrix accessor
	GridMatrix accessor(Matrix &&matrix) = delete;

	/**
	 * Creates and returns a matrix builder. This should be used when the
	 * approximate number of nonzero entries per row (if the underlying
	 * implementation is row-major) is not known ahead of time. The example
	 * below builds the discrete Laplace operator.
	 * \code{.cpp}
	 * // Create a 1-dimensional grid
	 * Grid2 G(X, Y);
	 * Index m = X.size();
	 * Index n = Y.size();
	 *
	 * auto builder = G.matrixBuilder();
	 *
	 * // Branching inside the loop is expensive and should be eliminated
	 * // in production code
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
	 * A call to matrix() should only occur once. The results are undefined
	 * otherwise. Similarly, nonzero entries should only be specified once.
	 * @return A matrix builder.
	 */
	LazyMatrixBuilder matrixBuilder() const {
		return LazyMatrixBuilder(*this);
	}

	/**
	 * Creates and returns a matrix builder.
	 * @param nonzeros The approximate total number of nonzero entries in
	 *                 the resulting matrix.
	 */
	LazyMatrixBuilder matrixBuilder(size_t nonzeros) const {
		return LazyMatrixBuilder(*this, nonzeros);
	}

	/**
	 * Creates and returns a matrix builder. When the number of nonzeros per
	 * row (if the underlying implementation is row-major) is known a
	 * priori, this method should be used.
	 * @param nonzeros An integer vector whose i-th element represents the
	 *                 number of elements in the i-th row (if the underlying
	 *                 implementation is row-major) of the matrix.
	 */
	template <typename IV>
	FastMatrixBuilder matrixBuilder(IV &&nonzeros) const {
		return FastMatrixBuilder(*this, std::forward(nonzeros));
	}

	/**
	 * @return A square matrix of order equal to the number of nodes on this
	 *         grid.
	 */
	Matrix matrix() const {
		return Matrix(vsize, vsize);
	}

	/**
	 * @return A zero vector of length equal to the number of nodes on this
	 *         grid.
	 */
	Vector zero() const {
		return Vector::Zero(vsize);
	}

	/**
	 * @return A vector of ones of length equal to the number of nodes on
	 *         this grid.
	 */
	Vector ones() const {
		return Vector::Ones(vsize);
	}

	/**
	 * @return A vector of length equal to the number of nodes on this grid.
	 *         The contents of the vector are not guaranteed to be any
	 *         particular values.
	 */
	Vector vector() const {
		return Vector(vsize);
	}

	/**
	 * Converts indices to coordinates on this grid.
	 * @param indices An array of indices of length equal to the dimension.
	 * @param coordinates An array of length equal to the dimension. The
	 *                    i-th element of this array is set to the
	 *                    corresponding tick on the i-th axis.
	 */
	void coordinates(const Index *indices, Real *coordinates) const {
		// Unroll loop
		for(Index i = 0; i < dim; i++) {
			coordinates[i] = axes[i][indices[i]];
		}
	}

	/**
	 * @param index An index.
	 * @return The index-th axis.
	 */
	const Axis &operator[](Index index) const {
		return axes[index];
	}

	/**
	 * @param function A function.
	 * @return The image of a function on this grid as a vector.
	 */
	template <typename F>
	Vector image(F &&function) const {
		Vector v = vector();

		Real coords[dim];
		for(auto node : accessor(v)) {
			coordinates(&node, coords);
			*node = QUANT_PDE_UNPACK_AND_CALL(
					std::forward<F>(function), coords, dim);
		}

		return v;
	}

	/**
	 * @return The total number of nodes on this grid.
	 */
	Index size() const {
		return vsize;
	}

	/**
	 * Prettifies and prints the grid.
	 * @param os The output stream.
	 * @param grid A grid.
	 */
	friend std::ostream &operator<<(std::ostream &os,
			const Grid<dim> &grid) {
		os << grid[0];
		for(Index i = 1; i < dim; i++) {
			os << " x " << grid[i];
		}
		return os;
	}

};

typedef Grid<1> Grid1;
typedef Grid<2> Grid2;
typedef Grid<3> Grid3;

}

#endif

