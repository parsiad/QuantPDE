#ifndef QUANT_PDE_GRID_HPP
#define QUANT_PDE_GRID_HPP

#include <algorithm>        // std::swap
#include <cstring>          // std::memcpy
#include <cassert>          // assert
#include <initializer_list> // std::initializer_list
#include <iostream>         // std::ostream

namespace QuantPDE {

class Grid;

/**
 * A rectilinear grid.
 */
class Grid {

	template <bool isConst>
	class VectorIndex {

		typedef typename std::conditional<isConst, const Vector,
				Vector>::type V;
		typedef typename std::conditional<isConst, Real,
				Real &>::type D;

		const Grid *grid;
		V *vector;
		Index index;
		Index *idxs;

	public:

		VectorIndex(const Grid &grid, V &vector, Index index)
				: grid(&grid), vector(&vector), index(index) {
			idxs = new Index[grid.size()];
			grid.roll(index, idxs);
		}

		VectorIndex(const VectorIndex &that) : grid(that.grid),
				vector(that.vector), index(that.index) {
			idxs = new Index[grid->size()];
			std::memcpy(idxs, that.idxs,
					sizeof(Index) * grid->size());
		}

		VectorIndex(VectorIndex &&that) : grid(that.grid),
				vector(that.vector), index(that.index) {
			idxs = that.idxs;
			that.idxs = NULL;
		}

		~VectorIndex() {
			delete [] idxs;
		}

		VectorIndex &operator=(VectorIndex that) {
			grid = that.grid;
			vector = that.vector;
			index = that.index;
			std::swap(idxs, that.idxs);
			return *this;
		}

		D operator*() const {
			return (*vector)(index);
		}

		const Index *operator&() const {
			return idxs;
		}

	};

	template <bool isConst>
	class VectorIterator : public std::iterator<std::forward_iterator_tag,
			VectorIndex<isConst>> {

		typedef typename std::conditional<isConst, const Vector,
				Vector>::type V;

		const Grid *grid;
		V *vector;
		Index index;

	public:

		VectorIterator(const Grid &grid, V &vector, Index index = 0)
				: grid(&grid), vector(&vector), index(index) {
		}

		VectorIterator(const VectorIterator &that) : grid(that.grid),
				vector(that.vector), index(that.index) {
		}

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
	class GridVector {

		typedef typename std::conditional<isConst, const Vector,
				Vector>::type V;
		typedef typename std::conditional<isConst, Real,
				Real &>::type D;

		const Grid &grid;
		V &vector;

		struct U {};

	public:

		GridVector(const Grid &grid, V &vector) : grid(grid),
				vector(vector) {
			assert(vector.size() == grid.size());
		}

		GridVector(const GridVector &that) : grid(that.grid),
				vector(that.vector) {
		}

		///////////////////////////////////////////////////////////////

		Real operator()(const Real *coordinates) const;

		template <typename ...ArgsT>
		Real operator()(Real c0, ArgsT ...coordinates) const {
			Real coords[] {c0, coordinates...};
			return (*this)(coords);
		}

		D operator()(const Index *indices) const {
			return vector(grid.unroll(indices));
		}

		template <typename ...ArgsT>
		D operator()(Index i0, ArgsT ...indices) const {
			Index idxs[] {i0, indices...};
			assert(sizeof(idxs) == grid.size() * sizeof(Index));
			return (*this)(idxs);
		}

		///////////////////////////////////////////////////////////////

		VectorIterator<isConst> begin() const {
			return VectorIterator<isConst>(grid, vector);
		}

		VectorIterator<isConst> end() const {
			return VectorIterator<isConst>(grid, vector,
					grid.nodes());
		}

		///////////////////////////////////////////////////////////////

		template <bool B>
		friend std::ostream &operator<<(std::ostream &,
				const GridVector<B> &);

	};

	typedef GridVector<true> GridVectorConst;
	typedef GridVector<false> GridVectorNonConst;

	class GridMatrix {

		const Grid &grid;
		Matrix &matrix;

	public:

		GridMatrix(const Grid &grid, Matrix &matrix) : grid(grid),
				matrix(matrix) {
			assert(matrix.rows() == grid.nodes());
			assert(matrix.cols() == grid.nodes());
		}

		GridMatrix(const GridMatrix &that) : grid(that.grid),
				matrix(that.matrix) {
		}

		///////////////////////////////////////////////////////////////

		Real &operator()(const Index *indices) const {
			return matrix.insert( grid.unroll(indices),
					grid.unroll(indices + grid.size()) );
		}

		template <typename ...ArgsT>
		Real &operator()(ArgsT ...indices) const {
			Index idxs[] {indices...};
			assert(sizeof(idxs) == grid.size() * sizeof(Index) * 2);
			return matrix.insert( grid.unroll(idxs),
					grid.unroll(idxs + grid.size()) );
		}

	};

	virtual void roll(Index, Index *) const = 0;

	virtual Index unroll(const Index *) const = 0;

public:

	/**
	 * Creates and returns an accessor to a vector.
	 * Accessors can be used to easily index into vectors. For example,
	 * \code{.cpp}
	 * // Create a 2-dimensional grid
	 * Axis X {0., .25, .5, .75, 1.};
	 * Axis Y {0., .5, 1.};
	 * Grid2d G(X, Y);
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
	 * @return An accessor to a vector.
	 */
	GridVectorConst accessor(const Vector &vector) const {
		return GridVectorConst(*this, vector);
	}

	/**
	 * @return An accessor to a constant vector.
	 * @see QuantPDE::Grid::accessor
	 */
	GridVectorNonConst accessor(Vector &vector) const {
		return GridVectorNonConst(*this, vector);
	}

	/**
	 * Creates and returns an accessor to a matrix.
	 * A matrix accessor can only be used to insert nonzero entries into
	 * a matrix. Furthermore, attempting to insert into a location twice
	 * will cause undefined behaviour.
	 * \code{.cpp}
	 * // Create a 2-dimensional grid
	 * Grid2d G(X, Y);
	 *
	 * // Create a matrix with 1 nonzero element reserved per row
	 * Matrix M = G.matrix();
	 * M.reserve(Eigen::VectorXi::Constant(n, 1));
	 *
	 * // A very roundabout way to make the identity
	 * auto M_G = G.accessor(M);
	 * for(Index i = 0; i < X.size(); i++) {
	 * 	for(Index j = 0; j < Y.size(); j++) {
	 * 		M_G(i, j, i, j) = 1.;
	 * 		// The first pair of indices is used to resolve the row
	 * 		// index of the matrix, while the second pair is used to
	 * 		// resolve the column index of the matrix
	 * 	}
	 * }
	 *
	 * // The following results in undefined behaviour since we have already
	 * // assigned a value to this position:
	 * // M_G(0, 0, 0, 0) = 1.;
	 *
	 * M.makeCompressed();
	 * \endcode
	 * @return An accessor to a matrix.
	 */
	GridMatrix accessor(Matrix &matrix) const {
		return GridMatrix(*this, matrix);
	}

	/**
	 * @return A square matrix of order equal to the number of nodes on this
	 *         grid.
	 */
	Matrix matrix() const {
		return Matrix(nodes(), nodes());
	}

	/**
	 * @return A zero vector of length equal to the number of nodes on this
	 *         grid.
	 */
	Vector zero() const {
		return Vector::Zero(nodes());
	}

	/**
	 * @return A vector of ones of length equal to the number of nodes on
	 *         this grid.
	 */
	Vector ones() const {
		return Vector::Ones(nodes());
	}

	/**
	 * @return A vector of length equal to the number of nodes on this grid.
	 *         The contents of the vector are not guaranteed to be any
	 *         particular values.
	 */
	Vector vector() const {
		return Vector(nodes());
	}

	/**
	 * Converts indices to coordinates on this grid.
	 * @param indices An array of indices of length size().
	 * @param coordinates An array of size size(). The i-th element of this
	 *                    array is set to the corresponding tick on the i-th
	 *                    axis.
	 */
	virtual void coordinates(const Index *indices, Real *coordinates)
			const = 0;

	/**
	 * @param index An index.
	 * @return The index-th axis.
	 */
	virtual const Axis &operator()(Index index) const = 0;

	/**
	 * @return The number of axes (i.e. the dimension of the domain on which
	 *         this rectilinear grid lives).
	 */
	virtual Index size() const = 0;

	/**
	 * @return The total number of nodes on this grid.
	 */
	virtual Index nodes() const = 0;

};

/**
 * Prettifies and prints a grid to an output stream.
 * @param os The output stream.
 * @param axis The axis.
 */
std::ostream &operator<<(std::ostream &os, const Grid &grid) {
	os << grid(0);
	for(Index i = 1; i < grid.size(); i++) {
		os << " \u00D7 " << grid(i);
	}
	return os;
}


/**
 * A rectilinear grid of arbitrary dimension. Use of this should be avoided
 * when the dimension of the grid is small (and known a priori).
 */
class GridXd final : public Grid {

	Axis **axes;
	Index dim, vsize;

	void initialize(std::initializer_list<Axis> list) {
		dim = list.size();
		assert(dim > 0);

		vsize = 1;
		Axis **p = axes = new Axis*[list.size()];
		for(Axis axis : list) {
			assert(axis.size() > 0);

			vsize *= axis.size();
			*(p++) = new Axis(axis);
		}
	}

	///////////////////////////////////////////////////////////////////////
	// (un)roll
	///////////////////////////////////////////////////////////////////////

	virtual void roll(Index index, Index *indices) const {
		Index m = 1;
		for(Index i = 0; i < size() - 1; i++) {
			m *= axes[i]->size();
		}

		for(Index i = size() - 1; i > 0; i++) {
			indices[i] = index / m;
			index -= indices[i] * m;
			m /= axes[i--]->size();
		}
		indices[0] = index;
	}

	virtual Index unroll(const Index *indices) const {
		Index
			index = 0,
			horner = 1
		;

		const Index *end = indices + dim;
		Axis **p = axes;
		index = *(indices++);
		while(indices != end) {
			horner *= (*(p++))->size();
			index += horner * (*(indices++));
		}

		return index;
	}

	///////////////////////////////////////////////////////////////////////

public:

	/**
	 * Constructor.
	 * @param axes The axes that make up the rectilinear grid.
	 */
	template <typename ...ArgsT>
	GridXd(ArgsT ...axes) {
		initialize({axes...});
	}

	/**
	 * Constructor.
	 * @param list A list of axes that make up the rectilinear grid.
	 */
	GridXd(std::initializer_list<Axis> list) {
		initialize(list);
	}

	/**
	 * Copy constructor.
	 */
	GridXd(const GridXd &that) {
		dim = that.dim;
		vsize = that.vsize;

		Axis
			**p = axes,
			**q = that.axes,
			**end = axes + dim
		;

		while(p != end) {
			*(p++) = *(q++);
		}
	}

	/**
	 * Move constructor.
	 */
	GridXd(GridXd &&that) {
		dim = that.dim;
		vsize = that.vsize;
		axes = that.axes;

		that.axes = NULL;
	}

	/**
	 * Destructor.
	 */
	~GridXd() {
		Axis
			**p = axes,
			**end = axes + dim
		;

		while(p != end) {
			delete *(p++);
		}

		delete [] axes;
	}

	/**
	 * Assignment operator.
	 */
	GridXd &operator=(GridXd that) {
		dim = that.dim;
		vsize = that.vsize;
		std::swap(axes, that.axes);
		return *this;
	}

	virtual void coordinates(const Index *indices, Real *coordinates)
			const {
		Axis **p = axes, **end = axes + dim;
		while(p != end) {
			*(coordinates++) = (**(p++))(*(indices++));
		}
	}

	virtual const Axis &operator()(Index index) const {
		return *axes[index];
	}

	virtual Index size() const {
		return dim;
	}

	virtual Index nodes() const {
		return vsize;
	}

};

/**
 * Prettifies and prints the contents of a vector along with the corresponding
 * grid indices.
 * @param os The output stream.
 * @param v_G A reference to a GridVector object.
 * @see QuantPDE::Grid::accessor
 */
template <bool isConst>
std::ostream &operator<<(std::ostream &os,
		const Grid::GridVector<isConst> &v_G) {
	Real *coordinates = new Real[v_G.grid.size()];
	for(auto v_indices : v_G) {
		v_G.grid.coordinates(&v_indices, coordinates);
		for(Index i = 0; i < v_G.grid.size(); i++) {
			os << coordinates[i] << '\t';
		}
		os << *v_indices << std::endl;
	}
	delete [] coordinates;
	return os;
}

// TODO: Grid1d, Grid2d, Grid3d
typedef GridXd Grid1d;
typedef GridXd Grid2d;
typedef GridXd Grid3d;

}

#endif

