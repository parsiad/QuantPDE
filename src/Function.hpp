#ifndef QUANT_PDE_FUNCTION_HPP
#define QUANT_PDE_FUNCTION_HPP

#include <cstring> // std::memcpy

namespace QuantPDE {

/**
 * Represents a function \f$f\f$ defined on some domain.
 */
class Function {

	virtual Real get(const Real *) const = 0;

public:

	/**
	 * Destructor.
	 */
	virtual ~Function() {
	}

	/**
	 * @param coordinates Coordinates \f$\left\{c_i\right\}\f$ (passed in
	 *                    as an array) specifying a point in \f$f\f$'s
	 *                    domain.
	 * @return The value of \f$f\f$ at the specified coordinates (i.e.
	 *         \f$f\left(c_0, \ldots, c_{n-1}\right)\f$).
	 */
	Real operator()(const Real *coordinates) const {
		return get(coordinates);
	}

	/**
	 * @param coordinates Coordinates \f$\left\{c_i\right\}\f$ specifying
	 *                    a point in \f$f\f$'s domain.
	 * @return The value of \f$f\f$ at the specified coordinates (i.e.
	 *         \f$f\left(c_0, \ldots, c_{n-1}\right)\f$).
	 */
	template <typename ...ArgsT>
	Real operator()(ArgsT ...coordinates) const {
		Real coords[] {coordinates...};
		return get(coords);
	}

	/**
	 * @param grid A grid.
	 * @return The image of this function on the set of grid nodes as a
	 *         vector.
	 */
	virtual Vector image(const Grid &grid) const {
		Vector v = grid.vector();

		Real *coordinates = new Real[grid.size()];
		for(auto node : grid.accessor(v)) {
			grid.coordinates(&node, coordinates);
			*node = get(coordinates);
		}
		delete [] coordinates;

		return v;
	}

};

/**
 * Represents a function defined on a one-dimensional domain.
 */
class Function1d : public Function {

	virtual Real get(const Real *coordinates) const {
		return get(*coordinates);
	}

	virtual Real get(Real) const = 0;

public:

	virtual Vector image(const Grid &grid) const {
		Vector v = grid.vector();

		Real coordinates[1];
		for(auto node : grid.accessor(v)) {
			grid.coordinates(&node, coordinates);
			*node = get(*coordinates);
		}

		return v;
	}

};

// TODO: Function2d, Function3d, Function4d

/**
 * A function \f$f\f$ whose value is constant (i.e. \f$f\left(x\right)=c\f$ for
 * all \f$x\f$ in the domain of \f$f\f$).
 */
class Constant : public Function {

	Real c;

	virtual Real get(const Real *coordinates) const {
		return c;
	}

public:

	/**
	 * Constructor.
	 * @param c The value of the function everywhere.
	 */
	Constant(Real c) : c(c) {
	}

	/**
	 * Copy constructor.
	 */
	Constant(const Constant &that) : c(that.c) {
	}

};

/**
 * A function \f$f\f$ defined piecewise whose pieces are affine functions.
 * Suppose \f$f\f$ is defined on a grid composed of \f$n\f$ axes, and that the
 * \f$i\f$-th axis is a partition consisting of \f$m_i\f$ nodes.
 * Querying \f$f\f$ at a point involves a binary search on each axis.
 * This procedure has worst-case complexity
 * \f$\sim \sum_i \lg m_i \leq n \lg m\f$, where
 * \f$m \equiv \max\left\{m_i\right\}\f$.
 */
class PiecewiseLinear : public Function {

	const Grid &grid;
	const Vector &vector;

	Real interpolate(Index *indices, Real *weights, Index n = 0) const {
		const Index dim = grid.size();

		// Base case
		if(n == dim) {
			return grid.accessor(vector)(indices);
		}

		Index *stride = indices + ( (dim - n) * dim );
		std::memcpy(stride, indices, sizeof(Index) * dim);
		stride[dim - n]++;

		return weights[n] * interpolate(indices, weights, n + 1)
				+ (1 - weights[n]) * interpolate(stride,
				weights, n + 1);
	}

	virtual Real get(const Real *coordinates) const {
		const Index dim = grid.size();

		Real *weights = new Real[dim];

		// We need 2^dim arrays of size dim to store indices for the
		// interpolation routine. Initialize all of this at once.
		Index *indices;
		{
			Index power = 1;
			for(Index i = 0; i < dim; i++) {
				power *= 2;
			}

			indices = new Index[power * dim];
		}

		// For the i-th coordinate, find the ticks on the i-th axis that
		// it lies between along with the distance from the leftmost
		// tick
		for(Index i = 0; i < dim; i++, coordinates++) {
			const Axis &x = grid(i);
			Index length = x.size();

			if(*coordinates <= x(0)) {
				indices[i] = 0;
				weights[i] = 1.;
				continue;
			}

			if(*coordinates >= x(length - 1)) {
				indices[i] = length - 2;
				weights[i] = 0.;
				continue;
			}

			// Binary search to find tick
			Index lo = 0, hi = length - 2, mid = 0;
			Real weight = 0.;
			while(lo <= hi) {
				mid = (lo + hi) / 2;
				if(*coordinates < x(mid)) {
					hi = mid - 1;
				} else if(*coordinates >= x(mid + 1)) {
					lo = mid + 1;
				} else {
					weight = (x(mid + 1)
							- *coordinates)
							/ ( x(mid + 1)
							- x(mid) );
					break;
				}
			}

			indices[i] = mid;
			weights[i] = weight;
		}

		const Real v = interpolate(indices, weights);

		delete [] indices;
		delete [] weights;

		return v;
	}

public:

	/**
	 * Constructor.
	 * @param grid A grid.
	 * @param vector A vector containing the values of the function at the
	 *               grid nodes.
	 */
	PiecewiseLinear(const Grid &grid, const Vector &vector) : grid(grid),
			vector(vector) {
	}

	/**
	 * Copy constructor.
	 */
	PiecewiseLinear(const PiecewiseLinear &that) : grid(that.grid),
			vector(that.vector) {
	}

};

template <bool isConst>
Real Grid::GridVector<isConst>::operator()(const Real *coordinates) const {
	return (PiecewiseLinear(grid, vector))(coordinates);
}

}

#endif

