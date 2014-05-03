#ifndef QUANT_PDE_FUNCTION_HPP
#define QUANT_PDE_FUNCTION_HPP

namespace QuantPDE {

/**
 * Represents a function \f$f\f$ defined on some domain.
 */
class Function {

	virtual double get(const double *) const = 0;

public:

	/**
	 * @param coordinates Coordinates \f$\left\{c_i\right\}\f$ (passed in
	 *                    as an array) specifying a point in \f$f\f$'s
	 *                    domain.
	 * @return The value of \f$f\f$ at the specified coordinates (i.e.
	 *         \f$f\left(c_0, \ldots, c_{n-1}\right)\f$).
	 */
	double operator()(const double *coordinates) const {
		return get(coordinates);
	}

	/**
	 * @param coordinates Coordinates \f$\left\{c_i\right\}\f$ specifying
	 *                    a point in \f$f\f$'s domain.
	 * @return The value of \f$f\f$ at the specified coordinates (i.e.
	 *         \f$f\left(c_0, \ldots, c_{n-1}\right)\f$).
	 */
	template <typename ...ArgsT>
	double operator()(ArgsT ...coordinates) const {
		double coords[] {coordinates...};
		return get(coords);
	}

	/**
	 * @param grid A grid.
	 * @return The image of this function on the set of grid nodes as a
	 *         vector.
	 */
	virtual Vector image(const Grid &grid) const {
		Vector v = grid.vector();

		double *coordinates = new double[grid.size()];
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

	virtual double get(const double *coordinates) const {
		return get(*coordinates);
	}

	virtual double get(double) const = 0;

public:

	virtual Vector image(const Grid &grid) const {
		Vector v = grid.vector();

		double coordinates[1];
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

	double c;

	virtual double get(const double *coordinates) const {
		return c;
	}

public:

	/**
	 * Constructor.
	 * @param c The value of the function everywhere.
	 */
	Constant(double c) : c(c) {
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

	double interpolate(const Index *indices, const double *start,
			const double *end) const {
		if(start == end) {
			return grid.accessor(vector)(indices);
		}

		Index dim = grid.size();
		Index *other = new Index[dim];
		memcpy(other, indices, sizeof(Index) * dim);
		other[dim - (end - start)]++;

		double v = (*start)
				* interpolate(indices, start + 1, end)
				+ (1 - *start) * interpolate(other, start + 1,
				end);

		delete [] other;

		return v;
	}

	virtual double get(const double *coordinates) const {
		Index dim = grid.size();

		Index *indices = new Index[dim];
		double *weights = new double[dim];

		{
			for(Index i = 0; i < grid.size(); i++, coordinates++) {
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

				// Binary search to find node
				Index lo = 0, hi = length - 2, mid = 0;
				double weight = 0.;
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
		}

		double v = interpolate(indices, weights, weights + dim);

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
double Grid::GridVector<isConst>::operator()(const double *coordinates) const {
	return (PiecewiseLinear(grid, vector))(coordinates);
}

}

#endif

