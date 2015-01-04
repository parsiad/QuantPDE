#ifndef QUANT_PDE_CORE_IMPULSE
#define QUANT_PDE_CORE_IMPULSE

#include <array>   // std::array
#include <cstdint> // std::intmax_t
#include <utility> // std::forward

namespace QuantPDE {

/**
 * Represents an impulse of the form \f$V - \mathcal{M}V\f$, where
 * \f$
 * \left[ \mathcal{M} V \right] \left( t, \mathbf{x}, \mathbf{q} \right) \equiv
 * V \left(
 * 	f_{1} \left( t, \mathbf{x}, \mathbf{q} \right),
 * 	\ldots,
 * 	f_{n} \left( t, \mathbf{x}, \mathbf{q} \right)
 * \right)
 * + g \left( t, \mathbf{x}, \mathbf{q} \right)
 * \f$
 * and \f$\mathbf{q}\f$ is a fixed control.
 *
 * The \f$f_i\f$s are referred to as the transition functions, and \f$g\f$ is
 * referred to as the flow.
 *
 * @tparam Dimension The dimension of the Euclidean space associated with
 *                   \f$\mathbf{x}\f$.
 * @tparam ControlDimension The dimension of the Euclidean space associated with
 *                          \f$\mathbf{q}\f$.
 * @tparam Negative If true, the sign is flipped so that the impulse is of the
 *                  form \f$\mathcal{M}V - V\f$.
 */
template <Index Dimension, Index ControlDimension, bool Negative = false>
class Impulse final : public RawControlledLinearSystem<Dimension,
		ControlDimension> {

	typedef IntegerPower<2, Dimension> TwoToTheDimension;

	const RectilinearGrid<Dimension> &grid;

	// t, x1, ..., xn, q1, ..., qn
	typedef Function<1 + Dimension + ControlDimension> F;

	F flow;
	std::array<F, Dimension> transitions;

public:

	/**
	 * Creates an impulse.
	 * @param grid The rectilinear grid associated with \f$\mathbf{x}\f$.
	 * @param flow The flow function.
	 * @param transitions The transition functions.
	 * @see QuantPDE::RectilinearGrid
	 */
	template <typename G, typename F1, typename ...F2>
	Impulse(
		G &grid,
		F1 &&flow,
		F2 &&...transitions
	) noexcept :
		grid(grid),
		flow( std::forward<F1>(flow) ),
		transitions( {{std::forward<F2>(transitions)...}} )
	{
	}

	// TODO: Allow for multiple types of interpolation

	virtual Matrix A(Real t) {
		Matrix M = grid.matrix();
		M.reserve(IntegerVector::Constant(
			grid.size(),
			TwoToTheDimension::value
		));

		Real args[Dimension + ControlDimension + 1];
		Real plus[Dimension];

		// TODO: Fix. Currying explodes in Clang; not sure why...
		/*
		Function<Dimension + ControlDimension> transitions_t[Dimension];
		for(Index i = 0; i < Dimension; ++i) {
			transitions_t[i] =
				curry<1 + Dimension + ControlDimension>(
					transitions[i],
					t
				)
			;
		}
		*/

		args[0] = t;

		Index k = 0;
		for(auto node : grid) {
			// Coordinates
			for(int i = 0; i < Dimension; ++i) {
				args[i + 1] = node[i];
			}

			// Control coordinates
			for(int i = 0; i < ControlDimension; ++i) {
				args[Dimension + i + 1]=(this->control(i))(k);
			}

			// Get new state
			for(int i = 0; i < Dimension; ++i) {
				plus[i] =
					packAndCall<
						  Dimension
						+ ControlDimension
						+ 1
					>(
						transitions[i],
						args
					)
				;
			}

			auto data = linearInterpolationData(grid, plus);

			Index idxs[Dimension];
			for(
				std::intmax_t i = 0;
				i < TwoToTheDimension::value;
				++i
			) {
				Real factor = 1.;

				for(Index j = 0; j < Dimension; ++j) {
					if(i & (1 << j)) {
						// j-th bit of i is 1
						idxs[j] = std::get<0>(data[j]);
						factor *= std::get<1>(data[j]);
					} else {
						// j-th bit of i is 0
						idxs[j] =  std::get<0>(data[j])
								+ 1 ;
						factor *= -std::get<1>(data[j])
								+ 1.;
					}
				}

				// n-Dimensional
				NaryMethodConst<
					Index,
					RectilinearGrid<Dimension>,
					Dimension, Index
				> tmp = &RectilinearGrid<Dimension>::index;
				Index j = packAndCall<Dimension>(grid, tmp,
						idxs);

				M.insert(k, j) = factor;
			};

			++k;
		}

		M.makeCompressed();

		if(Negative) {
			return M - grid.identity();
		}

		return grid.identity() - M;
	}

	virtual Vector b(Real t) {
		Vector b(grid.vector());

		// flow function at time t
		auto flow_t = curry<1 + Dimension + ControlDimension>(flow, t);

		Real args[Dimension + ControlDimension];

		Index k = 0;
		for(auto node : accessor(grid, b)) {

			// Coordinates
			auto coords = &node;
			for(int i = 0; i < Dimension; ++i) {
				args[i] = coords[i];
			}

			// Control coordinates
			for(int i = 0; i < ControlDimension; ++i) {
				args[Dimension + i] = (this->control(i))(k);
			}

			// Set value at node
			*node = (Negative ? -1. : 1.) * packAndCall<Dimension
					+ ControlDimension>(flow_t, args);

			++k;
		}

		return b;
	}

};

template <Index ControlDimension>
using Impulse1 = Impulse<1, ControlDimension>;

template <Index ControlDimension>
using Impulse2 = Impulse<2, ControlDimension>;

template <Index ControlDimension>
using Impulse3 = Impulse<3, ControlDimension>;

typedef Impulse1<1> Impulse1_1;
typedef Impulse1<2> Impulse1_2;
typedef Impulse1<3> Impulse1_3;

typedef Impulse2<1> Impulse2_1;
typedef Impulse2<2> Impulse2_2;
typedef Impulse2<3> Impulse2_3;

typedef Impulse3<1> Impulse3_1;
typedef Impulse3<2> Impulse3_2;
typedef Impulse3<3> Impulse3_3;

template <Index ControlDimension>
using NegativeImpulse1 = Impulse<1, ControlDimension, true>;

template <Index ControlDimension>
using NegativeImpulse2 = Impulse<2, ControlDimension, true>;

template <Index ControlDimension>
using NegativeImpulse3 = Impulse<3, ControlDimension, true>;

typedef NegativeImpulse1<1> NegativeImpulse1_1;
typedef NegativeImpulse1<2> NegativeImpulse1_2;
typedef NegativeImpulse1<3> NegativeImpulse1_3;

typedef NegativeImpulse2<1> NegativeImpulse2_1;
typedef NegativeImpulse2<2> NegativeImpulse2_2;
typedef NegativeImpulse2<3> NegativeImpulse2_3;

typedef NegativeImpulse3<1> NegativeImpulse3_1;
typedef NegativeImpulse3<2> NegativeImpulse3_2;
typedef NegativeImpulse3<3> NegativeImpulse3_3;

}

#endif
