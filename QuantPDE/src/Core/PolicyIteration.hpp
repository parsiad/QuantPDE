#ifndef QUANT_PDE_CORE_POLICY_ITERATION_HPP
#define QUANT_PDE_CORE_POLICY_ITERATION_HPP

#include <cassert>    // assert
#include <functional> // std::greater, std::less
#include <limits>     // std::numeric_limits

namespace QuantPDE {

// gcc complains about SignedInfinity for some reason in debug mode; disable it

#if 0
//#ifdef NDEBUG
/**
 * Positive or negative infinity.
 * @tparam Positive If true, value becomes positive infinity. Negative infinity
 *                  otherwise.
 */
template <bool Positive>
struct SignedInfinity {

	// Negative/positive infinity assert
	static_assert(std::numeric_limits<Real>::is_iec559,
			"IEEE 754 required");

	/**
	 * Positive or negative infinity.
	 */
	static constexpr Real value =
			+std::numeric_limits<Real>::infinity();

};

template <>
struct SignedInfinity<false> {

	static constexpr Real value =
			-std::numeric_limits<Real>::infinity();

};
//#endif
#endif

template <Index Dimension, Index ControlDimension, bool Max>
class PolicyIteration : public IterationNode {

	typedef typename std::conditional<Max, std::greater<Real>,
			std::less<Real>>::type Order;

	const Domain<Dimension> *domain;
	const Domain<ControlDimension> *controlDomain;
	ControlledLinearSystemBase *system;

	virtual void onIterationStart() {
		NaryMethodNonConst<void, ControlledLinearSystemBase,
				ControlDimension, Vector &&> setInputs =
				&ControlledLinearSystemBase::setInputs;

		// Optimal controls
		Vector optimal[ControlDimension];

		// We cannot template a virtual function, so we are stuck
		// initializing this array the old-fashioned way
		for(Index i = 0; i < ControlDimension; ++i) {
			optimal[i] = domain->vector();
		}

		// Best configuration; initialize to plus-minus infinity
		//#ifdef NDEBUG
		//Vector best = domain->ones() * SignedInfinity<!Max>::value;
		//#else
		Vector best = domain->ones()
				* std::numeric_limits<Real>::infinity();
		if(Max) {
			best *= -1;
		}
		//#endif

		for(auto node : *controlDomain) {
			Vector inputs[ControlDimension];
			for(Index i = 0; i < ControlDimension; ++i) {
				inputs[i] = domain->ones() * node[i];
			}
			packMoveAndCall<ControlDimension>(*system, setInputs,
					inputs);

			Vector bCandidate = system->b( nextTime() );
			Matrix ACandidate = system->A( nextTime() );

			// Compute A(q)x - b(q)
			Vector candidate = ACandidate * iterand(0) - bCandidate;

			for(Index i = 0; i < domain->size(); ++i) {
				if( Order()(candidate(i), best(i)) ) {
					best(i) = candidate(i);
					for(Index j = 0; j < ControlDimension;
							++j) {
						optimal[j](i) = node[j];
					}
				}
			}
		}

		packMoveAndCall<ControlDimension>(*system, setInputs, optimal);
	}

public:

	template <typename D, typename C>
	PolicyIteration(D &domain, C &controlDomain,
			ControlledLinearSystemBase &system) noexcept
			: domain(&domain), controlDomain(&controlDomain),
			system(&system) {
	}

	virtual Matrix A(Real t) {
		assert( t == nextTime() );
		return system->A(t);
	}

	virtual Vector b(Real t) {
		assert( t == nextTime() );
		return system->b(t);
	}

	// TODO: Explicit discretizations

};

////////////////////////////////////////////////////////////////////////////////

template <Index Dimension, Index ControlDimension>
using MinPolicyIteration = PolicyIteration<Dimension, ControlDimension, false>;

template <Index Dimension, Index ControlDimension>
using MaxPolicyIteration = PolicyIteration<Dimension, ControlDimension, true>;

template <Index ControlDimension>
using MinPolicyIteration1 = MinPolicyIteration<1, ControlDimension>;
template <Index ControlDimension>
using MinPolicyIteration2 = MinPolicyIteration<2, ControlDimension>;
template <Index ControlDimension>
using MinPolicyIteration3 = MinPolicyIteration<3, ControlDimension>;

template <Index ControlDimension>
using MaxPolicyIteration1 = MaxPolicyIteration<1, ControlDimension>;
template <Index ControlDimension>
using MaxPolicyIteration2 = MaxPolicyIteration<2, ControlDimension>;
template <Index ControlDimension>
using MaxPolicyIteration3 = MaxPolicyIteration<3, ControlDimension>;

typedef MinPolicyIteration1<1> MinPolicyIteration1_1;
typedef MinPolicyIteration1<2> MinPolicyIteration1_2;
typedef MinPolicyIteration1<3> MinPolicyIteration1_3;

typedef MinPolicyIteration2<1> MinPolicyIteration2_1;
typedef MinPolicyIteration2<2> MinPolicyIteration2_2;
typedef MinPolicyIteration2<3> MinPolicyIteration2_3;

typedef MinPolicyIteration3<1> MinPolicyIteration3_1;
typedef MinPolicyIteration3<2> MinPolicyIteration3_2;
typedef MinPolicyIteration3<3> MinPolicyIteration3_3;

typedef MaxPolicyIteration1<1> MaxPolicyIteration1_1;
typedef MaxPolicyIteration1<2> MaxPolicyIteration1_2;
typedef MaxPolicyIteration1<3> MaxPolicyIteration1_3;

typedef MaxPolicyIteration2<1> MaxPolicyIteration2_1;
typedef MaxPolicyIteration2<2> MaxPolicyIteration2_2;
typedef MaxPolicyIteration2<3> MaxPolicyIteration2_3;

typedef MaxPolicyIteration3<1> MaxPolicyIteration3_1;
typedef MaxPolicyIteration3<2> MaxPolicyIteration3_2;
typedef MaxPolicyIteration3<3> MaxPolicyIteration3_3;

}

#endif
