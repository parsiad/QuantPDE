#ifndef QUANT_PDE_CORE_PENALTY_METHOD_HPP
#define QUANT_PDE_CORE_PENALTY_METHOD_HPP

#include <vector> // std::vector

namespace QuantPDE {

/**
 * Used to solve a problem of the form
 * \f$\min\left( LV, L^\prime V \right)=0\f$.
 */
class PenaltyMethod : public IterationNode {

	const DomainBase *domain;
	LinearSystem *left, *right;
	Real scale;
	Real direct;

	Matrix P, Q;

	Matrix rA, lA;
	Vector rb, lb;

	#ifdef QUANT_PDE_MODULES_HJBQVI_ITERATED_OPTIMAL_STOPPING
	const bool right_explicit;
	#endif

	virtual void onIterationStart() {
		// Only need to do this once
		rA = right->A( nextTime() );
		rb = right->b( nextTime() );
		lA =  left->A( nextTime() );
		lb =  left->b( nextTime() );

		// Evaluate the predicate using the previous iterand
		Vector predicate = rA * iterand(0) - rb;

		Vector compare = domain->zero();
		if(direct) {
			compare = lA * iterand(0) - lb;
		}

		// Initialize penalty matrix
		P.setZero();
		P.reserve( IntegerVector::Constant( domain->size(), 1 ) );

		// Initialize other matrix
		Q.setZero();
		Q.reserve( IntegerVector::Constant( domain->size(), 1 ) );

		// Build penalty matrix
		for(Index i = 0; i < domain->size(); ++i) {
			const bool c = (scale * predicate(i)) < compare(i);
			if(c) { P.insert(i, i) = scale; }
			if(!direct || !c) { Q.insert(i, i) = 1.; }
		}
	}

public:

	template <typename D>
	PenaltyMethod(
		D &domain,
		LinearSystem &constraint, LinearSystem &penalizedConstraint,
		Real tolerance = QuantPDE::tolerance,
		bool direct = false
		#ifdef QUANT_PDE_MODULES_HJBQVI_ITERATED_OPTIMAL_STOPPING
		, bool right_explicit = false
		#endif
	) noexcept :
		domain(&domain),
		left(&constraint), right(&penalizedConstraint),
		scale( 1. / tolerance ),
		direct(direct),
		P( domain.size(), domain.size() ),
		Q( domain.size(), domain.size() )
		#ifdef QUANT_PDE_MODULES_HJBQVI_ITERATED_OPTIMAL_STOPPING
		, right_explicit(right_explicit)
		#endif
	{
		assert(tolerance > 0);
	}

	virtual Matrix A(Real t) {
		assert(t == nextTime());

		#ifdef QUANT_PDE_MODULES_HJBQVI_ITERATED_OPTIMAL_STOPPING
		if(right_explicit) {
			return Q * lA + P;
		}
		#endif

		return Q * lA + P * rA;
	}

	virtual Vector b(Real t) {
		assert(t == nextTime());

		#ifdef QUANT_PDE_MODULES_HJBQVI_ITERATED_OPTIMAL_STOPPING
		if(right_explicit) {
			Matrix M = -(rA - domain->identity());
			return Q * lb + P * (rb + M * iterand(0));
		}
		#endif

		return Q * lb + P * rb;
	}

	// Note: std::vector<bool> optimizes for space like std::bitset

	/**
	 * Creates a mask indicating where the right constraint is active.
	 * @return An std::vector of bools.
	 */
	std::vector<bool> constraintMask() {
		std::vector<bool> mask;
		mask.reserve(domain->size());

		auto v = rA * iterand(0) - rb;

		for(Index i = 0; i < domain->size(); ++i) {
			// Active constraint should be slightly negative
			mask.push_back( v(i) < 0. );
		}

		return mask;
	}

};

////////////////////////////////////////////////////////////////////////////////

/**
 * Used to solve a problem of the form
 * \f$\min\left( LV, V - V_0 \right)=0\f$ where V_0 is a function of time and
 * space.
 *
 * There are two ways to pass V_0:
 * \code{.cpp}
 * [] (Real time, Real x1, Real x2, ..., Real xDimension) {
 * 	// Return the value of the function at this point
 * }
 * \endcode
 * and
 * \code{.cpp}
 * [] (Real x1, Real x2, ..., Real xDimension) {
 * 	// Return the value of the function at this point
 * }
 * \endcode
 * The second form is independent of time.
 * @tparam Dimension The spatial dimension.
 */
template <Index Dimension>
class PenaltyMethodDifference : public PenaltyMethod {

	class DifferenceSystem : public LinearSystem {

		typedef Function<Dimension + 1> F;

		////////////////////////////////////////////////////////////////

		template <int ...Indices>
		inline Real packAndCall(
			Real time,
			const Real *array,
			Sequence<Indices...>
		) const {
			return function(time, array[Indices]...);
		}

		////////////////////////////////////////////////////////////////

		template <Index N, typename T, typename ...Ts>
		struct TimeWrapper : public TimeWrapper<N - 1, T, T, Ts...> {
		};

		template <typename T, typename ...Ts>
		struct TimeWrapper<0, T, Ts...> {
			static_assert(
				sizeof...(Ts) == Dimension,
				"The number of binding arguments in an event "
				"should be equal to the dimension"
			);

			template <typename F1>
			static inline Function<Dimension + 1> function(
					F1 &&function) {
				return [ QUANT_PDE_MOVE_CAPTURE(F1, function) ]
						(Real t, Ts ...coordinates) {
					return function(coordinates...);
				};
			}
		};

		////////////////////////////////////////////////////////////////

		const Domain<Dimension> *domain;
		F function;

	public:

		template <typename D>
		DifferenceSystem(
			D &domain,
			const Function<Dimension> &function
		) noexcept :
			domain(&domain),
			function( TimeWrapper<Dimension, Real>::function(
					function) )
		{
		}

		template <typename D>
		DifferenceSystem(
			D &domain,
			Function<Dimension> &&function
		) noexcept :
			domain(&domain),
			function( TimeWrapper<Dimension, Real>::function(
					std::move(function)) )
		{
		}

		template <typename D>
		DifferenceSystem(
			D &domain,
			const Function<Dimension + 1> &function
		) noexcept :
			domain(&domain),
			function(function)
		{
		}

		template <typename D>
		DifferenceSystem(
			D &domain,
			Function<Dimension + 1> &&function
		) noexcept :
			domain(&domain),
			function( std::move(function) )
		{
		}

		virtual Matrix A(Real) {
			return domain->identity();
		}

		virtual Vector b(Real t) {
			Vector v = domain->vector();
			for(auto node : accessor(*domain, v)) {
				*node = packAndCall(
					t,
					(&node).data(),
					GenerateSequence<Dimension>()
				);
			}
			return v;
		}

	};

	DifferenceSystem penalty;

public:

	template <typename D, typename F>
	PenaltyMethodDifference(
		D &domain,
		LinearSystem &constraint,
		F &&function,
		Real tolerance = QuantPDE::tolerance,
		bool direct = false
	) noexcept :
		PenaltyMethod(
			domain,
			constraint,
			penalty,
			tolerance,
			direct
		),
		penalty(
			domain,
			std::forward<F>(function)
		)
	{
	}

};

typedef PenaltyMethodDifference<1> PenaltyMethodDifference1;
typedef PenaltyMethodDifference<2> PenaltyMethodDifference2;
typedef PenaltyMethodDifference<3> PenaltyMethodDifference3;

} // QuantPDE

#endif
