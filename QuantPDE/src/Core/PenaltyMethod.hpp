#ifndef QUANT_PDE_CORE_PENALTY_METHOD_HPP
#define QUANT_PDE_CORE_PENALTY_METHOD_HPP

namespace QuantPDE {

/**
 * Used to solve a problem of the form
 * \f$\min\left( LV, L^\prime V \right)=0\f$.
 */
class PenaltyMethod : public IterationNode {

	DomainBase *domain;
	LinearSystem *left, *right;
	Real large;
	Matrix P;

	Matrix rA;
	Vector rb;

	virtual void onIterationStart() {
		// Only need to do this once
		rA = right->A( nextTime() );
		rb = right->b( nextTime() );

		// Evaluate the predicate using the previous iterand
		Vector predicate = rA * iterand(0) - rb;

		// Build penalty matrix
		P.setZero();
		for(Index i = 0; i < domain->size(); ++i) {
			if( predicate(i) < 0. ) {
				P.insert(i, i) = large;
			}
		}
	}

public:

	template <typename D>
	PenaltyMethod(
		D &domain,
		LinearSystem &constraint, LinearSystem &penalizedConstraint,
		Real tolerance = QuantPDE::tolerance
	) noexcept :
		domain(&domain),
		left(&constraint), right(&penalizedConstraint),
		large( 1. / tolerance ),
		P( domain.size(), domain.size() ) {
		P.reserve( IntegerVector::Constant( domain.size(), 1 ) );
		assert(tolerance > 0);
	}

	virtual Matrix A(Real t) {
		return left->A(t) + P * rA;
	}

	virtual Vector b(Real t) {
		return left->b(t) + P * rb;
	}

	// TODO: Explicit discretizations

};

////////////////////////////////////////////////////////////////////////////////

/**
 * Used to solve a problem of the form
 * \f$\min\left( LV, V - V_0 \right)=0\f$ where V_0 is a function of time and space.
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
		LinearSystem &constraint, F &&function,
		Real tolerance = QuantPDE::tolerance
	) noexcept :
		PenaltyMethod( domain, constraint, penalty, tolerance ),
		penalty( domain, std::forward<F>(function) ) {
	}

};

typedef PenaltyMethodDifference<1> PenaltyMethodDifference1;
typedef PenaltyMethodDifference<2> PenaltyMethodDifference2;
typedef PenaltyMethodDifference<3> PenaltyMethodDifference3;

} // QuantPDE

#endif

