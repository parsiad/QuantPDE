#ifndef QUANT_PDE_CORE_PENALTY_METHOD
#define QUANT_PDE_CORE_PENALTY_METHOD

// TODO: Cache A if constant

namespace QuantPDE {

/**
 * Used to solve a problem of the form
 * \f$\min( LV, L^\prime V \right)=0\f$.
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
		for(Index i = 0; i < domain->size(); i++) {
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
		Real tolerance = 1e-6
	) noexcept :
		domain(&domain),
		left(&constraint), right(&penalizedConstraint),
		large( 1. / tolerance ),
		P( domain.size(), domain.size() ) {
		P.reserve( IntegerVector::Constant( domain.size(), 1 ) );
		assert(tolerance > 0);
	}

	virtual Matrix A() {
		return left->A( nextTime() ) + P * rA;
	}

	virtual Vector b() {
		return left->b( nextTime() ) + P * rb;
	}

	virtual int minimumLookback() const {
		return 1;
	}

	// TODO: Explicit discretizations

};

////////////////////////////////////////////////////////////////////////////////

/**
 * Used to solve a problem of the form
 * \f$\min( LV, V - V_0 \right)=0\f$ where V_0 is a function of time and space.
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

		/*
		template <typename R, typename ...Ts>
		using Target = R (const Domain<Dimension> &, const Vector &,
				Real, Ts...);

		typedef std::function<
			Metafunctions::NaryFunctionSignatureHelpers::Type<
				Target,
				Real,
				Dimension, Real
			>
		> Predicate;
		*/

		typedef Function<Dimension + 1> F;

		////////////////////////////////////////////////////////////////

		template <int ...Indices>
		static inline Real packAndCall(
			const F &V_0,
			//const Domain<Dimension> &domain,
			//const Vector &iterand,
			Real time,
			const Real *array,
			Metafunctions::Sequence<Indices...>
		) {
			return V_0( time, array[Indices]... );
		}

		template <int N>
		static inline Real packAndCall(
			const F &V_0,
			//const Domain<Dimension> &domain,
			//const Vector &iterand,
			Real time,
			const Real *array
		) {
			return packAndCall(
				V_0,
				time,
				array,
				GenerateSequence<N>()
			);
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
				F1 &&function
			) {
				// TODO: Use C++1y move-capture in the future;
				//       currently we perform an unconditional
				//       capture-by-value

				// [ function(std::forward<F1>(function)) ]

				return [function] (Real t, Ts ...coordinates) {
					return function(coordinates...);
				};
			}
		};

		////////////////////////////////////////////////////////////////

		const PenaltyMethodDifference *parent;
		const Domain<Dimension> *domain;
		F function;

	public:

		template <typename P, typename D>
		DifferenceSystem(
			P &parent,
			D &domain,
			const Function<Dimension> &function
		) noexcept :
			parent(&parent),
			domain(&domain),
			function( TimeWrapper<Dimension, Real>::function(
					function) )
		{
		}

		template <typename P, typename D>
		DifferenceSystem(
			P &parent,
			D &domain,
			Function<Dimension> &&function
		) noexcept :
			parent(&parent),
			domain(&domain),
			function( TimeWrapper<Dimension, Real>::function(
					std::move(function)) )
		{
		}

		template <typename P, typename D>
		DifferenceSystem(
			P &parent,
			D &domain,
			const Function<Dimension + 1> &function
		) noexcept :
			parent(&parent),
			domain(&domain),
			function(function)
		{
		}

		template <typename P, typename D>
		DifferenceSystem(
			P &parent,
			D &domain,
			Function<Dimension + 1> &&function
		) noexcept :
			parent(&parent),
			domain(&domain),
			function( std::move(function) )
		{
		}

		virtual Matrix A(Real) {
			return domain->identity();
		}

		virtual Vector b(Real) {
			Vector v = domain->vector();
			for(auto node : accessor(*domain, v)) {
				*node = packAndCall<Dimension> (
					function,
					parent->nextTime(),
					(&node).data()
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
		Real tolerance = 1e-6
	) noexcept :
		PenaltyMethod( domain, constraint, penalty, tolerance ),
		penalty( *this, domain, std::forward<F>(function) ) {
	}

};

typedef PenaltyMethodDifference<1> PenaltyMethodDifference1;
typedef PenaltyMethodDifference<2> PenaltyMethodDifference2;
typedef PenaltyMethodDifference<3> PenaltyMethodDifference3;

} // QuantPDE

#endif

