#ifndef QUANT_PDE_CORE_PENALTY_METHOD
#define QUANT_PDE_CORE_PENALTY_METHOD

namespace QuantPDE {

/**
 * Used to solve a problem of the form
 * \f$\min( LV, L^\prime V \right)=0\f$.
 */
class PenaltyMethod : public LinearSystemIteration {

	DomainBase *domain;
	LinearSystem *left, *right;
	Real large;
	Matrix P;

	Matrix rA;
	Vector rb;

	virtual void onIterationStart() {
		// Only need to do this once
		rA = right->A( this->nextTime() );
		rb = right->b( this->nextTime() );

		// Evaluate the predicate using the previous iterand
		Vector predicate = rA * this->iterands()[0] - rb;

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

};

////////////////////////////////////////////////////////////////////////////////

/**
 * Used to solve a problem of the form
 * \f$\min( LV, V - V_0 \right)=0\f$ where V_0 is a function of time and space.
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

		typedef Function<Dimension + 1> Difference;

		template <int ...Indices>
		static inline Real packAndCall(
			const Difference &predicate,
			//const Domain<Dimension> &domain,
			//const Vector &iterand,
			Real time,
			const Real *array,
			Metafunctions::GenerateSequenceHelpers::Sequence<
					Indices...>
		) {
			return predicate( time, array[Indices]... );
		}

		template <int N>
		static inline Real packAndCall(
			const Difference &predicate,
			//const Domain<Dimension> &domain,
			//const Vector &iterand,
			Real time,
			const Real *array
		) {
			return DifferenceSystem::packAndCall( predicate,
					time, array, GenerateSequence<N>() );
		}

		const PenaltyMethodDifference *parent;
		const Domain<Dimension> *domain;
		Difference function;

	public:

		template <typename P, typename D, typename F>
		DifferenceSystem(P &parent, D &domain, F &&function) noexcept
				: parent(&parent), domain(&domain),
				function( std::forward<F>(function) ) {
		}

		virtual Matrix A(Real) {
			return domain->identity();
		}

		virtual Vector b(Real) {
			Vector v = domain->vector();
			for(auto node : domain->accessor(v)) {
				*node = DifferenceSystem::packAndCall<
						Dimension>(
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

