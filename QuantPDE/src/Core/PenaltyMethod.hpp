#ifndef QUANT_PDE_CORE_PENALTY_METHOD
#define QUANT_PDE_CORE_PENALTY_METHOD

namespace QuantPDE {

/*
template <Index Dimension>
class PredicatePenaltyMethod : public Linearizer {

	template <typename R, typename ...Ts>
	using Target = R (const Domain<Dimension> &, const Vector &iterand,
			Real, Ts...);

	typedef std::function<
		Metafunctions::NaryFunctionSignatureHelpers::Type<
			Target,
			bool,
			Dimension, Real
		>
	> Predicate;

	#define QUANT_PDE_TMP predicate( domain, iterand, time, \
			array[Indices]... )
	template <int ...Indices>
	static inline auto packAndCall(
		const Predicate &predicate,
		const Domain<Dimension> &domain,
		const Vector &iterand,
		Real time,
		const Real *array,
		Metafunctions::GenerateSequenceHelpers::Sequence<Indices...>
	) -> decltype( QUANT_PDE_TMP ) {
		return QUANT_PDE_TMP;
	}
	#undef QUANT_PDE_TMP

	#define QUANT_PDE_TMP packAndCall( predicate, domain, iterand, time, \
			array, GenerateSequence<N>() )
	template <int N>
	static inline auto packAndCall(
		const Predicate &predicate,
		const Domain<Dimension> &domain,
		const Vector &iterand,
		Real time,
		const Real *array
	) -> decltype( QUANT_PDE_TMP ) {
		return QUANT_PDE_TMP;
	}
	#undef QUANT_PDE_TMP

	const Domain<Dimension> *domain;
	Linearizer *left, *right;

	Predicate predicate;
	Real large;
	Matrix P;

	virtual void onIterationStart() {
		// Build penalty matrix using predicate
		// TODO: Initialize P

		for(Index i = 0; i < domain->size(); i++) {
			if( PenaltyMethod::packAndCall<Dimension>(
				predicate,
				*domain,
				this->iterands()[0], // Most recent iterand
				this->nextTime(), // Current time
				domain->coordinates(i).data() // Coordinates
			) ) {
				// Predicate evaluated to true; need to penalize
				P.insert(i, i) = large;
			}
		}
	}

protected:

	virtual bool doesAChange() const {
		return left->doesAChange() || right->doesAChange();
	}

	virtual Matrix A() {
		return left->A() + P * right->A();
	}

	virtual Vector b() {
		return left->b() + P * right->b();
	}

public:

	template <typename D, typename F1>
	PenaltyMethod(
		D &domain,
		Linearizer &constraint, Linearizer &penalizedConstraint,
		F1 &&predicate,
		Real tolerance = 1e-6
	) noexcept :
		domain(&domain),
		left(&left), right(&penalizedConstraint),
		predicate( std::forward<F1>(predicate) ),
		large( 1. / tolerance ),
		P( domain.size(), domain.size() )
	{
		assert(large > 0);
	}

};
*/

/**
 * Used to solve a problem of the form
 * \f$\min( LV, L^\prime V \right)=0\f$.
 */
class PenaltyMethod : public Linearizer {

	DomainBase *domain;
	Linearizer *left, *right;
	Real large;
	Matrix P;

	Matrix rA;
	Vector rb;

	virtual void onIterationStart() {
		// Only need to do this once
		rA = right->A();
		rb = right->b();

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
		Linearizer &constraint, Linearizer &penalizedConstraint,
		Real tolerance = 1e-6
	) noexcept :
		domain(&domain),
		left(&constraint), right(&penalizedConstraint),
		large( 1. / tolerance ),
		P( domain.size(), domain.size() )
	{
		P.reserve( IntegerVector::Constant( domain.size(), 1 ) );
		assert(tolerance > 0);
	}

	virtual bool doesAChange() const {
		// Presumably, yes, since the penalty matrix will most likely
		// change (otherwise, why are we using a penalty method?)
		return true;
	}

	virtual Matrix A() {
		return left->A() + P * rA;
	}

	virtual Vector b() {
		return left->b() + P * rb;
	}

};

/**
 * Used to solve a problem of the form
 * \f$\min( LV, V - V_0 \right)=0\f$ where V_0 is a simple function of space
 * alone.
 */
template <Index Dimension>
class SimplePenaltyMethod : public PenaltyMethod {

	class SimplePenaltyConstraint : public Linearizer {

		Domain<Dimension> *domain;
		Function<Dimension> function;

	public:

		template <typename D, typename F>
		SimplePenaltyConstraint(D &domain, F &&function) noexcept
				: domain(&domain),
				function( std::forward<F>(function) ) {
		}

		virtual Matrix A() {
			return domain->identity();
		}

		virtual Vector b() {
			return domain->image(function);
		}

	};

	SimplePenaltyConstraint penaltyConstraint;

public:

	template <typename D, typename F>
	SimplePenaltyMethod(
		D &domain,
		Linearizer &constraint, F &&function,
		Real tolerance = 1e-6
	) noexcept :
		PenaltyMethod(domain, constraint, penaltyConstraint, tolerance),
		penaltyConstraint( domain, std::forward<F>(function) )
	{
	}

};

typedef SimplePenaltyMethod<1> SimplePenaltyMethod1;
typedef SimplePenaltyMethod<2> SimplePenaltyMethod2;
typedef SimplePenaltyMethod<3> SimplePenaltyMethod3;

} // QuantPDE

#endif

