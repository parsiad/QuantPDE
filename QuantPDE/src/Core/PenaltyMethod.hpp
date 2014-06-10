#ifndef QUANT_PDE_CORE_PENALTY_METHOD
#define QUANT_PDE_CORE_PENALTY_METHOD

namespace QuantPDE {

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
		P( domain.size(), domain.size() ) {
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

////////////////////////////////////////////////////////////////////////////////

/**
 * Used to solve a problem of the form
 * \f$\min( LV, V - V_0 \right)=0\f$ where V_0 is a function of time and space.
 */
template <Index Dimension>
class PenaltyMethodDifference : public PenaltyMethod {

	class DifferenceLinearizer : public Linearizer {

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
			return DifferenceLinearizer::packAndCall( predicate,
					time, array, GenerateSequence<N>() );
		}

		const Domain<Dimension> *domain;
		Difference function;

	public:

		template <typename D, typename F>
		DifferenceLinearizer(D &domain, F &&function) noexcept
				: domain(&domain),
				function( std::forward<F>(function) ) {
		}

		virtual Matrix A() {
			return domain->identity();
		}

		virtual Vector b() {
			Vector v = domain->vector();
			for(auto node : domain->accessor(v)) {
				*node = DifferenceLinearizer::packAndCall<
						Dimension>(
					function,
					//*domain,
					//this->iterands()[0],
					this->nextTime(),
					(&node).data()
				);
			}
			return v;
		}

	};

	DifferenceLinearizer linearizer;

public:

	template <typename D, typename F>
	PenaltyMethodDifference(
		D &domain,
		Linearizer &constraint, F &&function,
		Real tolerance = 1e-6
	) noexcept :
		PenaltyMethod( domain, constraint, linearizer, tolerance ),
		linearizer( domain, std::forward<F>(function) ) {
	}

	virtual void setIteration(Iteration &iteration) {
		PenaltyMethod::setIteration(iteration);
		linearizer.setIteration(iteration);
	}

};

typedef PenaltyMethodDifference<1> PenaltyMethodDifference1;
typedef PenaltyMethodDifference<2> PenaltyMethodDifference2;
typedef PenaltyMethodDifference<3> PenaltyMethodDifference3;

} // QuantPDE

#endif

