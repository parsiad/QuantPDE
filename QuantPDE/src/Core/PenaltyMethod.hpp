#ifndef QUANT_PDE_CORE_PENALTY_METHOD
#define QUANT_PDE_CORE_PENALTY_METHOD

namespace QuantPDE {

template <Index Dimension>
class PenaltyMethod : public Linearizer {

	template <typename R, typename ...Ts>
	using Target = R (const Domain<Dimension> &, Real, Ts...);

	typedef std::function<
		Metafunctions::NaryFunctionSignatureHelpers::Type<
			Target,
			bool,
			Dimension, Real
		>
	> F;

	#define QUANT_PDE_TMP predicate( domain, time, array[Indices]... )
	template <int ...Indices>
	static inline auto packAndCall(
		const F &predicate,
		const Domain<Dimension> &domain,
		Real time,
		const Real *array,
		Metafunctions::GenerateSequenceHelpers::Sequence<Indices...>
	) -> decltype( QUANT_PDE_TMP ) {
		return QUANT_PDE_TMP;
	}
	#undef QUANT_PDE_TMP

	#define QUANT_PDE_TMP packAndCall( predicate, array, \
			GenerateSequence<N>() )
	template <int N>
	static inline auto packAndCall(
		const F &predicate,
		const Domain<Dimension> &domain,
		Real time,
		const Real *array
	) -> decltype( QUANT_PDE_TMP ) {
		return QUANT_PDE_TMP;
	}
	#undef QUANT_PDE_TMP

	const Domain<Dimension> *domain;
	Linearizer *left, *right;

	F predicate;
	Real large;
	Matrix P;

public:

	template <typename D, typename F1>
	PenaltyMethod(
		D &domain,
		Linearizer &left, Linearizer &right,
		F1 &&predicate,
		Real tolerance = 1e-6
	) noexcept :
		domain(&domain),
		left(&left), right(&right),
		predicate( std::forward<F1>(predicate) ),
		large( 1. / tolerance )
	{
		assert(large > 0);
	}

	virtual void onIterationStart() {
		P = domain->matrix();

		for(Index i = 0; i < domain->size(); i++) {

			bool penalize = PenaltyMethod::packAndCall<Dimension>(
				predicate,
				*domain,
				this->nextTime(),
				domain->coordinates(i)
			);

			if(penalize) {
				P.insert(i, i) = large;
			}
		}
	}

	virtual Matrix A() {
		return left->A() + P * right->A();
	}

	virtual Vector b() {
		return left->b() + P * right->b();
	}

};

}

#endif

