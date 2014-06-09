#ifndef QUANT_PDE_CORE_CRANK_NICOLSON
#define QUANT_PDE_CORE_CRANK_NICOLSON

namespace QuantPDE {

template <bool Forward>
class CrankNicolson : public Linearizer {

	const DomainBase *domain;
	const LinearOperator *op;

	inline Real dt() {
		const Real
			t1 = this->nextTime(),
			t0 = this->times()[0]
		;

		return Forward ? t1 - t0 : t0 - t1;
	}

public:

	virtual bool doesAChange() const {
		return op->isConstantInTime();
	}

	virtual Matrix A() {
		const Real t1 = this->nextTime();

		return
			domain->identity()
			+ op->discretize(t1) * dt() / 2.
		;
	}

	virtual Vector b() {
		const Real    t0 = this->times()[0];
		const Vector &v0 = this->iterands()[0];

		return (
			domain->identity()
			- op->discretize(t0) * dt() / 2.
		) * v0;
	}

	template <typename D, typename L>
	CrankNicolson(D &domain, L &op) noexcept : domain(&domain), op(&op) {
	}

};

typedef CrankNicolson<false> ReverseCrankNicolson;
typedef CrankNicolson<true > ForwardCrankNicolson;

} // QuantPDE

#endif

