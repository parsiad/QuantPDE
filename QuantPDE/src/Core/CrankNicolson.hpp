#ifndef QUANT_PDE_CORE_CRANK_NICOLSON_HPP
#define QUANT_PDE_CORE_CRANK_NICOLSON_HPP

namespace QuantPDE {

template <bool Forward>
class CrankNicolson : public IterationNode {

	const DomainBase *domain;
	LinearSystem *op;

	inline Real dt() {
		const Real
			t1 = nextTime(),
			t0 = time(0)
		;

		return Forward ? t1 - t0 : t0 - t1;
	}

	virtual Matrix A(Real t1) {
		return
			domain->identity()
			+ op->A(t1) * dt() / 2.
		;
	}

	virtual Vector b(Real t1) {
		const Real    t0 = time(0);
		const Vector &v0 = iterand(0);

		return (
			domain->identity()
			- op->A(t0) * dt() / 2.
		) * v0 + ( op->b(t1) + op->b(t0) ) / 2.;
	}

public:

	virtual bool isATheSame() const {
		return isTimestepTheSame() && op->isATheSame();
	}

	template <typename D>
	CrankNicolson(D &domain, LinearSystem &op) noexcept : domain(&domain),
			op(&op) {
	}

};

typedef CrankNicolson<false> ReverseCrankNicolson;
typedef CrankNicolson<true > ForwardCrankNicolson;

} // QuantPDE

#endif

