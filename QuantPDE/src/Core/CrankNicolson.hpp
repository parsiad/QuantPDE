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

	virtual bool isAConstant() const {
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

template <bool Forward>
class Rannacher : public Linearizer {

	const DomainBase *domain;
	const LinearOperator *op;

	bool   (Rannacher::*_isAConstant)() const;
	Matrix (Rannacher::*_A )();
	Vector (Rannacher::*_b )();
	void   (Rannacher::*_onIterationEnd)();

	Real t1, t0, h0;

	inline Real difference(Real t1, Real t0) {
		return Forward ? t1 - t0 : t0 - t1;
	}

	bool _isAConstant1() const {
		return op->isConstantInTime();
	}

	bool _isAConstant2() const {
		return true;
	}

	Matrix _A1() {
		t1 = this->nextTime();
		t0 = this->times()[0];

		h0 = difference(t1, t0);

		return
			domain->identity()
			+ op->discretize(t1) * h0;
	}

	Vector _b1() {
		const Vector &v0 = this->iterands()[0];
		return v0;
	}

	Matrix _A2() {
		t1 = this->nextTime();
		t0 = this->times()[0];

		return
			domain->identity()
			+ op->discretize(t1) * h0 / 2.
		;
	}

	Vector _b2() {
		const Vector &v0 = this->iterands()[0];

		return (
			domain->identity()
			- op->discretize(t0) * h0 / 2.
		) * v0;
	}

	void _onIterationEnd1() {
		_isAConstant = &Rannacher::_isAConstant2;
		_onIterationEnd = &Rannacher::_onIterationEnd2;
	}

	void _onIterationEnd2() {
		_isAConstant = &Rannacher::_isAConstant1;
		_A = &Rannacher::_A2;
		_b = &Rannacher::_b2;
		_onIterationEnd = &Rannacher::_onIterationEnd3;
	}

	void _onIterationEnd3() {
	}

public:

	virtual bool isAConstant() const {
		return (this->*_isAConstant)();
	}

	virtual Matrix A() {
		return (this->*_A)();
	}

	virtual Vector b() {
		return (this->*_b)();
	}

	virtual void clear() {
		_isAConstant = &Rannacher::_isAConstant1;
		_A = &Rannacher::_A1;
		_b = &Rannacher::_b1;
		_onIterationEnd = &Rannacher::_onIterationEnd3;
	}

	virtual void onIterationEnd() {
		(this->*_onIterationEnd)();
	}

	template <typename D, typename L>
	Rannacher(D &domain, L &op) noexcept : domain(&domain), op(&op) {
	}

};

typedef Rannacher<false> ReverseRannacher;
typedef Rannacher<true> ForwardRannacher;

} // QuantPDE

#endif

