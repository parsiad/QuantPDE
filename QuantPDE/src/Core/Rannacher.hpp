#ifndef QUANT_PDE_CORE_RANNACHER_HPP
#define QUANT_PDE_CORE_RANNACHER_HPP

namespace QuantPDE {

template <Index Dimension, bool Forward>
class Rannacher : public Discretization<Dimension> {

	LinearSystem &op;

	bool   (Rannacher::*_isATheSame)() const;
	Matrix (Rannacher::*_A)(Real);
	Vector (Rannacher::*_b)(Real);
	void   (Rannacher::*_onIterationEnd)();

	inline Real difference(Real t1, Real t0) const {
		const Real dt = Forward ? t1 - t0 : t0 - t1;
		assert(dt > epsilon);
		return dt;
	}

	bool _isATheSame1() const {
		return this->isTimestepTheSame() && op.isATheSame();
	}

	bool _isATheSame2() const {
		// Switching from fully-implicit to Crank-Nicolson
		return false;
	}

	Matrix _A1(Real t1) {
		const Real t0 = this->time(0);
		const Real h0 = difference(t1, t0);

		return
			this->domain.identity()
			+ op.A(t1) * h0;
	}

	Vector _b1(Real t1) {
		const Real t0 = this->time(0);
		const Real h0 = difference(t1, t0);

		const Vector &v0 = this->iterand(0);

		return v0 + op.b(t1) * h0;
	}

	Matrix _A2(Real t1) {
		const Real t0 = this->time(0);
		const Real h0 = difference(t1, t0);

		return
			this->domain.identity()
			+ op.A(t1) * h0 / 2.
		;
	}

	Vector _b2(Real t1) {
		const Real t0 = this->time(0);
		const Real h0 = difference(t1, t0);

		const Vector &v0 = this->iterand(0);

		return (
			this->domain.identity()
			- op.A(t0) * h0 / 2.
		) * v0 + ( op.b(t1) + op.b(t0) ) / 2.;
	}

	void _onIterationEnd1() {
		_onIterationEnd = &Rannacher::_onIterationEnd2;
	}

	void _onIterationEnd2() {
		_isATheSame = &Rannacher::_isATheSame2;
		_A = &Rannacher::_A2;
		_b = &Rannacher::_b2;
		_onIterationEnd = &Rannacher::_onIterationEnd3;
	}

	void _onIterationEnd3() {
		_isATheSame = &Rannacher::_isATheSame1;
		_onIterationEnd = &Rannacher::_onIterationEnd4;
	}

	void _onIterationEnd4() {
	}

	virtual Matrix A(Real t) {
		return (this->*_A)(t);
	}

	virtual Vector bd(Real t) {
		return (this->*_b)(t);
	}

public:

	virtual bool isATheSame() const {
		return (this->*_isATheSame)();
	}

	virtual void clear() {
		_isATheSame = &Rannacher::_isATheSame1;
		_A = &Rannacher::_A1;
		_b = &Rannacher::_b1;
		_onIterationEnd = &Rannacher::_onIterationEnd1;
	}

	virtual void onAfterEvent() {
		// Do nothing; assume smoothness is preserved
	}

	virtual void onIterationEnd() {
		(this->*_onIterationEnd)();
	}

	template <typename D>
	Rannacher(
		D &domain,
		LinearSystem &op
	) noexcept :
		Discretization<Dimension>(domain),
		op(op),
		_isATheSame(nullptr),
		_A(nullptr),
		_b(nullptr),
		_onIterationEnd(nullptr)
	{
	}

};

template <Index Dimension>
using ReverseRannacher = Rannacher<Dimension, false>;

template <Index Dimension>
using ForwardRannacher = Rannacher<Dimension, true>;

typedef ReverseRannacher<1> ReverseRannacher1;
typedef ReverseRannacher<2> ReverseRannacher2;
typedef ReverseRannacher<3> ReverseRannacher3;

typedef ForwardRannacher<1> ForwardRannacher1;
typedef ForwardRannacher<2> ForwardRannacher2;
typedef ForwardRannacher<3> ForwardRannacher3;

} // QuantPDE

#endif

