#ifndef QUANT_PDE_CORE_LINEAR_BDF
#define QUANT_PDE_CORE_LINEAR_BDF

namespace QuantPDE {

template <size_t Lookback, bool Forward>
class LinearBDFBase : public Linearizer<Lookback> {

	const DomainBase *domain;
	const LinearOperator *op;

	inline Real dt() {
		const Real
			t0 = std::get<0>( this->iterands()[0] ),
			t1 = this->nextTime()
		;

		return Forward ? t1 - t0 : t0 - t1;
	}

protected:

	inline Matrix A1() {
		const Real t = this->nextTime();

		return
			domain->identity()
			+ op->discretize(t) * dt();
	}

	inline Vector b1() {
		const Vector &v_n0 = std::get<1>( this->iterands()[0] );
		return v_n0;
	}

	inline Matrix A2() {
		const Real t = this->nextTime();

		return
			domain->identity()
			+ 2. / 3. * op->discretize(t) * dt();
	}

	inline Vector b2() {
		const Vector
			&v_n1 = std::get<1>( this->iterands()[ 0] ),
			&v_n0 = std::get<1>( this->iterands()[-1] )
		;

		return (
			  4. * v_n1
			- 1. * v_n0
		) / 3.;
	}

	inline Matrix A3() {
		const Real t = this->nextTime();

		return
			domain->identity()
			+ 6. / 11. * op->discretize(t) * dt();
	}

	inline Vector b3() {
		const Vector
			&v_n2 = std::get<1>( this->iterands()[ 0] ),
			&v_n1 = std::get<1>( this->iterands()[-1] ),
			&v_n0 = std::get<1>( this->iterands()[-2] )
		;

		return (
			  18. * v_n2
			- 9.  * v_n1
			+ 2.  * v_n0
		) / 11.;
	}

	inline Matrix A4() {
		const Real t = this->nextTime();
		return
			domain->identity()
			+ 12. / 25. * op->discretize(t) * dt();
	}

	inline Vector b4() {
		const Vector
			&v_n3 = std::get<1>( this->iterands()[ 0] ),
			&v_n2 = std::get<1>( this->iterands()[-1] ),
			&v_n1 = std::get<1>( this->iterands()[-2] ),
			&v_n0 = std::get<1>( this->iterands()[-3] )
		;

		return (
			  48. * v_n3
			- 36. * v_n2
			+ 16. * v_n1
			- 3.  * v_n0
		) / 25.;
	}

	inline Matrix A5() {
		const Real t = this->nextTime();

		return
			domain->identity()
			+ 60. / 137. * op->discretize(t) * dt();
	}

	inline Vector b5() {
		const Vector
			&v_n4 = std::get<1>( this->iterands()[ 0] ),
			&v_n3 = std::get<1>( this->iterands()[-1] ),
			&v_n2 = std::get<1>( this->iterands()[-2] ),
			&v_n1 = std::get<1>( this->iterands()[-3] ),
			&v_n0 = std::get<1>( this->iterands()[-4] )
		;

		return (
			  300. * v_n4
			- 300. * v_n3
			+ 200. * v_n2
			- 75.  * v_n1
			+ 12.  * v_n0
		) / 137.;
	}

	inline Matrix A6() {
		const Real t = this->nextTime();

		return
			domain->identity()
			+ 60. / 147. * op->discretize(t) * dt();
	}

	inline Vector b6() {
		const Vector
			&v_n5 = std::get<1>( this->iterands()[ 0] ),
			&v_n4 = std::get<1>( this->iterands()[-1] ),
			&v_n3 = std::get<1>( this->iterands()[-2] ),
			&v_n2 = std::get<1>( this->iterands()[-3] ),
			&v_n1 = std::get<1>( this->iterands()[-4] ),
			&v_n0 = std::get<1>( this->iterands()[-5] )
		;

		return (
			  360. * v_n5
			- 450. * v_n4
			+ 400. * v_n3
			- 225. * v_n2
			+ 72.  * v_n1
			- 10.  * v_n0
		) / 147.;
	}

public:

	template <typename D, typename L>
	LinearBDFBase(D &domain, L &op) noexcept
			: domain(&domain), op(&op) {
	}

};

template <bool Forward = false>
class LinearBDFOne : public LinearBDFBase<1, Forward> {

public:

	template <typename D, typename L>
	LinearBDFOne(D &domain, L &op) noexcept
			: LinearBDFBase<1, Forward>(domain, op) {
	}

	virtual Matrix A() {
		return this->A1();
	}

	virtual Vector b() {
		return this->b1();
	}

};

// Convenient alias
template <bool Forward = false>
using ImplicitMethod = LinearBDFOne<Forward>;

////////////////////////////////////////////////////////////////////////////////

template <bool Forward = false>
class LinearBDFTwo : public LinearBDFBase<2, Forward> {

	Matrix (LinearBDFTwo<Forward>::*AA)();
	Vector (LinearBDFTwo<Forward>::*bb)();

	////////////////////////////////////////////////////////////////////////

	Matrix _A1() {
		AA = &LinearBDFTwo::A2;
		return this->A1();
	}

	Vector _b1() {
		bb = &LinearBDFTwo::b2;
		return this->b1();
	}

public:

	template <typename D, typename L>
	LinearBDFTwo(D &domain, L &op) noexcept
			: LinearBDFBase<2, Forward>(domain, op) {
	}

	virtual void clear() {
		AA = &LinearBDFTwo::_A1;
		bb = &LinearBDFTwo::_b1;
	}

	virtual Matrix A() {
		return (this->*AA)();
	}

	virtual Vector b() {
		return (this->*bb)();
	}

};

////////////////////////////////////////////////////////////////////////////////

template <bool Forward = false>
class LinearBDFThree : public LinearBDFBase<3, Forward> {

	Matrix (LinearBDFThree<Forward>::*AA)();
	Vector (LinearBDFThree<Forward>::*bb)();

	////////////////////////////////////////////////////////////////////////

	Matrix _A1() {
		AA = &LinearBDFThree::_A2;
		return this->A1();
	}

	Vector _b1() {
		bb = &LinearBDFThree::_b2;
		return this->b1();
	}

	Matrix _A2() {
		AA = &LinearBDFThree::A3;
		return this->A2();
	}

	Vector _b2() {
		bb = &LinearBDFThree::b3;
		return this->b2();
	}

public:

	template <typename D, typename L>
	LinearBDFThree(D &domain, L &op) noexcept
			: LinearBDFBase<3, Forward>(domain, op) {
	}

	virtual void clear() {
		AA = &LinearBDFThree::_A1;
		bb = &LinearBDFThree::_b1;
	}

	virtual Matrix A() {
		return (this->*AA)();
	}

	virtual Vector b() {
		return (this->*bb)();
	}

};

////////////////////////////////////////////////////////////////////////////////

template <bool Forward = false>
class LinearBDFFour : public LinearBDFBase<4, Forward> {

	Matrix (LinearBDFFour<Forward>::*AA)();
	Vector (LinearBDFFour<Forward>::*bb)();

	////////////////////////////////////////////////////////////////////////

	Matrix _A1() {
		AA = &LinearBDFFour::_A2;
		return this->A1();
	}

	Vector _b1() {
		bb = &LinearBDFFour::_b2;
		return this->b1();
	}

	Matrix _A2() {
		AA = &LinearBDFFour::A3;
		return this->A2();
	}

	Vector _b2() {
		bb = &LinearBDFFour::b3;
		return this->b2();
	}

	Matrix _A3() {
		AA = &LinearBDFFour::A4;
		return this->A3();
	}

	Vector _b3() {
		bb = &LinearBDFFour::b4;
		return this->b3();
	}

public:

	template <typename D, typename L>
	LinearBDFFour(D &domain, L &op) noexcept
			: LinearBDFBase<4, Forward>(domain, op) {
	}

	virtual void clear() {
		AA = &LinearBDFFour::_A1;
		bb = &LinearBDFFour::_b1;
	}

	virtual Matrix A() {
		return (this->*AA)();
	}

	virtual Vector b() {
		return (this->*bb)();
	}

};

////////////////////////////////////////////////////////////////////////////////

template <bool Forward = false>
class LinearBDFFive : public LinearBDFBase<5, Forward> {

	Matrix (LinearBDFFive<Forward>::*AA)();
	Vector (LinearBDFFive<Forward>::*bb)();

	////////////////////////////////////////////////////////////////////////

	Matrix _A1() {
		AA = &LinearBDFFive::_A2;
		return this->A1();
	}

	Vector _b1() {
		bb = &LinearBDFFive::_b2;
		return this->b1();
	}

	Matrix _A2() {
		AA = &LinearBDFFive::A3;
		return this->A2();
	}

	Vector _b2() {
		bb = &LinearBDFFive::b3;
		return this->b2();
	}

	Matrix _A3() {
		AA = &LinearBDFFive::_A4;
		return this->A3();
	}

	Vector _b3() {
		bb = &LinearBDFFive::_b4;
		return this->b3();
	}

	Matrix _A4() {
		AA = &LinearBDFFive::A5;
		return this->A4();
	}

	Vector _b4() {
		bb = &LinearBDFFive::b5;
		return this->b4();
	}

public:

	template <typename D, typename L>
	LinearBDFFive(D &domain, L &op) noexcept
			: LinearBDFBase<5, Forward>(domain, op) {
	}

	virtual void clear() {
		AA = &LinearBDFFive::_A1;
		bb = &LinearBDFFive::_b1;
	}

	virtual Matrix A() {
		return (this->*AA)();
	}

	virtual Vector b() {
		return (this->*bb)();
	}

};

////////////////////////////////////////////////////////////////////////////////

template <bool Forward = false>
class LinearBDFSix : public LinearBDFBase<6, Forward> {

	Matrix (LinearBDFSix<Forward>::*AA)();
	Vector (LinearBDFSix<Forward>::*bb)();

	////////////////////////////////////////////////////////////////////////

	Matrix _A1() {
		AA = &LinearBDFSix::_A2;
		return this->A1();
	}

	Vector _b1() {
		bb = &LinearBDFSix::_b2;
		return this->b1();
	}

	Matrix _A2() {
		AA = &LinearBDFSix::A3;
		return this->A2();
	}

	Vector _b2() {
		bb = &LinearBDFSix::b3;
		return this->b2();
	}

	Matrix _A3() {
		AA = &LinearBDFSix::_A4;
		return this->A3();
	}

	Vector _b3() {
		bb = &LinearBDFSix::_b4;
		return this->b3();
	}

	Matrix _A4() {
		AA = &LinearBDFSix::_A5;
		return this->A4();
	}

	Vector _b4() {
		bb = &LinearBDFSix::_b5;
		return this->b4();
	}

	Matrix _A5() {
		AA = &LinearBDFSix::A6;
		return this->A5();
	}

	Vector _b5() {
		bb = &LinearBDFSix::b6;
		return this->b6();
	}

public:

	template <typename D, typename L>
	LinearBDFSix(D &domain, L &op) noexcept
			: LinearBDFBase<6, Forward>(domain, op) {
	}

	virtual void clear() {
		AA = &LinearBDFSix::_A1;
		bb = &LinearBDFSix::_b1;
	}

	virtual Matrix A() {
		return (this->*AA)();
	}

	virtual Vector b() {
		return (this->*bb)();
	}

};

} // QuantPDE

#endif

