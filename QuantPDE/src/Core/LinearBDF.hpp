#ifndef QUANT_PDE_CORE_LINEAR_BDF
#define QUANT_PDE_CORE_LINEAR_BDF

namespace QuantPDE {

template <bool Forward>
class LinearBDFBase : public Linearizer {

	const DomainBase *domain;
	const LinearOperator *op;

	Real t6, t5, t4, t3, t2, t1, t0;
	Real     h5, h4, h3, h2, h1, h0;

	inline Real difference(Real t1, Real t0) {
		return Forward ? t1 - t0 : t0 - t1;
	}

protected:

	inline Matrix A1() {
		t1 = this->nextTime();
		t0 = this->times()[0];

		h0 = difference(t1, t0);

		return
			domain->identity()
			+ op->discretize(t1) * h0;
	}

	inline Vector b1() {
		const Vector &v0 = this->iterands()[0];
		return v0;
	}

	inline Matrix A2() {
		t2 = this->nextTime();
		t1 = this->times()[ 0];
		t0 = this->times()[-1];

		h1 = difference(t2, t1);
		h0 = difference(t2, t0);

		return
			domain->identity()
			+ (h0 * h1) / (h0 + h1) * op->discretize(t2)
		;

		// Constant timestep case:
		/*
		return
			domain->identity()
			+ 2. / 3. * op->discretize(t2) * dt()
		;
		*/
	}

	inline Vector b2() {
		const Vector
			&v1 = this->iterands()[ 0],
			&v0 = this->iterands()[-1]
		;

		return (
			  h0*h0 * v1
			- h1*h1 * v0
		) / ( (h0 + h1) * (h0 - h1) );

		// Constant timestep case:
		/*
		return (
			  4. * v1
			- 1. * v0
		) / 3.;
		*/
	}

	inline Matrix A3() {
		t3 = this->nextTime();
		t2 = this->times()[ 0];
		t1 = this->times()[-1];
		t0 = this->times()[-2];

		h2 = difference(t3, t2);
		h1 = difference(t3, t1);
		h0 = difference(t3, t0);

		return
			domain->identity()
			+ ( h0 * h1 * h2 ) / ( h0 * h1 + h0 * h2 + h1 * h2 )
					* op->discretize(t3)
		;

		// Constant timestep case:
		/*
		return
			domain->identity()
			+ 6. / 11. * op->discretize(t3) * dt()
		;
		*/
	}

	inline Vector b3() {
		const Vector
			&v2 = this->iterands()[ 0],
			&v1 = this->iterands()[-1],
			&v0 = this->iterands()[-2]
		;

		return (
			  (h0*h0 * h1*h1) / ((h0 - h2) * (h1 - h2)) * v2
			- (h0*h0 * h2*h2) / ((h0 - h1) * (h1 - h2)) * v1
			+ (h1*h1 * h2*h2) / ((h0 - h1) * (h0 - h2)) * v0
		) / (h0 * h1 + h0 * h2 + h1 * h2);

		// Constant timestep case:
		/*
		return (
			  18. * v2
			- 9.  * v1
			+ 2.  * v0
		) / 11.;
		*/
	}

	inline Matrix A4() {

		t4 = this->nextTime();
		t3 = this->times()[ 0];
		t2 = this->times()[-1];
		t1 = this->times()[-2];
		t0 = this->times()[-3];

		h3 = difference(t4, t3);
		h2 = difference(t4, t2);
		h1 = difference(t4, t1);
		h0 = difference(t4, t0);

		return domain->identity()
			+ (h0*h1*h2*h3) / (h0*h1*h2 + h0*h1*h3 + h0*h2*h3 + h1*h2*h3)
					* op->discretize(t4);


		// Constant timestep case:
		/*
		return
			domain->identity()
			+ 12. / 25. * op->discretize(t4) * dt();
		*/
	}

	inline Vector b4() {
		const Vector
			&v3 = this->iterands()[ 0],
			&v2 = this->iterands()[-1],
			&v1 = this->iterands()[-2],
			&v0 = this->iterands()[-3]
		;

		return (
			  (h0*h0 * h1*h1 * h2*h2) / ((h0 - h3)*(h1 - h3)*(h2 - h3)) * v3
			- (h0*h0 * h1*h1 * h3*h3) / ((h0 - h2)*(h1 - h2)*(h2 - h3)) * v2
			+ (h0*h0 * h2*h2 * h3*h3) / ((h0 - h1)*(h1 - h2)*(h1 - h3)) * v1
			- (h1*h1 * h2*h2 * h3*h3) / ((h0 - h1)*(h0 - h2)*(h0 - h3)) * v0
		) / (h0*h1*h2 + h0*h1*h3 + h0*h2*h3 + h1*h2*h3);

		// Constant timestep case:
		/*
		return (
			  48. * v3
			- 36. * v2
			+ 16. * v1
			- 3.  * v0
		) / 25.;
		*/
	}

	inline Matrix A5() {
		t5 = this->nextTime();
		t4 = this->times()[ 0];
		t3 = this->times()[-1];
		t2 = this->times()[-2];
		t1 = this->times()[-3];
		t0 = this->times()[-4];

		h4 = difference(t5, t4);
		h3 = difference(t5, t3);
		h2 = difference(t5, t2);
		h1 = difference(t5, t1);
		h0 = difference(t5, t0);

		return
			domain->identity()
			+ (h0*h1*h2*h3*h4) / (h0*h1*h2*h3 + h0*h1*h2*h4 + h0*h1*h3*h4 + h0*h2*h3*h4 + h1*h2*h3*h4)
					* op->discretize(t5)
		;

		// Constant timestep case:
		/*
		return
			domain->identity()
			+ 60. / 137. * op->discretize(t5) * dt()
		;
		*/
	}

	inline Vector b5() {
		const Vector
			&v4 = this->iterands()[ 0],
			&v3 = this->iterands()[-1],
			&v2 = this->iterands()[-2],
			&v1 = this->iterands()[-3],
			&v0 = this->iterands()[-4]
		;

		return (
			  (h0*h0 * h1*h1 * h2*h2 * h3*h3) / ((h0 - h4)*(h1 - h4)*(h2 - h4)*(h3 - h4)) * v4
			- (h0*h0 * h1*h1 * h2*h2 * h4*h4) / ((h0 - h3)*(h1 - h3)*(h2 - h3)*(h3 - h4)) * v3
			+ (h0*h0 * h1*h1 * h3*h3 * h4*h4) / ((h0 - h2)*(h1 - h2)*(h2 - h3)*(h2 - h4)) * v2
			- (h0*h0 * h2*h2 * h3*h3 * h4*h4) / ((h0 - h1)*(h1 - h2)*(h1 - h3)*(h1 - h4)) * v1
			+ (h1*h1 * h2*h2 * h3*h3 * h4*h4) / ((h0 - h1)*(h0 - h2)*(h0 - h3)*(h0 - h4)) * v0
		) / (h0*h1*h2*h3 + h0*h1*h2*h4 + h0*h1*h3*h4 + h0*h2*h3*h4 + h1*h2*h3*h4);

		// Constant timestep case:
		/*
		return (
			  300. * v4
			- 300. * v3
			+ 200. * v2
			- 75.  * v1
			+ 12.  * v0
		) / 137.;
		*/
	}

	inline Matrix A6() {
		t6 = this->nextTime();
		t5 = this->times()[ 0];
		t4 = this->times()[-1];
		t3 = this->times()[-2];
		t2 = this->times()[-3];
		t1 = this->times()[-4];
		t0 = this->times()[-5];

		h5 = difference(t6, t5);
		h4 = difference(t6, t4);
		h3 = difference(t6, t3);
		h2 = difference(t6, t2);
		h1 = difference(t6, t1);
		h0 = difference(t6, t0);

		return
			domain->identity()
			+ (h0*h1*h2*h3*h4*h5) / (h0*h1*h2*h3*h4 + h0*h1*h2*h3*h5 + h0*h1*h2*h4*h5 + h0*h1*h3*h4*h5 + h0*h2*h3*h4*h5 + h1*h2*h3*h4*h5)
					* op->discretize(t6)
		;

		// Constant timestep case:
		/*
		return
			domain->identity()
			+ 60. / 147. * op->discretize(t) * dt()
		;
		*/
	}

	inline Vector b6() {
		const Vector
			&v5 = this->iterands()[ 0],
			&v4 = this->iterands()[-1],
			&v3 = this->iterands()[-2],
			&v2 = this->iterands()[-3],
			&v1 = this->iterands()[-4],
			&v0 = this->iterands()[-5]
		;

		return (
			  (h0*h0 * h1*h1 * h2*h2 * h3*h3 * h4*h4) / ((h0 - h5)*(h1 - h5)*(h2 - h5)*(h3 - h5)*(h4 - h5)) * v5
			- (h0*h0 * h1*h1 * h2*h2 * h3*h3 * h5*h5) / ((h0 - h4)*(h1 - h4)*(h2 - h4)*(h3 - h4)*(h4 - h5)) * v4
			+ (h0*h0 * h1*h1 * h2*h2 * h4*h4 * h5*h5) / ((h0 - h3)*(h1 - h3)*(h2 - h3)*(h3 - h4)*(h3 - h5)) * v3
			- (h0*h0 * h1*h1 * h3*h3 * h4*h4 * h5*h5) / ((h0 - h2)*(h1 - h2)*(h2 - h3)*(h2 - h4)*(h2 - h5)) * v2
			+ (h0*h0 * h2*h2 * h3*h3 * h4*h4 * h5*h5) / ((h0 - h1)*(h1 - h2)*(h1 - h3)*(h1 - h4)*(h1 - h5)) * v1
			- (h1*h1 * h2*h2 * h3*h3 * h4*h4 * h5*h5) / ((h0 - h1)*(h0 - h2)*(h0 - h3)*(h0 - h4)*(h0 - h5)) * v0
		) / (h0*h1*h2*h3*h4 + h0*h1*h2*h3*h5 + h0*h1*h2*h4*h5 + h0*h1*h3*h4*h5 + h0*h2*h3*h4*h5 + h1*h2*h3*h4*h5);

		// Constant timestep case:
		/*
		return (
			  360. * v5
			- 450. * v4
			+ 400. * v3
			- 225. * v2
			+ 72.  * v1
			- 10.  * v0
		) / 147.;
		*/
	}

public:

	template <typename D, typename L>
	LinearBDFBase(D &domain, L &op) noexcept : domain(&domain), op(&op) {
	}

};

////////////////////////////////////////////////////////////////////////////////

// TODO: This should be in a different file

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

protected:

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

public:

	template <typename D, typename L>
	CrankNicolson(D &domain, L &op) noexcept : domain(&domain), op(&op) {
	}

};

typedef CrankNicolson<false> ReverseCrankNicolson;
typedef CrankNicolson<true > ForwardCrankNicolson;

////////////////////////////////////////////////////////////////////////////////

template <bool Forward>
class LinearBDFOne : public LinearBDFBase<Forward> {

protected:

	virtual Matrix A() {
		return this->A1();
	}

	virtual Vector b() {
		return this->b1();
	}

public:

	template <typename D, typename L>
	LinearBDFOne(D &domain, L &op) noexcept
			: LinearBDFBase<Forward>(domain, op) {
	}

};

typedef LinearBDFOne<false> ReverseLinearBDFOne;
typedef LinearBDFOne<true > ForwardLinearBDFOne;

// More commonly referred to as the implicit Euler method
typedef ReverseLinearBDFOne ReverseImplicitEuler;
typedef ForwardLinearBDFOne ForwardImplicitEuler;

////////////////////////////////////////////////////////////////////////////////

template <bool Forward>
class LinearBDFTwo : public LinearBDFBase<Forward> {

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

	virtual void clear() {
		AA = &LinearBDFTwo::_A1;
		bb = &LinearBDFTwo::_b1;
	}

protected:

	virtual Matrix A() {
		return (this->*AA)();
	}

	virtual Vector b() {
		return (this->*bb)();
	}

public:

	template <typename D, typename L>
	LinearBDFTwo(D &domain, L &op) noexcept
			: LinearBDFBase<Forward>(domain, op) {
	}

};

typedef LinearBDFTwo<false> ReverseLinearBDFTwo;
typedef LinearBDFTwo<true > ForwardLinearBDFTwo;

////////////////////////////////////////////////////////////////////////////////

template <bool Forward>
class LinearBDFThree : public LinearBDFBase<Forward> {

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

	virtual void clear() {
		AA = &LinearBDFThree::_A1;
		bb = &LinearBDFThree::_b1;
	}

protected:

	virtual Matrix A() {
		return (this->*AA)();
	}

	virtual Vector b() {
		return (this->*bb)();
	}

public:

	template <typename D, typename L>
	LinearBDFThree(D &domain, L &op) noexcept
			: LinearBDFBase<Forward>(domain, op) {
	}

};

typedef LinearBDFThree<false> ReverseLinearBDFThree;
typedef LinearBDFThree<true > ForwardLinearBDFThree;

////////////////////////////////////////////////////////////////////////////////

template <bool Forward>
class LinearBDFFour : public LinearBDFBase<Forward> {

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
		AA = &LinearBDFFour::_A3;
		return this->A2();
	}

	Vector _b2() {
		bb = &LinearBDFFour::_b3;
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

	virtual void clear() {
		AA = &LinearBDFFour::_A1;
		bb = &LinearBDFFour::_b1;
	}

protected:

	virtual Matrix A() {
		return (this->*AA)();
	}

	virtual Vector b() {
		return (this->*bb)();
	}

public:

	template <typename D, typename L>
	LinearBDFFour(D &domain, L &op) noexcept
			: LinearBDFBase<Forward>(domain, op) {
	}

};

typedef LinearBDFFour<false> ReverseLinearBDFFour;
typedef LinearBDFFour<true > ForwardLinearBDFFour;

////////////////////////////////////////////////////////////////////////////////

template <bool Forward>
class LinearBDFFive : public LinearBDFBase<Forward> {

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
		AA = &LinearBDFFive::_A3;
		return this->A2();
	}

	Vector _b2() {
		bb = &LinearBDFFive::_b3;
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

	virtual void clear() {
		AA = &LinearBDFFive::_A1;
		bb = &LinearBDFFive::_b1;
	}

protected:

	virtual Matrix A() {
		return (this->*AA)();
	}

	virtual Vector b() {
		return (this->*bb)();
	}

public:

	template <typename D, typename L>
	LinearBDFFive(D &domain, L &op) noexcept
			: LinearBDFBase<Forward>(domain, op) {
	}

};

typedef LinearBDFFive<false> ReverseLinearBDFFive;
typedef LinearBDFFive<true > ForwardLinearBDFFive;

////////////////////////////////////////////////////////////////////////////////

template <bool Forward>
class LinearBDFSix : public LinearBDFBase<Forward> {

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
		AA = &LinearBDFSix::_A3;
		return this->A2();
	}

	Vector _b2() {
		bb = &LinearBDFSix::_b3;
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
		return this->b5();
	}

	virtual void clear() {
		AA = &LinearBDFSix::_A1;
		bb = &LinearBDFSix::_b1;
	}

protected:

	virtual Matrix A() {
		return (this->*AA)();
	}

	virtual Vector b() {
		return (this->*bb)();
	}

public:

	template <typename D, typename L>
	LinearBDFSix(D &domain, L &op) noexcept
			: LinearBDFBase<Forward>(domain, op) {
	}

};

typedef LinearBDFSix<false> ReverseLinearBDFSix;
typedef LinearBDFSix<true > ForwardLinearBDFSix;

} // QuantPDE

#endif

