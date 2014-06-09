#ifndef QUANT_PDE_CORE_LINEAR_BDF
#define QUANT_PDE_CORE_LINEAR_BDF

namespace QuantPDE {

template <bool Forward, size_t Lookback>
class LinearBDFBase : public Linearizer {

	const DomainBase *domain;

	Real t[Lookback + 1];
	Real h[Lookback];

	inline Real difference(Real t1, Real t0) {
		return Forward ? t1 - t0 : t0 - t1;
	}

protected:

	const LinearOperator *op;

	inline Matrix A1() {
		t[1] = this->nextTime();
		t[0] = this->times()[0];

		h[0] = difference(t[1], t[0]);

		return
			domain->identity()
			+ op->discretize(t[1]) * h[0];
	}

	inline Vector b1() {
		const Vector &v0 = this->iterands()[0];
		return v0;
	}

	inline Matrix A2() {
		t[2] = this->nextTime();
		t[1] = this->times()[ 0];
		t[0] = this->times()[-1];

		h[1] = difference(t[2], t[1]);
		h[0] = difference(t[2], t[0]);

		return
			domain->identity()
			+ (h[0] * h[1]) / (h[0] + h[1]) * op->discretize(t[2])
		;

		// Constant timestep case:
		/*
		return
			domain->identity()
			+ 2. / 3. * op->discretize(t[2]) * dt()
		;
		*/
	}

	inline Vector b2() {
		const Vector
			&v1 = this->iterands()[ 0],
			&v0 = this->iterands()[-1]
		;

		return (
			  h[0]*h[0] * v1
			- h[1]*h[1] * v0
		) / ( (h[0] + h[1]) * (h[0] - h[1]) );

		// Constant timestep case:
		/*
		return (
			  4. * v1
			- 1. * v0
		) / 3.;
		*/
	}

	inline Matrix A3() {
		t[3] = this->nextTime();
		t[2] = this->times()[ 0];
		t[1] = this->times()[-1];
		t[0] = this->times()[-2];

		h[2] = difference(t[3], t[2]);
		h[1] = difference(t[3], t[1]);
		h[0] = difference(t[3], t[0]);

		return
			domain->identity()
			+ ( h[0] * h[1] * h[2] ) / ( h[0] * h[1] + h[0] * h[2] + h[1] * h[2] )
					* op->discretize(t[3])
		;

		// Constant timestep case:
		/*
		return
			domain->identity()
			+ 6. / 11. * op->discretize(t[3]) * dt()
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
			  (h[0]*h[0] * h[1]*h[1]) / ((h[0] - h[2]) * (h[1] - h[2])) * v2
			- (h[0]*h[0] * h[2]*h[2]) / ((h[0] - h[1]) * (h[1] - h[2])) * v1
			+ (h[1]*h[1] * h[2]*h[2]) / ((h[0] - h[1]) * (h[0] - h[2])) * v0
		) / (h[0] * h[1] + h[0] * h[2] + h[1] * h[2]);

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

		t[4] = this->nextTime();
		t[3] = this->times()[ 0];
		t[2] = this->times()[-1];
		t[1] = this->times()[-2];
		t[0] = this->times()[-3];

		h[3] = difference(t[4], t[3]);
		h[2] = difference(t[4], t[2]);
		h[1] = difference(t[4], t[1]);
		h[0] = difference(t[4], t[0]);

		return domain->identity()
			+ (h[0]*h[1]*h[2]*h[3]) / (h[0]*h[1]*h[2] + h[0]*h[1]*h[3] + h[0]*h[2]*h[3] + h[1]*h[2]*h[3])
					* op->discretize(t[4]);


		// Constant timestep case:
		/*
		return
			domain->identity()
			+ 12. / 25. * op->discretize(t[4]) * dt();
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
			  (h[0]*h[0] * h[1]*h[1] * h[2]*h[2]) / ((h[0] - h[3])*(h[1] - h[3])*(h[2] - h[3])) * v3
			- (h[0]*h[0] * h[1]*h[1] * h[3]*h[3]) / ((h[0] - h[2])*(h[1] - h[2])*(h[2] - h[3])) * v2
			+ (h[0]*h[0] * h[2]*h[2] * h[3]*h[3]) / ((h[0] - h[1])*(h[1] - h[2])*(h[1] - h[3])) * v1
			- (h[1]*h[1] * h[2]*h[2] * h[3]*h[3]) / ((h[0] - h[1])*(h[0] - h[2])*(h[0] - h[3])) * v0
		) / (h[0]*h[1]*h[2] + h[0]*h[1]*h[3] + h[0]*h[2]*h[3] + h[1]*h[2]*h[3]);

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
		t[5] = this->nextTime();
		t[4] = this->times()[ 0];
		t[3] = this->times()[-1];
		t[2] = this->times()[-2];
		t[1] = this->times()[-3];
		t[0] = this->times()[-4];

		h[4] = difference(t[5], t[4]);
		h[3] = difference(t[5], t[3]);
		h[2] = difference(t[5], t[2]);
		h[1] = difference(t[5], t[1]);
		h[0] = difference(t[5], t[0]);

		return
			domain->identity()
			+ (h[0]*h[1]*h[2]*h[3]*h[4]) / (h[0]*h[1]*h[2]*h[3] + h[0]*h[1]*h[2]*h[4] + h[0]*h[1]*h[3]*h[4] + h[0]*h[2]*h[3]*h[4] + h[1]*h[2]*h[3]*h[4])
					* op->discretize(t[5])
		;

		// Constant timestep case:
		/*
		return
			domain->identity()
			+ 60. / 137. * op->discretize(t[5]) * dt()
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
			  (h[0]*h[0] * h[1]*h[1] * h[2]*h[2] * h[3]*h[3]) / ((h[0] - h[4])*(h[1] - h[4])*(h[2] - h[4])*(h[3] - h[4])) * v4
			- (h[0]*h[0] * h[1]*h[1] * h[2]*h[2] * h[4]*h[4]) / ((h[0] - h[3])*(h[1] - h[3])*(h[2] - h[3])*(h[3] - h[4])) * v3
			+ (h[0]*h[0] * h[1]*h[1] * h[3]*h[3] * h[4]*h[4]) / ((h[0] - h[2])*(h[1] - h[2])*(h[2] - h[3])*(h[2] - h[4])) * v2
			- (h[0]*h[0] * h[2]*h[2] * h[3]*h[3] * h[4]*h[4]) / ((h[0] - h[1])*(h[1] - h[2])*(h[1] - h[3])*(h[1] - h[4])) * v1
			+ (h[1]*h[1] * h[2]*h[2] * h[3]*h[3] * h[4]*h[4]) / ((h[0] - h[1])*(h[0] - h[2])*(h[0] - h[3])*(h[0] - h[4])) * v0
		) / (h[0]*h[1]*h[2]*h[3] + h[0]*h[1]*h[2]*h[4] + h[0]*h[1]*h[3]*h[4] + h[0]*h[2]*h[3]*h[4] + h[1]*h[2]*h[3]*h[4]);

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
		t[6] = this->nextTime();
		t[5] = this->times()[ 0];
		t[4] = this->times()[-1];
		t[3] = this->times()[-2];
		t[2] = this->times()[-3];
		t[1] = this->times()[-4];
		t[0] = this->times()[-5];

		h[5] = difference(t[6], t[5]);
		h[4] = difference(t[6], t[4]);
		h[3] = difference(t[6], t[3]);
		h[2] = difference(t[6], t[2]);
		h[1] = difference(t[6], t[1]);
		h[0] = difference(t[6], t[0]);

		return
			domain->identity()
			+ (h[0]*h[1]*h[2]*h[3]*h[4]*h[5]) / (h[0]*h[1]*h[2]*h[3]*h[4] + h[0]*h[1]*h[2]*h[3]*h[5] + h[0]*h[1]*h[2]*h[4]*h[5] + h[0]*h[1]*h[3]*h[4]*h[5] + h[0]*h[2]*h[3]*h[4]*h[5] + h[1]*h[2]*h[3]*h[4]*h[5])
					* op->discretize(t[6])
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
			  (h[0]*h[0] * h[1]*h[1] * h[2]*h[2] * h[3]*h[3] * h[4]*h[4]) / ((h[0] - h[5])*(h[1] - h[5])*(h[2] - h[5])*(h[3] - h[5])*(h[4] - h[5])) * v5
			- (h[0]*h[0] * h[1]*h[1] * h[2]*h[2] * h[3]*h[3] * h[5]*h[5]) / ((h[0] - h[4])*(h[1] - h[4])*(h[2] - h[4])*(h[3] - h[4])*(h[4] - h[5])) * v4
			+ (h[0]*h[0] * h[1]*h[1] * h[2]*h[2] * h[4]*h[4] * h[5]*h[5]) / ((h[0] - h[3])*(h[1] - h[3])*(h[2] - h[3])*(h[3] - h[4])*(h[3] - h[5])) * v3
			- (h[0]*h[0] * h[1]*h[1] * h[3]*h[3] * h[4]*h[4] * h[5]*h[5]) / ((h[0] - h[2])*(h[1] - h[2])*(h[2] - h[3])*(h[2] - h[4])*(h[2] - h[5])) * v2
			+ (h[0]*h[0] * h[2]*h[2] * h[3]*h[3] * h[4]*h[4] * h[5]*h[5]) / ((h[0] - h[1])*(h[1] - h[2])*(h[1] - h[3])*(h[1] - h[4])*(h[1] - h[5])) * v1
			- (h[1]*h[1] * h[2]*h[2] * h[3]*h[3] * h[4]*h[4] * h[5]*h[5]) / ((h[0] - h[1])*(h[0] - h[2])*(h[0] - h[3])*(h[0] - h[4])*(h[0] - h[5])) * v0
		) / (h[0]*h[1]*h[2]*h[3]*h[4] + h[0]*h[1]*h[2]*h[3]*h[5] + h[0]*h[1]*h[2]*h[4]*h[5] + h[0]*h[1]*h[3]*h[4]*h[5] + h[0]*h[2]*h[3]*h[4]*h[5] + h[1]*h[2]*h[3]*h[4]*h[5]);

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

	virtual bool doesAChange() const {
		// TODO: Implement this correctly
		return true;
	}

};

////////////////////////////////////////////////////////////////////////////////

template <bool Forward>
class LinearBDFOne : public LinearBDFBase<Forward, 1> {

public:

	template <typename D, typename L>
	LinearBDFOne(D &domain, L &op) noexcept
			: LinearBDFBase<Forward, 1>(domain, op) {
	}

	virtual bool doesAChange() const {
		return this->op->isConstantInTime();
	}

	virtual Matrix A() {
		return this->A1();
	}

	virtual Vector b() {
		return this->b1();
	}

};

typedef LinearBDFOne<false> ReverseLinearBDFOne;
typedef LinearBDFOne<true > ForwardLinearBDFOne;

// More commonly referred to as the implicit Euler method
typedef ReverseLinearBDFOne ReverseImplicitEuler;
typedef ForwardLinearBDFOne ForwardImplicitEuler;

////////////////////////////////////////////////////////////////////////////////

template <bool Forward>
class LinearBDFTwo : public LinearBDFBase<Forward, 2> {

	Matrix (LinearBDFTwo<Forward>::*AA )();
	Vector (LinearBDFTwo<Forward>::*bb )();
	void   (LinearBDFTwo<Forward>::*end)();

	////////////////////////////////////////////////////////////////////////

	void _end1() {
		AA  = &LinearBDFTwo::A2;
		bb  = &LinearBDFTwo::b2;
		end = &LinearBDFTwo::_end2;
	}

	void _end2() {
	}

	virtual void clear() {
		AA  = &LinearBDFTwo::A1;
		bb  = &LinearBDFTwo::b1;
		end = &LinearBDFTwo::_end1;
	}

	virtual void onIterationEnd() {
		(this->*end)();
	}

public:

	template <typename D, typename L>
	LinearBDFTwo(D &domain, L &op) noexcept
			: LinearBDFBase<Forward, 2>(domain, op) {
	}

	virtual Matrix A() {
		return (this->*AA)();
	}

	virtual Vector b() {
		return (this->*bb)();
	}

};

typedef LinearBDFTwo<false> ReverseLinearBDFTwo;
typedef LinearBDFTwo<true > ForwardLinearBDFTwo;

////////////////////////////////////////////////////////////////////////////////

template <bool Forward>
class LinearBDFThree : public LinearBDFBase<Forward, 3> {

	Matrix (LinearBDFThree<Forward>::*AA )();
	Vector (LinearBDFThree<Forward>::*bb )();
	void   (LinearBDFThree<Forward>::*end)();

	////////////////////////////////////////////////////////////////////////

	void _end1() {
		AA  = &LinearBDFThree::A2;
		bb  = &LinearBDFThree::b2;
		end = &LinearBDFThree::_end2;
	}

	void _end2() {
		AA  = &LinearBDFThree::A3;
		bb  = &LinearBDFThree::b3;
		end = &LinearBDFThree::_end3;
	}

	void _end3() {
	}

	virtual void clear() {
		AA  = &LinearBDFThree::A1;
		bb  = &LinearBDFThree::b1;
		end = &LinearBDFThree::_end1;
	}

	virtual void onIterationEnd() {
		(this->*end)();
	}

public:

	template <typename D, typename L>
	LinearBDFThree(D &domain, L &op) noexcept
			: LinearBDFBase<Forward, 3>(domain, op) {
	}

	virtual Matrix A() {
		return (this->*AA)();
	}

	virtual Vector b() {
		return (this->*bb)();
	}

};

typedef LinearBDFThree<false> ReverseLinearBDFThree;
typedef LinearBDFThree<true > ForwardLinearBDFThree;

////////////////////////////////////////////////////////////////////////////////

template <bool Forward>
class LinearBDFFour : public LinearBDFBase<Forward, 4> {

	Matrix (LinearBDFFour<Forward>::*AA )();
	Vector (LinearBDFFour<Forward>::*bb )();
	void   (LinearBDFFour<Forward>::*end)();

	////////////////////////////////////////////////////////////////////////

	void _end1() {
		AA  = &LinearBDFFour::A2;
		bb  = &LinearBDFFour::b2;
		end = &LinearBDFFour::_end2;
	}

	void _end2() {
		AA  = &LinearBDFFour::A3;
		bb  = &LinearBDFFour::b3;
		end = &LinearBDFFour::_end3;
	}

	void _end3() {
		AA  = &LinearBDFFour::A4;
		bb  = &LinearBDFFour::b4;
		end = &LinearBDFFour::_end4;
	}

	void _end4() {
	}

	virtual void clear() {
		AA  = &LinearBDFFour::A1;
		bb  = &LinearBDFFour::b1;
		end = &LinearBDFFour::_end1;
	}

	virtual void onIterationEnd() {
		(this->*end)();
	}

public:

	template <typename D, typename L>
	LinearBDFFour(D &domain, L &op) noexcept
			: LinearBDFBase<Forward, 4>(domain, op) {
	}

	virtual Matrix A() {
		return (this->*AA)();
	}

	virtual Vector b() {
		return (this->*bb)();
	}

};

typedef LinearBDFFour<false> ReverseLinearBDFFour;
typedef LinearBDFFour<true > ForwardLinearBDFFour;

////////////////////////////////////////////////////////////////////////////////

template <bool Forward>
class LinearBDFFive : public LinearBDFBase<Forward, 5> {

	Matrix (LinearBDFFive<Forward>::*AA )();
	Vector (LinearBDFFive<Forward>::*bb )();
	void   (LinearBDFFive<Forward>::*end)();

	////////////////////////////////////////////////////////////////////////

	void _end1() {
		AA  = &LinearBDFFive::A2;
		bb  = &LinearBDFFive::b2;
		end = &LinearBDFFive::_end2;
	}

	void _end2() {
		AA  = &LinearBDFFive::A3;
		bb  = &LinearBDFFive::b3;
		end = &LinearBDFFive::_end3;
	}

	void _end3() {
		AA  = &LinearBDFFive::A4;
		bb  = &LinearBDFFive::b4;
		end = &LinearBDFFive::_end4;
	}

	void _end4() {
		AA  = &LinearBDFFive::A5;
		bb  = &LinearBDFFive::b5;
		end = &LinearBDFFive::_end5;
	}

	void _end5() {
	}

	virtual void clear() {
		AA  = &LinearBDFFive::A1;
		bb  = &LinearBDFFive::b1;
		end = &LinearBDFFive::_end1;
	}

	virtual void onIterationEnd() {
		(this->*end)();
	}

public:

	template <typename D, typename L>
	LinearBDFFive(D &domain, L &op) noexcept
			: LinearBDFBase<Forward, 5>(domain, op) {
	}

	virtual Matrix A() {
		return (this->*AA)();
	}

	virtual Vector b() {
		return (this->*bb)();
	}

};

typedef LinearBDFFive<false> ReverseLinearBDFFive;
typedef LinearBDFFive<true > ForwardLinearBDFFive;

////////////////////////////////////////////////////////////////////////////////

template <bool Forward>
class LinearBDFSix : public LinearBDFBase<Forward, 6> {

	Matrix (LinearBDFSix<Forward>::*AA )();
	Vector (LinearBDFSix<Forward>::*bb )();
	void   (LinearBDFSix<Forward>::*end)();

	////////////////////////////////////////////////////////////////////////

	void _end1() {
		AA  = &LinearBDFSix::A2;
		bb  = &LinearBDFSix::b2;
		end = &LinearBDFSix::_end2;
	}

	void _end2() {
		AA  = &LinearBDFSix::A3;
		bb  = &LinearBDFSix::b3;
		end = &LinearBDFSix::_end3;
	}

	void _end3() {
		AA  = &LinearBDFSix::A4;
		bb  = &LinearBDFSix::b4;
		end = &LinearBDFSix::_end4;
	}

	void _end4() {
		AA  = &LinearBDFSix::A5;
		bb  = &LinearBDFSix::b5;
		end = &LinearBDFSix::_end5;
	}

	void _end5() {
		AA  = &LinearBDFSix::A6;
		bb  = &LinearBDFSix::b6;
		end = &LinearBDFSix::_end6;
	}

	void _end6() {
	}

	virtual void clear() {
		AA  = &LinearBDFSix::A1;
		bb  = &LinearBDFSix::b1;
		end = &LinearBDFSix::_end1;
	}

	virtual void onIterationEnd() {
		(this->*end)();
	}

public:

	template <typename D, typename L>
	LinearBDFSix(D &domain, L &op) noexcept
			: LinearBDFBase<Forward, 6>(domain, op) {
	}

	virtual Matrix A() {
		return (this->*AA)();
	}

	virtual Vector b() {
		return (this->*bb)();
	}

};

typedef LinearBDFSix<false> ReverseLinearBDFSix;
typedef LinearBDFSix<true > ForwardLinearBDFSix;

} // QuantPDE

#endif

