#ifndef QUANT_PDE_CORE_LINEAR_BDF
#define QUANT_PDE_CORE_LINEAR_BDF

// TODO: Cache A if constant

namespace QuantPDE {

template <bool Forward, size_t Lookback>
class LinearBDFBase : public LinearSystemIteration {

	const DomainBase *domain;

	Real t[Lookback + 1];
	Real h[Lookback];
	Real hh;

	inline Real difference(Real t1, Real t0) {
		return Forward ? t1 - t0 : t0 - t1;
	}

protected:

	LinearSystem *op;

	inline bool _isATheSame1() const {
		return false;
	}

	inline bool _isATheSame2() const {
		return isTimestepTheSame() && op->isATheSame();
	}

	inline Matrix _A1() {
		t[1] = nextTime();
		t[0] = time(0);

		h[0] = difference(t[1], t[0]);

		return
			domain->identity()
			+ op->A(t[1]) * h[0];
	}

	inline Vector _b1() {
		const Vector &v0 = iterand(0);
		return v0 + h[0] * op->b(t[1]);
	}

	inline Matrix _A2() {
		t[2] = nextTime();
		t[1] = time(0);
		t[0] = time(1);

		h[1] = difference(t[2], t[1]);
		h[0] = difference(t[2], t[0]);

		hh = (h[0] * h[1]) / (h[0] + h[1]);

		return
			domain->identity()
			+ hh * op->A(t[2])
		;

		// Constant timestep case:
		/*
		return
			domain->identity()
			- 2. / 3. * op->A(t[2]) * dt()
		;
		*/
	}

	inline Vector _b2() {
		const Vector
			&v1 = iterand(0),
			&v0 = iterand(1)
		;

		return
			(
				  h[0]*h[0] * v1
				- h[1]*h[1] * v0
			) / ( (h[0] + h[1]) * (h[0] - h[1]) )
			+ hh * op->b(t[2])
		;

		// Constant timestep case:
		/*
		return (
			  4. * v1
			- 1. * v0
		) / 3.;
		*/
	}

	inline Matrix _A3() {
		t[3] = nextTime();
		t[2] = time(0);
		t[1] = time(1);
		t[0] = time(2);

		h[2] = difference(t[3], t[2]);
		h[1] = difference(t[3], t[1]);
		h[0] = difference(t[3], t[0]);

		hh = ( h[0] * h[1] * h[2] ) / ( h[0] * h[1] + h[0] * h[2] + h[1] * h[2] );

		return
			domain->identity()
			+ hh * op->A(t[3])
		;

		// Constant timestep case:
		/*
		return
			domain->identity()
			- 6. / 11. * op->A(t[3]) * dt()
		;
		*/
	}

	inline Vector _b3() {
		const Vector
			&v2 = iterand(0),
			&v1 = iterand(1),
			&v0 = iterand(2)
		;

		return
			(
				  (h[0]*h[0] * h[1]*h[1]) / ((h[0] - h[2]) * (h[1] - h[2])) * v2
				- (h[0]*h[0] * h[2]*h[2]) / ((h[0] - h[1]) * (h[1] - h[2])) * v1
				+ (h[1]*h[1] * h[2]*h[2]) / ((h[0] - h[1]) * (h[0] - h[2])) * v0
			) / (h[0] * h[1] + h[0] * h[2] + h[1] * h[2])
			+ hh * op->b(t[3])
		;

		// Constant timestep case:
		/*
		return (
			  18. * v2
			- 9.  * v1
			+ 2.  * v0
		) / 11.;
		*/
	}

	inline Matrix _A4() {

		t[4] = nextTime();
		t[3] = time(0);
		t[2] = time(1);
		t[1] = time(2);
		t[0] = time(3);

		h[3] = difference(t[4], t[3]);
		h[2] = difference(t[4], t[2]);
		h[1] = difference(t[4], t[1]);
		h[0] = difference(t[4], t[0]);

		hh = (h[0]*h[1]*h[2]*h[3]) / (h[0]*h[1]*h[2] + h[0]*h[1]*h[3] + h[0]*h[2]*h[3] + h[1]*h[2]*h[3]);

		return
			domain->identity()
			+ hh * op->A(t[4])
		;


		// Constant timestep case:
		/*
		return
			domain->identity()
			- 12. / 25. * op->A(t[4]) * dt();
		*/
	}

	inline Vector _b4() {
		const Vector
			&v3 = iterand(0),
			&v2 = iterand(1),
			&v1 = iterand(2),
			&v0 = iterand(3)
		;

		return
			(
				  (h[0]*h[0] * h[1]*h[1] * h[2]*h[2]) / ((h[0] - h[3])*(h[1] - h[3])*(h[2] - h[3])) * v3
				- (h[0]*h[0] * h[1]*h[1] * h[3]*h[3]) / ((h[0] - h[2])*(h[1] - h[2])*(h[2] - h[3])) * v2
				+ (h[0]*h[0] * h[2]*h[2] * h[3]*h[3]) / ((h[0] - h[1])*(h[1] - h[2])*(h[1] - h[3])) * v1
				- (h[1]*h[1] * h[2]*h[2] * h[3]*h[3]) / ((h[0] - h[1])*(h[0] - h[2])*(h[0] - h[3])) * v0
			) / (h[0]*h[1]*h[2] + h[0]*h[1]*h[3] + h[0]*h[2]*h[3] + h[1]*h[2]*h[3])
			+ hh * op->b(t[4])
		;

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

	inline Matrix _A5() {
		t[5] = nextTime();
		t[4] = time(0);
		t[3] = time(1);
		t[2] = time(2);
		t[1] = time(3);
		t[0] = time(4);

		h[4] = difference(t[5], t[4]);
		h[3] = difference(t[5], t[3]);
		h[2] = difference(t[5], t[2]);
		h[1] = difference(t[5], t[1]);
		h[0] = difference(t[5], t[0]);

		hh = (h[0]*h[1]*h[2]*h[3]*h[4]) / (h[0]*h[1]*h[2]*h[3] + h[0]*h[1]*h[2]*h[4] + h[0]*h[1]*h[3]*h[4] + h[0]*h[2]*h[3]*h[4] + h[1]*h[2]*h[3]*h[4]);

		return
			domain->identity()
			+ hh * op->A(t[5])
		;

		// Constant timestep case:
		/*
		return
			domain->identity()
			- 60. / 137. * op->A(t[5]) * dt()
		;
		*/
	}

	inline Vector _b5() {
		const Vector
			&v4 = iterand(0),
			&v3 = iterand(1),
			&v2 = iterand(2),
			&v1 = iterand(3),
			&v0 = iterand(4)
		;

		return
			(
				  (h[0]*h[0] * h[1]*h[1] * h[2]*h[2] * h[3]*h[3]) / ((h[0] - h[4])*(h[1] - h[4])*(h[2] - h[4])*(h[3] - h[4])) * v4
				- (h[0]*h[0] * h[1]*h[1] * h[2]*h[2] * h[4]*h[4]) / ((h[0] - h[3])*(h[1] - h[3])*(h[2] - h[3])*(h[3] - h[4])) * v3
				+ (h[0]*h[0] * h[1]*h[1] * h[3]*h[3] * h[4]*h[4]) / ((h[0] - h[2])*(h[1] - h[2])*(h[2] - h[3])*(h[2] - h[4])) * v2
				- (h[0]*h[0] * h[2]*h[2] * h[3]*h[3] * h[4]*h[4]) / ((h[0] - h[1])*(h[1] - h[2])*(h[1] - h[3])*(h[1] - h[4])) * v1
				+ (h[1]*h[1] * h[2]*h[2] * h[3]*h[3] * h[4]*h[4]) / ((h[0] - h[1])*(h[0] - h[2])*(h[0] - h[3])*(h[0] - h[4])) * v0
			) / (h[0]*h[1]*h[2]*h[3] + h[0]*h[1]*h[2]*h[4] + h[0]*h[1]*h[3]*h[4] + h[0]*h[2]*h[3]*h[4] + h[1]*h[2]*h[3]*h[4])
			+ hh * op->b(t[5])
		;

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

	inline Matrix _A6() {
		t[6] = nextTime();
		t[5] = time(0);
		t[4] = time(1);
		t[3] = time(2);
		t[2] = time(3);
		t[1] = time(4);
		t[0] = time(5);

		h[5] = difference(t[6], t[5]);
		h[4] = difference(t[6], t[4]);
		h[3] = difference(t[6], t[3]);
		h[2] = difference(t[6], t[2]);
		h[1] = difference(t[6], t[1]);
		h[0] = difference(t[6], t[0]);

		hh = (h[0]*h[1]*h[2]*h[3]*h[4]*h[5]) / (h[0]*h[1]*h[2]*h[3]*h[4] + h[0]*h[1]*h[2]*h[3]*h[5] + h[0]*h[1]*h[2]*h[4]*h[5] + h[0]*h[1]*h[3]*h[4]*h[5] + h[0]*h[2]*h[3]*h[4]*h[5] + h[1]*h[2]*h[3]*h[4]*h[5]);

		return
			domain->identity()
			+ hh * op->A(t[6])
		;

		// Constant timestep case:
		/*
		return
			domain->identity()
			- 60. / 147. * op->A(t) * dt()
		;
		*/
	}

	inline Vector _b6() {
		const Vector
			&v5 = iterand(0),
			&v4 = iterand(1),
			&v3 = iterand(2),
			&v2 = iterand(3),
			&v1 = iterand(4),
			&v0 = iterand(5)
		;

		return
			(
				  (h[0]*h[0] * h[1]*h[1] * h[2]*h[2] * h[3]*h[3] * h[4]*h[4]) / ((h[0] - h[5])*(h[1] - h[5])*(h[2] - h[5])*(h[3] - h[5])*(h[4] - h[5])) * v5
				- (h[0]*h[0] * h[1]*h[1] * h[2]*h[2] * h[3]*h[3] * h[5]*h[5]) / ((h[0] - h[4])*(h[1] - h[4])*(h[2] - h[4])*(h[3] - h[4])*(h[4] - h[5])) * v4
				+ (h[0]*h[0] * h[1]*h[1] * h[2]*h[2] * h[4]*h[4] * h[5]*h[5]) / ((h[0] - h[3])*(h[1] - h[3])*(h[2] - h[3])*(h[3] - h[4])*(h[3] - h[5])) * v3
				- (h[0]*h[0] * h[1]*h[1] * h[3]*h[3] * h[4]*h[4] * h[5]*h[5]) / ((h[0] - h[2])*(h[1] - h[2])*(h[2] - h[3])*(h[2] - h[4])*(h[2] - h[5])) * v2
				+ (h[0]*h[0] * h[2]*h[2] * h[3]*h[3] * h[4]*h[4] * h[5]*h[5]) / ((h[0] - h[1])*(h[1] - h[2])*(h[1] - h[3])*(h[1] - h[4])*(h[1] - h[5])) * v1
				- (h[1]*h[1] * h[2]*h[2] * h[3]*h[3] * h[4]*h[4] * h[5]*h[5]) / ((h[0] - h[1])*(h[0] - h[2])*(h[0] - h[3])*(h[0] - h[4])*(h[0] - h[5])) * v0
			) / (h[0]*h[1]*h[2]*h[3]*h[4] + h[0]*h[1]*h[2]*h[3]*h[5] + h[0]*h[1]*h[2]*h[4]*h[5] + h[0]*h[1]*h[3]*h[4]*h[5] + h[0]*h[2]*h[3]*h[4]*h[5] + h[1]*h[2]*h[3]*h[4]*h[5])
			+ hh * op->b(t[6])
		;

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

	template <typename D>
	LinearBDFBase(D &domain, LinearSystem &op) noexcept : domain(&domain),
			op(&op) {
	}

};

////////////////////////////////////////////////////////////////////////////////

template <bool Forward>
class LinearBDFOne : public LinearBDFBase<Forward, 1> {

	virtual Matrix A() {
		return this->_A1();
	}

	virtual Vector b() {
		return this->_b1();
	}

public:

	template <typename D>
	LinearBDFOne(D &domain, LinearSystem &op) noexcept
			: LinearBDFBase<Forward, 1>(domain, op) {
	}

	virtual bool isATheSame() const {
		return this->_isATheSame2();
	}

	virtual int minimumLookback() const {
		return 1;
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

	bool   (LinearBDFTwo<Forward>::*_isATheSame)() const;
	Matrix (LinearBDFTwo<Forward>::*_A)();
	Vector (LinearBDFTwo<Forward>::*_b)();
	void   (LinearBDFTwo<Forward>::*_onIterationEnd)();

	////////////////////////////////////////////////////////////////////////

	void _onIterationEnd1() {
		_A = &LinearBDFTwo::_A2;
		_b = &LinearBDFTwo::_b2;
		_onIterationEnd = &LinearBDFTwo::_onIterationEnd2;
	}

	void _onIterationEnd2() {
		_isATheSame = &LinearBDFTwo::_isATheSame2;
		_onIterationEnd = &LinearBDFTwo::_onIterationEnd3;
	}

	void _onIterationEnd3() {
	}

	virtual void clear() {
		_isATheSame = &LinearBDFTwo::_isATheSame1;
		_A = &LinearBDFTwo::_A1;
		_b = &LinearBDFTwo::_b1;
		_onIterationEnd = &LinearBDFTwo::_onIterationEnd1;
	}

	virtual void onIterationEnd() {
		(this->*_onIterationEnd)();
	}

	virtual Matrix A() {
		return (this->*_A)();
	}

	virtual Vector b() {
		return (this->*_b)();
	}

	virtual int minimumLookback() const {
		return 2;
	}

public:

	template <typename D>
	LinearBDFTwo(D &domain, LinearSystem &op) noexcept
			: LinearBDFBase<Forward, 2>(domain, op) {
	}

	virtual bool isATheSame() const {
		return (this->*_isATheSame)();
	}

};

typedef LinearBDFTwo<false> ReverseLinearBDFTwo;
typedef LinearBDFTwo<true > ForwardLinearBDFTwo;

////////////////////////////////////////////////////////////////////////////////

template <bool Forward>
class LinearBDFThree : public LinearBDFBase<Forward, 3> {

	bool   (LinearBDFThree<Forward>::*_isATheSame)() const;
	Matrix (LinearBDFThree<Forward>::*_A)();
	Vector (LinearBDFThree<Forward>::*_b)();
	void   (LinearBDFThree<Forward>::*_onIterationEnd)();

	////////////////////////////////////////////////////////////////////////

	void _onIterationEnd1() {
		_A = &LinearBDFThree::_A2;
		_b = &LinearBDFThree::_b2;
		_onIterationEnd = &LinearBDFThree::_onIterationEnd2;
	}

	void _onIterationEnd2() {
		_A = &LinearBDFThree::_A3;
		_b = &LinearBDFThree::_b3;
		_onIterationEnd = &LinearBDFThree::_onIterationEnd3;
	}

	void _onIterationEnd3() {
		_isATheSame = &LinearBDFThree::_isATheSame2;
		_onIterationEnd = &LinearBDFThree::_onIterationEnd4;
	}

	void _onIterationEnd4() {
	}

	virtual void clear() {
		_isATheSame = &LinearBDFThree::_isATheSame1;
		_A = &LinearBDFThree::_A1;
		_b = &LinearBDFThree::_b1;
		_onIterationEnd = &LinearBDFThree::_onIterationEnd1;
	}

	virtual void onIterationEnd() {
		(this->*_onIterationEnd)();
	}

	virtual Matrix A() {
		return (this->*_A)();
	}

	virtual Vector b() {
		return (this->*_b)();
	}

	virtual int minimumLookback() const {
		return 3;
	}

public:

	template <typename D>
	LinearBDFThree(D &domain, LinearSystem &op) noexcept
			: LinearBDFBase<Forward, 3>(domain, op) {
	}

	virtual bool isATheSame() const {
		return (this->*_isATheSame)();
	}

};

typedef LinearBDFThree<false> ReverseLinearBDFThree;
typedef LinearBDFThree<true > ForwardLinearBDFThree;

////////////////////////////////////////////////////////////////////////////////

template <bool Forward>
class LinearBDFFour : public LinearBDFBase<Forward, 4> {

	bool   (LinearBDFFour<Forward>::*_isATheSame)() const;
	Matrix (LinearBDFFour<Forward>::*_A )();
	Vector (LinearBDFFour<Forward>::*_b )();
	void   (LinearBDFFour<Forward>::*_onIterationEnd)();

	////////////////////////////////////////////////////////////////////////

	void _onIterationEnd1() {
		_A = &LinearBDFFour::_A2;
		_b = &LinearBDFFour::_b2;
		_onIterationEnd = &LinearBDFFour::_onIterationEnd2;
	}

	void _onIterationEnd2() {
		_A = &LinearBDFFour::_A3;
		_b = &LinearBDFFour::_b3;
		_onIterationEnd = &LinearBDFFour::_onIterationEnd3;
	}

	void _onIterationEnd3() {
		_A = &LinearBDFFour::_A4;
		_b = &LinearBDFFour::_b4;
		_onIterationEnd = &LinearBDFFour::_onIterationEnd4;
	}

	void _onIterationEnd4() {
		_isATheSame = &LinearBDFFour::_isATheSame2;
		_onIterationEnd = &LinearBDFFour::_onIterationEnd5;
	}

	void _onIterationEnd5() {
	}

	virtual void clear() {
		_isATheSame = &LinearBDFFour::_isATheSame1;
		_A = &LinearBDFFour::_A1;
		_b = &LinearBDFFour::_b1;
		_onIterationEnd = &LinearBDFFour::_onIterationEnd1;
	}

	virtual void onIterationEnd() {
		(this->*_onIterationEnd)();
	}

	virtual Matrix A() {
		return (this->*_A)();
	}

	virtual Vector b() {
		return (this->*_b)();
	}

	virtual int minimumLookback() const {
		return 4;
	}

public:

	template <typename D>
	LinearBDFFour(D &domain, LinearSystem &op) noexcept
			: LinearBDFBase<Forward, 4>(domain, op) {
	}

	virtual bool isATheSame() const {
		return (this->*_isATheSame)();
	}

};

typedef LinearBDFFour<false> ReverseLinearBDFFour;
typedef LinearBDFFour<true > ForwardLinearBDFFour;

////////////////////////////////////////////////////////////////////////////////

template <bool Forward>
class LinearBDFFive : public LinearBDFBase<Forward, 5> {

	bool   (LinearBDFFive<Forward>::*_isATheSame)() const;
	Matrix (LinearBDFFive<Forward>::*_A )();
	Vector (LinearBDFFive<Forward>::*_b )();
	void   (LinearBDFFive<Forward>::*_onIterationEnd)();

	////////////////////////////////////////////////////////////////////////

	void _onIterationEnd1() {
		_A = &LinearBDFFive::_A2;
		_b = &LinearBDFFive::_b2;
		_onIterationEnd = &LinearBDFFive::_onIterationEnd2;
	}

	void _onIterationEnd2() {
		_A = &LinearBDFFive::_A3;
		_b = &LinearBDFFive::_b3;
		_onIterationEnd = &LinearBDFFive::_onIterationEnd3;
	}

	void _onIterationEnd3() {
		_A = &LinearBDFFive::_A4;
		_b = &LinearBDFFive::_b4;
		_onIterationEnd = &LinearBDFFive::_onIterationEnd4;
	}

	void _onIterationEnd4() {
		_A = &LinearBDFFive::_A5;
		_b = &LinearBDFFive::_b5;
		_onIterationEnd = &LinearBDFFive::_onIterationEnd5;
	}

	void _onIterationEnd5() {
		_isATheSame = &LinearBDFFive::_isATheSame2;
		_onIterationEnd = &LinearBDFFive::_onIterationEnd6;
	}

	void _onIterationEnd6() {
	}

	virtual void clear() {
		_isATheSame = &LinearBDFFive::_isATheSame1;
		_A = &LinearBDFFive::_A1;
		_b = &LinearBDFFive::_b1;
		_onIterationEnd = &LinearBDFFive::_onIterationEnd1;
	}

	virtual void onIterationEnd() {
		(this->*_onIterationEnd)();
	}

	virtual Matrix A() {
		return (this->*_A)();
	}

	virtual Vector b() {
		return (this->*_b)();
	}

	virtual int minimumLookback() const {
		return 5;
	}

public:

	template <typename D>
	LinearBDFFive(D &domain, LinearSystem &op) noexcept
			: LinearBDFBase<Forward, 5>(domain, op) {
	}

	virtual bool isATheSame() const {
		return (this->*_isATheSame)();
	}

};

typedef LinearBDFFive<false> ReverseLinearBDFFive;
typedef LinearBDFFive<true > ForwardLinearBDFFive;

////////////////////////////////////////////////////////////////////////////////

template <bool Forward>
class LinearBDFSix : public LinearBDFBase<Forward, 6> {

	bool   (LinearBDFSix<Forward>::*_isATheSame)() const;
	Matrix (LinearBDFSix<Forward>::*_A)();
	Vector (LinearBDFSix<Forward>::*_b)();
	void   (LinearBDFSix<Forward>::*_onIterationEnd)();

	////////////////////////////////////////////////////////////////////////

	void _onIterationEnd1() {
		_A = &LinearBDFSix::_A2;
		_b = &LinearBDFSix::_b2;
		_onIterationEnd = &LinearBDFSix::_onIterationEnd2;
	}

	void _onIterationEnd2() {
		_A = &LinearBDFSix::_A3;
		_b = &LinearBDFSix::_b3;
		_onIterationEnd = &LinearBDFSix::_onIterationEnd3;
	}

	void _onIterationEnd3() {
		_A = &LinearBDFSix::_A4;
		_b = &LinearBDFSix::_b4;
		_onIterationEnd = &LinearBDFSix::_onIterationEnd4;
	}

	void _onIterationEnd4() {
		_A = &LinearBDFSix::_A5;
		_b = &LinearBDFSix::_b5;
		_onIterationEnd = &LinearBDFSix::_onIterationEnd5;
	}

	void _onIterationEnd5() {
		_A = &LinearBDFSix::_A6;
		_b = &LinearBDFSix::_b6;
		_onIterationEnd = &LinearBDFSix::_onIterationEnd6;
	}

	void _onIterationEnd6() {
		_isATheSame = &LinearBDFSix::_isATheSame2;
		_onIterationEnd = &LinearBDFSix::_onIterationEnd7;
	}

	void _onIterationEnd7() {
	}

	virtual void clear() {
		_isATheSame = &LinearBDFSix::_isATheSame1;
		_A = &LinearBDFSix::_A1;
		_b = &LinearBDFSix::_b1;
		_onIterationEnd = &LinearBDFSix::_onIterationEnd1;
	}

	virtual void onIterationEnd() {
		(this->*_onIterationEnd)();
	}

	virtual Matrix A() {
		return (this->*_A)();
	}

	virtual Vector b() {
		return (this->*_b)();
	}

	virtual int minimumLookback() const {
		return 6;
	}

public:

	template <typename D>
	LinearBDFSix(D &domain, LinearSystem &op) noexcept
			: LinearBDFBase<Forward, 6>(domain, op) {
	}

	virtual bool isATheSame() const {
		return (this->*_isATheSame)();
	}

};

typedef LinearBDFSix<false> ReverseLinearBDFSix;
typedef LinearBDFSix<true > ForwardLinearBDFSix;

} // QuantPDE

#endif

