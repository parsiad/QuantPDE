#ifndef QUANT_PDE_CORE_LINEAR_BDF
#define QUANT_PDE_CORE_LINEAR_BDF

namespace QuantPDE {

template <bool Forward, size_t Lookback>
class LinearBDFBase : public IterationNode {

	const DomainBase *domain;

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

#define QUANT_PDE_TMP \
	const Real t0 = time(0); \
	const Real h0 = difference(t1, t0);

	inline Matrix _A1(Real t1) {
		QUANT_PDE_TMP;

		return
			domain->identity()
			+ op->A(t1) * h0;
	}

	inline Vector _b1(Real t1) {
		QUANT_PDE_TMP;

		const Vector &v0 = iterand(0);

		return v0 + h0 * op->b(t1);
	}

#undef QUANT_PDE_TMP

/*
>>  A = [1 -1 -1;0 h1 h0;0 h1^2 h0^2]

A =

[ 1,   -1,   -1]
[ 0,   h1,   h0]
[ 0, h1^2, h0^2]

>> simplify(A\[0;1;0])

ans =

  (h0 + h1)/(h0*h1)
  h0/(h1*(h0 - h1))
 -h1/(h0*(h0 - h1))
*/

#define QUANT_PDE_TMP \
	const Real t1 = time(0); \
	const Real t0 = time(1); \
	const Real h1 = difference(t2, t1); \
	const Real h0 = difference(t2, t0); \
	const Real hh = (h0*h1)/(h0 + h1); \

	inline Matrix _A2(Real t2) {
		QUANT_PDE_TMP;

		return
			domain->identity()
			+ hh * op->A(t2)
		;

		// Constant timestep case:
		/*
		return
			domain->identity()
			- 2. / 3. * op->A(t2) * dt()
		;
		*/
	}

	inline Vector _b2(Real t2) {
		QUANT_PDE_TMP;

		const Vector
			&v1 = iterand(0),
			&v0 = iterand(1)
		;

		return
			(
				(
					  h0 / h1 * v1
					- h1 / h0 * v0
				) / ( h0 - h1 )
				+ op->b(t2)
			) * hh
		;

		// Constant timestep case:
		/*
		return (
			  4. * v1
			- 1. * v0
		) / 3.;
		*/
	}

#undef QUANT_PDE_TMP

/*
>> A = [1 -1 -1 -1;0 h2 h1 h0;0 h2^2 h1^2 h0^2;0 h2^3 h1^3 h0^3]

A =

[ 1,   -1,   -1,   -1]
[ 0,   h2,   h1,   h0]
[ 0, h2^2, h1^2, h0^2]
[ 0, h2^3, h1^3, h0^3]

>> simplify(A\[0;1;0;0])

ans =

 (h0*h1 + h0*h2 + h1*h2)/(h0*h1*h2)
   (h0*h1)/(h2*(h0 - h2)*(h1 - h2))
  -(h0*h2)/(h1*(h0 - h1)*(h1 - h2))
   (h1*h2)/(h0*(h0 - h1)*(h0 - h2))
*/

#define QUANT_PDE_TMP \
	const Real t2 = time(0); \
	const Real t1 = time(1); \
	const Real t0 = time(2); \
	const Real h2 = difference(t3, t2); \
	const Real h1 = difference(t3, t1); \
	const Real h0 = difference(t3, t0); \
	const Real hh = (h0*h1*h2)/(h0*h1 + h0*h2 + h1*h2);

	inline Matrix _A3(Real t3) {
		QUANT_PDE_TMP;

		return
			domain->identity()
			+ hh * op->A(t3)
		;

		// Constant timestep case:
		/*
		return
			domain->identity()
			- 6. / 11. * op->A(t3) * dt()
		;
		*/
	}

	inline Vector _b3(Real t3) {
		QUANT_PDE_TMP;

		const Vector
			&v2 = iterand(0),
			&v1 = iterand(1),
			&v0 = iterand(2)
		;

		return
			(
				  (h0*h1)/(h2*(h0 - h2)*(h1 - h2)) * v2
				- (h0*h2)/(h1*(h0 - h1)*(h1 - h2)) * v1
				+ (h1*h2)/(h0*(h0 - h1)*(h0 - h2)) * v0
				+ op->b(t3)
			) * hh
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

#undef QUANT_PDE_TMP

/*
>> A = [1 -1 -1 -1 -1;0 h3 h2 h1 h0;0 h3^2 h2^2 h1^2 h0^2;0 h3^3 h2^3 h1^3 h0^3;0 h3^4 h2^4 h1^4 h0^4]

A =

[ 1,   -1,   -1,   -1,   -1]
[ 0,   h3,   h2,   h1,   h0]
[ 0, h3^2, h2^2, h1^2, h0^2]
[ 0, h3^3, h2^3, h1^3, h0^3]
[ 0, h3^4, h2^4, h1^4, h0^4]

>> simplify(A\[0;1;0;0;0])

ans =

 (h0*h1*h2 + h0*h1*h3 + h0*h2*h3 + h1*h2*h3)/(h0*h1*h2*h3)
             (h0*h1*h2)/(h3*(h0 - h3)*(h1 - h3)*(h2 - h3))
            -(h0*h1*h3)/(h2*(h0 - h2)*(h1 - h2)*(h2 - h3))
             (h0*h2*h3)/(h1*(h0 - h1)*(h1 - h2)*(h1 - h3))
            -(h1*h2*h3)/(h0*(h0 - h1)*(h0 - h2)*(h0 - h3))
*/

#define QUANT_PDE_TMP \
	const Real t3 = time(0); \
	const Real t2 = time(1); \
	const Real t1 = time(2); \
	const Real t0 = time(3); \
	const Real h3 = difference(t4, t3); \
	const Real h2 = difference(t4, t2); \
	const Real h1 = difference(t4, t1); \
	const Real h0 = difference(t4, t0); \
	const Real hh = (h0*h1*h2*h3)/(h0*h1*h2 + h0*h1*h3 + h0*h2*h3 + h1*h2*h3);

	inline Matrix _A4(Real t4) {
		QUANT_PDE_TMP;

		return
			domain->identity()
			+ hh * op->A(t4)
		;


		// Constant timestep case:
		/*
		return
			domain->identity()
			- 12. / 25. * op->A(t4) * dt();
		*/
	}

	inline Vector _b4(Real t4) {
		QUANT_PDE_TMP;

		const Vector
			&v3 = iterand(0),
			&v2 = iterand(1),
			&v1 = iterand(2),
			&v0 = iterand(3)
		;

		return
			(
				  (h0*h1*h2)/(h3*(h0 - h3)*(h1 - h3)*(h2 - h3)) * v3
				- (h0*h1*h3)/(h2*(h0 - h2)*(h1 - h2)*(h2 - h3)) * v2
				+ (h0*h2*h3)/(h1*(h0 - h1)*(h1 - h2)*(h1 - h3)) * v1
				- (h1*h2*h3)/(h0*(h0 - h1)*(h0 - h2)*(h0 - h3)) * v0
				+ op->b(t4)
			) * hh
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

#undef QUANT_PDE_TMP

/*
>> A = [1 -1 -1 -1 -1 -1;0 h4 h3 h2 h1 h0;0 h4^2 h3^2 h2^2 h1^2 h0^2;0 h4^3 h3^3 h2^3 h1^3 h0^3;0 h4^4 h3^4 h2^4 h1^4 h0^4;0 h4^5 h3^5 h2^5 h1^5 h0^5]

A =

[ 1,   -1,   -1,   -1,   -1,   -1]
[ 0,   h4,   h3,   h2,   h1,   h0]
[ 0, h4^2, h3^2, h2^2, h1^2, h0^2]
[ 0, h4^3, h3^3, h2^3, h1^3, h0^3]
[ 0, h4^4, h3^4, h2^4, h1^4, h0^4]
[ 0, h4^5, h3^5, h2^5, h1^5, h0^5]

>> simplify(A\[0;1;0;0;0;0])

ans =

 (h0*h1*h2*h3 + h0*h1*h2*h4 + h0*h1*h3*h4 + h0*h2*h3*h4 + h1*h2*h3*h4)/(h0*h1*h2*h3*h4)
                             (h0*h1*h2*h3)/(h4*(h0 - h4)*(h1 - h4)*(h2 - h4)*(h3 - h4))
                            -(h0*h1*h2*h4)/(h3*(h0 - h3)*(h1 - h3)*(h2 - h3)*(h3 - h4))
                             (h0*h1*h3*h4)/(h2*(h0 - h2)*(h1 - h2)*(h2 - h3)*(h2 - h4))
                            -(h0*h2*h3*h4)/(h1*(h0 - h1)*(h1 - h2)*(h1 - h3)*(h1 - h4))
                             (h1*h2*h3*h4)/(h0*(h0 - h1)*(h0 - h2)*(h0 - h3)*(h0 - h4))
*/

#define QUANT_PDE_TMP \
	const Real t4 = time(0); \
	const Real t3 = time(1); \
	const Real t2 = time(2); \
	const Real t1 = time(3); \
	const Real t0 = time(4); \
	const Real h4 = difference(t5, t4); \
	const Real h3 = difference(t5, t3); \
	const Real h2 = difference(t5, t2); \
	const Real h1 = difference(t5, t1); \
	const Real h0 = difference(t5, t0); \
	const Real hh = (h0*h1*h2*h3*h4)/(h0*h1*h2*h3 + h0*h1*h2*h4 + h0*h1*h3*h4 + h0*h2*h3*h4 + h1*h2*h3*h4);

	inline Matrix _A5(Real t5) {
		QUANT_PDE_TMP;

		return
			domain->identity()
			+ hh * op->A(t5)
		;

		// Constant timestep case:
		/*
		return
			domain->identity()
			- 60. / 137. * op->A(t5) * dt()
		;
		*/
	}

	inline Vector _b5(Real t5) {
		QUANT_PDE_TMP;

		const Vector
			&v4 = iterand(0),
			&v3 = iterand(1),
			&v2 = iterand(2),
			&v1 = iterand(3),
			&v0 = iterand(4)
		;

		return
			(
				  (h0*h1*h2*h3)/(h4*(h0 - h4)*(h1 - h4)*(h2 - h4)*(h3 - h4)) * v4
				- (h0*h1*h2*h4)/(h3*(h0 - h3)*(h1 - h3)*(h2 - h3)*(h3 - h4)) * v3
				+ (h0*h1*h3*h4)/(h2*(h0 - h2)*(h1 - h2)*(h2 - h3)*(h2 - h4)) * v2
				- (h0*h2*h3*h4)/(h1*(h0 - h1)*(h1 - h2)*(h1 - h3)*(h1 - h4)) * v1
				+ (h1*h2*h3*h4)/(h0*(h0 - h1)*(h0 - h2)*(h0 - h3)*(h0 - h4)) * v0
				+ op->b(t5)
			) * hh
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

#undef QUANT_PDE_TMP

/*
>> A = [1 -1 -1 -1 -1 -1 -1;0 h5 h4 h3 h2 h1 h0;0 h5^2 h4^2 h3^2 h2^2 h1^2 h0^2;0 h5^3 h4^3 h3^3 h2^3 h1^3 h0^3;0 h5^4 h4^4 h3^4 h2^4 h1^4 h0^4;0 h5^5 h4^5 h3^5 h2^5 h1^5 h0^5;0 h5^6 h4^6 h3^6 h2^6 h1^6 h0^6]

A =

[ 1,   -1,   -1,   -1,   -1,   -1,   -1]
[ 0,   h5,   h4,   h3,   h2,   h1,   h0]
[ 0, h5^2, h4^2, h3^2, h2^2, h1^2, h0^2]
[ 0, h5^3, h4^3, h3^3, h2^3, h1^3, h0^3]
[ 0, h5^4, h4^4, h3^4, h2^4, h1^4, h0^4]
[ 0, h5^5, h4^5, h3^5, h2^5, h1^5, h0^5]
[ 0, h5^6, h4^6, h3^6, h2^6, h1^6, h0^6]

>> simplify(A\[0;1;0;0;0;0;0])

ans =

 (h0*h1*h2*h3*h4 + h0*h1*h2*h3*h5 + h0*h1*h2*h4*h5 + h0*h1*h3*h4*h5 + h0*h2*h3*h4*h5 + h1*h2*h3*h4*h5)/(h0*h1*h2*h3*h4*h5)
                                                   (h0*h1*h2*h3*h4)/(h5*(h0 - h5)*(h1 - h5)*(h2 - h5)*(h3 - h5)*(h4 - h5))
                                                  -(h0*h1*h2*h3*h5)/(h4*(h0 - h4)*(h1 - h4)*(h2 - h4)*(h3 - h4)*(h4 - h5))
                                                   (h0*h1*h2*h4*h5)/(h3*(h0 - h3)*(h1 - h3)*(h2 - h3)*(h3 - h4)*(h3 - h5))
                                                  -(h0*h1*h3*h4*h5)/(h2*(h0 - h2)*(h1 - h2)*(h2 - h3)*(h2 - h4)*(h2 - h5))
                                                   (h0*h2*h3*h4*h5)/(h1*(h0 - h1)*(h1 - h2)*(h1 - h3)*(h1 - h4)*(h1 - h5))
                                                  -(h1*h2*h3*h4*h5)/(h0*(h0 - h1)*(h0 - h2)*(h0 - h3)*(h0 - h4)*(h0 - h5))
*/

#define QUANT_PDE_TMP \
	const Real t5 = time(0); \
	const Real t4 = time(1); \
	const Real t3 = time(2); \
	const Real t2 = time(3); \
	const Real t1 = time(4); \
	const Real t0 = time(5); \
	const Real h5 = difference(t6, t5); \
	const Real h4 = difference(t6, t4); \
	const Real h3 = difference(t6, t3); \
	const Real h2 = difference(t6, t2); \
	const Real h1 = difference(t6, t1); \
	const Real h0 = difference(t6, t0); \
	const Real hh = (h0*h1*h2*h3*h4*h5)/(h0*h1*h2*h3*h4 + h0*h1*h2*h3*h5 + h0*h1*h2*h4*h5 + h0*h1*h3*h4*h5 + h0*h2*h3*h4*h5 + h1*h2*h3*h4*h5);

	inline Matrix _A6(Real t6) {
		QUANT_PDE_TMP;

		return
			domain->identity()
			+ hh * op->A(t6)
		;

		// Constant timestep case:
		/*
		return
			domain->identity()
			- 60. / 147. * op->A(t) * dt()
		;
		*/
	}

	inline Vector _b6(Real t6) {
		QUANT_PDE_TMP;

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
				  (h0*h1*h2*h3*h4)/(h5*(h0 - h5)*(h1 - h5)*(h2 - h5)*(h3 - h5)*(h4 - h5)) * v5
				- (h0*h1*h2*h3*h5)/(h4*(h0 - h4)*(h1 - h4)*(h2 - h4)*(h3 - h4)*(h4 - h5)) * v4
				+ (h0*h1*h2*h4*h5)/(h3*(h0 - h3)*(h1 - h3)*(h2 - h3)*(h3 - h4)*(h3 - h5)) * v3
				- (h0*h1*h3*h4*h5)/(h2*(h0 - h2)*(h1 - h2)*(h2 - h3)*(h2 - h4)*(h2 - h5)) * v2
				+ (h0*h2*h3*h4*h5)/(h1*(h0 - h1)*(h1 - h2)*(h1 - h3)*(h1 - h4)*(h1 - h5)) * v1
				- (h1*h2*h3*h4*h5)/(h0*(h0 - h1)*(h0 - h2)*(h0 - h3)*(h0 - h4)*(h0 - h5)) * v0
				+ op->b(t6)
			) * hh
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

#undef QUANT_PDE_TMP

public:

	template <typename D>
	LinearBDFBase(D &domain, LinearSystem &op) noexcept : domain(&domain),
			op(&op) {
	}

};

////////////////////////////////////////////////////////////////////////////////

template <bool Forward>
class LinearBDFOne : public LinearBDFBase<Forward, 1> {

	virtual Matrix A(Real t) {
		return this->_A1(t);
	}

	virtual Vector b(Real t) {
		return this->_b1(t);
	}

public:

	template <typename D>
	LinearBDFOne(D &domain, LinearSystem &op) noexcept
			: LinearBDFBase<Forward, 1>(domain, op) {
	}

	virtual bool isATheSame() const {
		return this->_isATheSame2();
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
	Matrix (LinearBDFTwo<Forward>::*_A)(Real);
	Vector (LinearBDFTwo<Forward>::*_b)(Real);
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

	virtual Matrix A(Real t) {
		return (this->*_A)(t);
	}

	virtual Vector b(Real t) {
		return (this->*_b)(t);
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
	Matrix (LinearBDFThree<Forward>::*_A)(Real);
	Vector (LinearBDFThree<Forward>::*_b)(Real);
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

	virtual Matrix A(Real t) {
		return (this->*_A)(t);
	}

	virtual Vector b(Real t) {
		return (this->*_b)(t);
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
	Matrix (LinearBDFFour<Forward>::*_A)(Real);
	Vector (LinearBDFFour<Forward>::*_b)(Real);
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

	virtual Matrix A(Real t) {
		return (this->*_A)(t);
	}

	virtual Vector b(Real t) {
		return (this->*_b)(t);
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
	Matrix (LinearBDFFive<Forward>::*_A)(Real);
	Vector (LinearBDFFive<Forward>::*_b)(Real);
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

	virtual Matrix A(Real t) {
		return (this->*_A)(t);
	}

	virtual Vector b(Real t) {
		return (this->*_b)(t);
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
	Matrix (LinearBDFSix<Forward>::*_A)(Real);
	Vector (LinearBDFSix<Forward>::*_b)(Real);
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

	virtual Matrix A(Real t) {
		return (this->*_A)(t);
	}

	virtual Vector b(Real t) {
		return (this->*_b)(t);
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

