#ifndef QUANT_PDE_CORE_BDF_HPP
#define QUANT_PDE_CORE_BDF_HPP

namespace QuantPDE {

template <Index Dimension, bool Forward, size_t Lookback>
class BDFBase : public Discretization<Dimension> {

	inline Real difference(Real t1, Real t0) const {
		const Real dt = Forward ? t1 - t0 : t0 - t1;
		assert(dt > epsilon);
		return dt;
	}

protected:

	LinearSystem &op;

	inline bool _isATheSame1() const {
		return false;
	}

	inline bool _isATheSame2() const {
		return this->isTimestepTheSame() && op.isATheSame();
	}

#define QUANT_PDE_TMP \
	const Real t0 = this->time(0); \
	const Real h0 = difference(t1, t0);

	inline Matrix _A1(Real t1) {
		QUANT_PDE_TMP;

		return
			this->domain.identity()
			+ op.A(t1) * h0;
	}

	inline Vector _b1(Real t1) {
		QUANT_PDE_TMP;

		const Vector &v0 = this->iterand(0);

		return v0 + h0 * op.b(t1);
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
	const Real t1 = this->time(0); \
	const Real t0 = this->time(1); \
	const Real h1 = difference(t2, t1); \
	const Real h0 = difference(t2, t0); \
	const Real hh = (h0*h1)/(h0 + h1); \

	inline Matrix _A2(Real t2) {
		QUANT_PDE_TMP;

		return
			this->domain.identity()
			+ hh * op.A(t2)
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
			&v1 = this->iterand(0),
			&v0 = this->iterand(1)
		;

		return
			(
				(
					  h0 / h1 * v1
					- h1 / h0 * v0
				) / ( h0 - h1 )
				+ op.b(t2)
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
	const Real t2 = this->time(0); \
	const Real t1 = this->time(1); \
	const Real t0 = this->time(2); \
	const Real h2 = difference(t3, t2); \
	const Real h1 = difference(t3, t1); \
	const Real h0 = difference(t3, t0); \
	const Real hh = (h0*h1*h2)/(h0*h1 + h0*h2 + h1*h2);

	inline Matrix _A3(Real t3) {
		QUANT_PDE_TMP;

		return
			this->domain.identity()
			+ hh * op.A(t3)
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
			&v2 = this->iterand(0),
			&v1 = this->iterand(1),
			&v0 = this->iterand(2)
		;

		return
			(
				  (h0*h1)/(h2*(h0 - h2)*(h1 - h2)) * v2
				- (h0*h2)/(h1*(h0 - h1)*(h1 - h2)) * v1
				+ (h1*h2)/(h0*(h0 - h1)*(h0 - h2)) * v0
				+ op.b(t3)
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
	const Real t3 = this->time(0); \
	const Real t2 = this->time(1); \
	const Real t1 = this->time(2); \
	const Real t0 = this->time(3); \
	const Real h3 = difference(t4, t3); \
	const Real h2 = difference(t4, t2); \
	const Real h1 = difference(t4, t1); \
	const Real h0 = difference(t4, t0); \
	const Real hh = (h0*h1*h2*h3)/(h0*h1*h2 + h0*h1*h3 + h0*h2*h3 + h1*h2*h3);

	inline Matrix _A4(Real t4) {
		QUANT_PDE_TMP;

		return
			this->domain.identity()
			+ hh * op.A(t4)
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
			&v3 = this->iterand(0),
			&v2 = this->iterand(1),
			&v1 = this->iterand(2),
			&v0 = this->iterand(3)
		;

		return
			(
				  (h0*h1*h2)/(h3*(h0 - h3)*(h1 - h3)*(h2 - h3)) * v3
				- (h0*h1*h3)/(h2*(h0 - h2)*(h1 - h2)*(h2 - h3)) * v2
				+ (h0*h2*h3)/(h1*(h0 - h1)*(h1 - h2)*(h1 - h3)) * v1
				- (h1*h2*h3)/(h0*(h0 - h1)*(h0 - h2)*(h0 - h3)) * v0
				+ op.b(t4)
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
	const Real t4 = this->time(0); \
	const Real t3 = this->time(1); \
	const Real t2 = this->time(2); \
	const Real t1 = this->time(3); \
	const Real t0 = this->time(4); \
	const Real h4 = difference(t5, t4); \
	const Real h3 = difference(t5, t3); \
	const Real h2 = difference(t5, t2); \
	const Real h1 = difference(t5, t1); \
	const Real h0 = difference(t5, t0); \
	const Real hh = (h0*h1*h2*h3*h4)/(h0*h1*h2*h3 + h0*h1*h2*h4 + h0*h1*h3*h4 + h0*h2*h3*h4 + h1*h2*h3*h4);

	inline Matrix _A5(Real t5) {
		QUANT_PDE_TMP;

		return
			this->domain.identity()
			+ hh * op.A(t5)
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
			&v4 = this->iterand(0),
			&v3 = this->iterand(1),
			&v2 = this->iterand(2),
			&v1 = this->iterand(3),
			&v0 = this->iterand(4)
		;

		return
			(
				  (h0*h1*h2*h3)/(h4*(h0 - h4)*(h1 - h4)*(h2 - h4)*(h3 - h4)) * v4
				- (h0*h1*h2*h4)/(h3*(h0 - h3)*(h1 - h3)*(h2 - h3)*(h3 - h4)) * v3
				+ (h0*h1*h3*h4)/(h2*(h0 - h2)*(h1 - h2)*(h2 - h3)*(h2 - h4)) * v2
				- (h0*h2*h3*h4)/(h1*(h0 - h1)*(h1 - h2)*(h1 - h3)*(h1 - h4)) * v1
				+ (h1*h2*h3*h4)/(h0*(h0 - h1)*(h0 - h2)*(h0 - h3)*(h0 - h4)) * v0
				+ op.b(t5)
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
	const Real t5 = this->time(0); \
	const Real t4 = this->time(1); \
	const Real t3 = this->time(2); \
	const Real t2 = this->time(3); \
	const Real t1 = this->time(4); \
	const Real t0 = this->time(5); \
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
			this->domain.identity()
			+ hh * op.A(t6)
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
			&v5 = this->iterand(0),
			&v4 = this->iterand(1),
			&v3 = this->iterand(2),
			&v2 = this->iterand(3),
			&v1 = this->iterand(4),
			&v0 = this->iterand(5)
		;

		return
			(
				  (h0*h1*h2*h3*h4)/(h5*(h0 - h5)*(h1 - h5)*(h2 - h5)*(h3 - h5)*(h4 - h5)) * v5
				- (h0*h1*h2*h3*h5)/(h4*(h0 - h4)*(h1 - h4)*(h2 - h4)*(h3 - h4)*(h4 - h5)) * v4
				+ (h0*h1*h2*h4*h5)/(h3*(h0 - h3)*(h1 - h3)*(h2 - h3)*(h3 - h4)*(h3 - h5)) * v3
				- (h0*h1*h3*h4*h5)/(h2*(h0 - h2)*(h1 - h2)*(h2 - h3)*(h2 - h4)*(h2 - h5)) * v2
				+ (h0*h2*h3*h4*h5)/(h1*(h0 - h1)*(h1 - h2)*(h1 - h3)*(h1 - h4)*(h1 - h5)) * v1
				- (h1*h2*h3*h4*h5)/(h0*(h0 - h1)*(h0 - h2)*(h0 - h3)*(h0 - h4)*(h0 - h5)) * v0
				+ op.b(t6)
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
	BDFBase(
		D &domain,
		LinearSystem &op
	) noexcept :
		Discretization<Dimension>(domain),
		op(op)
	{
	}

};

////////////////////////////////////////////////////////////////////////////////

template <Index Dimension, bool Forward>
class BDFOne : public BDFBase<Dimension, Forward, 1> {

	virtual Matrix A(Real t) {
		return this->_A1(t);
	}

	virtual Vector bd(Real t) {
		return this->_b1(t);
	}

public:

	template <typename D>
	BDFOne(
		D &domain,
		LinearSystem &op
	) noexcept :
		BDFBase<Dimension, Forward, 1>(
			domain,
			op
		)
	{
	}

	virtual bool isATheSame() const {
		return this->_isATheSame2();
	}

};

template <Index Dimension>
using ReverseBDFOne = BDFOne<Dimension, false>;

template <Index Dimension>
using ForwardBDFOne = BDFOne<Dimension, true>;

typedef ReverseBDFOne<1> ReverseBDFOne1;
typedef ReverseBDFOne<2> ReverseBDFOne2;
typedef ReverseBDFOne<3> ReverseBDFOne3;

typedef ForwardBDFOne<1> ForwardBDFOne1;
typedef ForwardBDFOne<2> ForwardBDFOne2;
typedef ForwardBDFOne<3> ForwardBDFOne3;

// More commonly referred to as the implicit Euler method

typedef ReverseBDFOne1 ReverseImplicitEuler1;
typedef ReverseBDFOne2 ReverseImplicitEuler2;
typedef ReverseBDFOne3 ReverseImplicitEuler3;

typedef ForwardBDFOne1 ForwardImplicitEuler1;
typedef ForwardBDFOne2 ForwardImplicitEuler2;
typedef ForwardBDFOne3 ForwardImplicitEuler3;

////////////////////////////////////////////////////////////////////////////////

template <Index Dimension, bool Forward>
class BDFTwo : public BDFBase<Dimension, Forward, 2> {

	bool   (BDFTwo<Dimension, Forward>::*_isATheSame)() const;
	Matrix (BDFTwo<Dimension, Forward>::*_A)(Real);
	Vector (BDFTwo<Dimension, Forward>::*_b)(Real);
	void   (BDFTwo<Dimension, Forward>::*_onIterationEnd)();

	////////////////////////////////////////////////////////////////////////

	void _onIterationEnd1() {
		_A = &BDFTwo::_A2;
		_b = &BDFTwo::_b2;
		_onIterationEnd = &BDFTwo::_onIterationEnd2;
	}

	void _onIterationEnd2() {
		_isATheSame = &BDFTwo::_isATheSame2;
		_onIterationEnd = &BDFTwo::_onIterationEnd3;
	}

	void _onIterationEnd3() {
	}

	virtual void clear() {
		_isATheSame = &BDFTwo::_isATheSame1;
		_A = &BDFTwo::_A1;
		_b = &BDFTwo::_b1;
		_onIterationEnd = &BDFTwo::_onIterationEnd1;
	}

	virtual void onIterationEnd() {
		(this->*_onIterationEnd)();
	}

	virtual Matrix A(Real t) {
		return (this->*_A)(t);
	}

	virtual Vector bd(Real t) {
		return (this->*_b)(t);
	}

	virtual int minimumLookback() const {
		return 2;
	}

public:

	template <typename D>
	BDFTwo(
		D &domain,
		LinearSystem &op
	) noexcept :
		BDFBase<Dimension, Forward, 2>(
			domain,
			op
		),
		_isATheSame(nullptr),
		_A(nullptr),
		_b(nullptr),
		_onIterationEnd(nullptr)
	{
	}

	virtual bool isATheSame() const {
		return (this->*_isATheSame)();
	}

};

template <Index Dimension>
using ReverseBDFTwo = BDFTwo<Dimension, false>;

template <Index Dimension>
using ForwardBDFTwo = BDFTwo<Dimension, true>;

typedef ReverseBDFTwo<1> ReverseBDFTwo1;
typedef ReverseBDFTwo<2> ReverseBDFTwo2;
typedef ReverseBDFTwo<3> ReverseBDFTwo3;

typedef ForwardBDFTwo<1> ForwardBDFTwo1;
typedef ForwardBDFTwo<2> ForwardBDFTwo2;
typedef ForwardBDFTwo<3> ForwardBDFTwo3;

////////////////////////////////////////////////////////////////////////////////

template <Index Dimension, bool Forward>
class BDFThree : public BDFBase<Dimension, Forward, 3> {

	bool   (BDFThree<Dimension, Forward>::*_isATheSame)() const;
	Matrix (BDFThree<Dimension, Forward>::*_A)(Real);
	Vector (BDFThree<Dimension, Forward>::*_b)(Real);
	void   (BDFThree<Dimension, Forward>::*_onIterationEnd)();

	////////////////////////////////////////////////////////////////////////

	void _onIterationEnd1() {
		_A = &BDFThree::_A2;
		_b = &BDFThree::_b2;
		_onIterationEnd = &BDFThree::_onIterationEnd2;
	}

	void _onIterationEnd2() {
		_A = &BDFThree::_A3;
		_b = &BDFThree::_b3;
		_onIterationEnd = &BDFThree::_onIterationEnd3;
	}

	void _onIterationEnd3() {
		_isATheSame = &BDFThree::_isATheSame2;
		_onIterationEnd = &BDFThree::_onIterationEnd4;
	}

	void _onIterationEnd4() {
	}

	virtual void clear() {
		_isATheSame = &BDFThree::_isATheSame1;
		_A = &BDFThree::_A1;
		_b = &BDFThree::_b1;
		_onIterationEnd = &BDFThree::_onIterationEnd1;
	}

	virtual void onIterationEnd() {
		(this->*_onIterationEnd)();
	}

	virtual Matrix A(Real t) {
		return (this->*_A)(t);
	}

	virtual Vector bd(Real t) {
		return (this->*_b)(t);
	}

	virtual int minimumLookback() const {
		return 3;
	}

public:

	template <typename D>
	BDFThree(
		D &domain,
		LinearSystem &op
	) noexcept :
		BDFBase<Dimension, Forward, 3>(
			domain,
			op
		),
		_isATheSame(nullptr),
		_A(nullptr),
		_b(nullptr),
		_onIterationEnd(nullptr)
	{
	}

	virtual bool isATheSame() const {
		return (this->*_isATheSame)();
	}

};

template <Index Dimension>
using ReverseBDFThree = BDFThree<Dimension, false>;

template <Index Dimension>
using ForwardBDFThree = BDFThree<Dimension, true>;

typedef ReverseBDFThree<1> ReverseBDFThree1;
typedef ReverseBDFThree<2> ReverseBDFThree2;
typedef ReverseBDFThree<3> ReverseBDFThree3;

typedef ForwardBDFThree<1> ForwardBDFThree1;
typedef ForwardBDFThree<2> ForwardBDFThree2;
typedef ForwardBDFThree<3> ForwardBDFThree3;

////////////////////////////////////////////////////////////////////////////////

template <Index Dimension, bool Forward>
class BDFFour : public BDFBase<Dimension, Forward, 4> {

	bool   (BDFFour<Dimension, Forward>::*_isATheSame)() const;
	Matrix (BDFFour<Dimension, Forward>::*_A)(Real);
	Vector (BDFFour<Dimension, Forward>::*_b)(Real);
	void   (BDFFour<Dimension, Forward>::*_onIterationEnd)();

	////////////////////////////////////////////////////////////////////////

	void _onIterationEnd1() {
		_A = &BDFFour::_A2;
		_b = &BDFFour::_b2;
		_onIterationEnd = &BDFFour::_onIterationEnd2;
	}

	void _onIterationEnd2() {
		_A = &BDFFour::_A3;
		_b = &BDFFour::_b3;
		_onIterationEnd = &BDFFour::_onIterationEnd3;
	}

	void _onIterationEnd3() {
		_A = &BDFFour::_A4;
		_b = &BDFFour::_b4;
		_onIterationEnd = &BDFFour::_onIterationEnd4;
	}

	void _onIterationEnd4() {
		_isATheSame = &BDFFour::_isATheSame2;
		_onIterationEnd = &BDFFour::_onIterationEnd5;
	}

	void _onIterationEnd5() {
	}

	virtual void clear() {
		_isATheSame = &BDFFour::_isATheSame1;
		_A = &BDFFour::_A1;
		_b = &BDFFour::_b1;
		_onIterationEnd = &BDFFour::_onIterationEnd1;
	}

	virtual void onIterationEnd() {
		(this->*_onIterationEnd)();
	}

	virtual Matrix A(Real t) {
		return (this->*_A)(t);
	}

	virtual Vector bd(Real t) {
		return (this->*_b)(t);
	}

	virtual int minimumLookback() const {
		return 4;
	}

public:

	template <typename D>
	BDFFour(
		D &domain,
		LinearSystem &op
	) noexcept :
		BDFBase<Dimension, Forward, 4>(
			domain,
			op
		),
		_isATheSame(nullptr),
		_A(nullptr),
		_b(nullptr),
		_onIterationEnd(nullptr)
	{
	}

	virtual bool isATheSame() const {
		return (this->*_isATheSame)();
	}

};

template <Index Dimension>
using ReverseBDFFour = BDFFour<Dimension, false>;

template <Index Dimension>
using ForwardBDFFour = BDFFour<Dimension, true>;

typedef ReverseBDFFour<1> ReverseBDFFour1;
typedef ReverseBDFFour<2> ReverseBDFFour2;
typedef ReverseBDFFour<3> ReverseBDFFour3;

typedef ForwardBDFFour<1> ForwardBDFFour1;
typedef ForwardBDFFour<2> ForwardBDFFour2;
typedef ForwardBDFFour<3> ForwardBDFFour3;

////////////////////////////////////////////////////////////////////////////////

template <Index Dimension, bool Forward>
class BDFFive : public BDFBase<Dimension, Forward, 5> {

	bool   (BDFFive<Dimension, Forward>::*_isATheSame)() const;
	Matrix (BDFFive<Dimension, Forward>::*_A)(Real);
	Vector (BDFFive<Dimension, Forward>::*_b)(Real);
	void   (BDFFive<Dimension, Forward>::*_onIterationEnd)();

	////////////////////////////////////////////////////////////////////////

	void _onIterationEnd1() {
		_A = &BDFFive::_A2;
		_b = &BDFFive::_b2;
		_onIterationEnd = &BDFFive::_onIterationEnd2;
	}

	void _onIterationEnd2() {
		_A = &BDFFive::_A3;
		_b = &BDFFive::_b3;
		_onIterationEnd = &BDFFive::_onIterationEnd3;
	}

	void _onIterationEnd3() {
		_A = &BDFFive::_A4;
		_b = &BDFFive::_b4;
		_onIterationEnd = &BDFFive::_onIterationEnd4;
	}

	void _onIterationEnd4() {
		_A = &BDFFive::_A5;
		_b = &BDFFive::_b5;
		_onIterationEnd = &BDFFive::_onIterationEnd5;
	}

	void _onIterationEnd5() {
		_isATheSame = &BDFFive::_isATheSame2;
		_onIterationEnd = &BDFFive::_onIterationEnd6;
	}

	void _onIterationEnd6() {
	}

	virtual void clear() {
		_isATheSame = &BDFFive::_isATheSame1;
		_A = &BDFFive::_A1;
		_b = &BDFFive::_b1;
		_onIterationEnd = &BDFFive::_onIterationEnd1;
	}

	virtual void onIterationEnd() {
		(this->*_onIterationEnd)();
	}

	virtual Matrix A(Real t) {
		return (this->*_A)(t);
	}

	virtual Vector bd(Real t) {
		return (this->*_b)(t);
	}

	virtual int minimumLookback() const {
		return 5;
	}

public:

	template <typename D>
	BDFFive(
		D &domain,
		LinearSystem &op
	) noexcept :
		BDFBase<Dimension, Forward, 5>(
			domain,
			op
		),
		_isATheSame(nullptr),
		_A(nullptr),
		_b(nullptr),
		_onIterationEnd(nullptr)
	{
	}

	virtual bool isATheSame() const {
		return (this->*_isATheSame)();
	}

};

template <Index Dimension>
using ReverseBDFFive = BDFFive<Dimension, false>;

template <Index Dimension>
using ForwardBDFFive = BDFFive<Dimension, true>;

typedef ReverseBDFFive<1> ReverseBDFFive1;
typedef ReverseBDFFive<2> ReverseBDFFive2;
typedef ReverseBDFFive<3> ReverseBDFFive3;

typedef ForwardBDFFive<1> ForwardBDFFive1;
typedef ForwardBDFFive<2> ForwardBDFFive2;
typedef ForwardBDFFive<3> ForwardBDFFive3;

////////////////////////////////////////////////////////////////////////////////

template <Index Dimension, bool Forward>
class BDFSix : public BDFBase<Dimension, Forward, 6> {

	bool   (BDFSix<Dimension, Forward>::*_isATheSame)() const;
	Matrix (BDFSix<Dimension, Forward>::*_A)(Real);
	Vector (BDFSix<Dimension, Forward>::*_b)(Real);
	void   (BDFSix<Dimension, Forward>::*_onIterationEnd)();

	////////////////////////////////////////////////////////////////////////

	void _onIterationEnd1() {
		_A = &BDFSix::_A2;
		_b = &BDFSix::_b2;
		_onIterationEnd = &BDFSix::_onIterationEnd2;
	}

	void _onIterationEnd2() {
		_A = &BDFSix::_A3;
		_b = &BDFSix::_b3;
		_onIterationEnd = &BDFSix::_onIterationEnd3;
	}

	void _onIterationEnd3() {
		_A = &BDFSix::_A4;
		_b = &BDFSix::_b4;
		_onIterationEnd = &BDFSix::_onIterationEnd4;
	}

	void _onIterationEnd4() {
		_A = &BDFSix::_A5;
		_b = &BDFSix::_b5;
		_onIterationEnd = &BDFSix::_onIterationEnd5;
	}

	void _onIterationEnd5() {
		_A = &BDFSix::_A6;
		_b = &BDFSix::_b6;
		_onIterationEnd = &BDFSix::_onIterationEnd6;
	}

	void _onIterationEnd6() {
		_isATheSame = &BDFSix::_isATheSame2;
		_onIterationEnd = &BDFSix::_onIterationEnd7;
	}

	void _onIterationEnd7() {
	}

	virtual void clear() {
		_isATheSame = &BDFSix::_isATheSame1;
		_A = &BDFSix::_A1;
		_b = &BDFSix::_b1;
		_onIterationEnd = &BDFSix::_onIterationEnd1;
	}

	virtual void onIterationEnd() {
		(this->*_onIterationEnd)();
	}

	virtual Matrix A(Real t) {
		return (this->*_A)(t);
	}

	virtual Vector bd(Real t) {
		return (this->*_b)(t);
	}

	virtual int minimumLookback() const {
		return 6;
	}

public:

	template <typename D>
	BDFSix(
		D &domain,
		LinearSystem &op
	) noexcept :
		BDFBase<Dimension, Forward, 6>(
			domain,
			op
		),
		_isATheSame(nullptr),
		_A(nullptr),
		_b(nullptr),
		_onIterationEnd(nullptr)
	{
	}

	virtual bool isATheSame() const {
		return (this->*_isATheSame)();
	}

};

template <Index Dimension>
using ReverseBDFSix = BDFSix<Dimension, false>;

template <Index Dimension>
using ForwardBDFSix = BDFSix<Dimension, true>;

typedef ReverseBDFSix<1> ReverseBDFSix1;
typedef ReverseBDFSix<2> ReverseBDFSix2;
typedef ReverseBDFSix<3> ReverseBDFSix3;

typedef ForwardBDFSix<1> ForwardBDFSix1;
typedef ForwardBDFSix<2> ForwardBDFSix2;
typedef ForwardBDFSix<3> ForwardBDFSix3;

} // QuantPDE

#endif

