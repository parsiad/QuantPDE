#ifndef QUANT_PDE_CORE_LINEAR_BDF
#define QUANT_PDE_CORE_LINEAR_BDF

namespace QuantPDE {

template <size_t Dimension, size_t Controls, size_t Lookback, bool Forward>
class LinearBDFBase : public Linearizer<Lookback> {

	const Domain<Dimension> *domain;
	const LinearOperator<Dimension, Controls> *op;

	inline Real dt(Real implicitTime) {
		if(Forward) {
			return implicitTime - std::get<0>(this->iterands()[0]);
		} else {
			return std::get<0>(this->iterands()[0]) - implicitTime;
		}
	}

protected:

	inline Matrix A1(Real implicitTime) {
		return domain->identity() + op->discretize(implicitTime)
				* dt(implicitTime);
	}

	inline Vector b1(Real implicitTime) {
		return std::get<1>( this->iterands()[0] );
	}

	inline Matrix A2(Real implicitTime) {
		return domain->identity() + 2. / 3. * op->discretize(
				implicitTime) * dt(implicitTime);
	}

	inline Vector b2(Real implicitTime) {
		const Vector
			&v_n1 = std::get<1>( this->iterands()[ 0] ),
			&v_n0 = std::get<1>( this->iterands()[-1] )
		;

		return (
			  4. * v_n1
			- 1. * v_n0
		) / 3.;
	}

	inline Matrix A3(Real implicitTime) {
		return domain->identity() + 6. / 11. * op->discretize(
				implicitTime) * dt(implicitTime);
	}

	inline Vector b3(Real implicitTime) {
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

	inline Matrix A4(Real implicitTime) {
		return domain->identity() + 12. / 25. * op->discretize(
				implicitTime) * dt(implicitTime);
	}

	inline Vector b4(Real implicitTime) {
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

	inline Matrix A5(Real implicitTime) {
		return domain->identity() + 60. / 137. * op->discretize(
				implicitTime) * dt(implicitTime);
	}

	inline Vector b5(Real implicitTime) {
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

	inline Matrix A6(Real implicitTime) {
		return domain->identity() + 60. / 147. * op->discretize(
				implicitTime) * dt(implicitTime);
	}

	inline Vector b6(Real implicitTime) {
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

	template <typename I, typename D, typename L>
	LinearBDFBase(I &iteration, D &domain, L &op) noexcept
			: Linearizer<Lookback>(iteration), domain(&domain),
			op(&op) {
	}

};

template <size_t Dimension, size_t Controls = 0, bool Forward = false>
class LinearBDFOne : public LinearBDFBase<Dimension, Controls, 1, Forward> {

public:

	template <typename I, typename D, typename L>
	LinearBDFOne(I &iteration, D &domain, L &op) noexcept
			: LinearBDFBase<Dimension, Controls, 1, Forward>(
			iteration, domain, op) {
	}

	virtual Matrix A(Real implicitTime) {
		return this->A1(implicitTime);
	}

	virtual Vector b(Real implicitTime) {
		return this->b1(implicitTime);
	}

};

// Convenient alias
template <size_t Dimension, size_t Controls = 0, bool Forward = false>
using ImplicitMethod = LinearBDFOne<Dimension, Controls, Forward>;

template <size_t Controls = 0, bool Forward = false>
using LinearBDFOne1 = LinearBDFOne<1, Controls, Forward>;

template <size_t Controls = 0, bool Forward = false>
using LinearBDFOne2 = LinearBDFOne<2, Controls, Forward>;

template <size_t Controls = 0, bool Forward = false>
using LinearBDFOne3 = LinearBDFOne<3, Controls, Forward>;

////////////////////////////////////////////////////////////////////////////////

template <size_t Dimension, size_t Controls = 0, bool Forward = false>
class LinearBDFTwo : public LinearBDFBase<Dimension, Controls, 2, Forward> {

	const Domain<Dimension> *domain;
	const LinearOperator<Dimension, Controls> *op;

	Matrix (LinearBDFTwo<Dimension, Controls, Forward>::*AA)(Real);
	Vector (LinearBDFTwo<Dimension, Controls, Forward>::*bb)(Real);

	////////////////////////////////////////////////////////////////////////

	Matrix _A1(Real implicitTime) {
		AA = &LinearBDFTwo::A2;
		return this->A1(implicitTime);
	}

	Vector _b1(Real implicitTime) {
		bb = &LinearBDFTwo::b2;
		return this->b1(implicitTime);
	}

public:

	template <typename I, typename D, typename L>
	LinearBDFTwo(I &iteration, D &domain, L &op) noexcept
			: LinearBDFBase<Dimension, Controls, 2, Forward>(
			iteration, domain, op) {
	}

	virtual void clear() {
		AA = &LinearBDFTwo::_A1;
		bb = &LinearBDFTwo::_b1;
	}

	virtual Matrix A(Real implicitTime) {
		return (this->*AA)(implicitTime);
	}

	virtual Vector b(Real implicitTime) {
		return (this->*bb)(implicitTime);
	}

};

template <size_t Controls = 0, bool Forward = false>
using LinearBDFTwo1 = LinearBDFTwo<1, Controls, Forward>;

template <size_t Controls = 0, bool Forward = false>
using LinearBDFTwo2 = LinearBDFTwo<2, Controls, Forward>;

template <size_t Controls = 0, bool Forward = false>
using LinearBDFTwo3 = LinearBDFTwo<3, Controls, Forward>;

////////////////////////////////////////////////////////////////////////////////

template <size_t Dimension, size_t Controls = 0, bool Forward = false>
class LinearBDFThree : public LinearBDFBase<Dimension, Controls, 3, Forward> {

	const Domain<Dimension> *domain;
	const LinearOperator<Dimension, Controls> *op;

	Matrix (LinearBDFThree<Dimension, Controls, Forward>::*AA)(Real);
	Vector (LinearBDFThree<Dimension, Controls, Forward>::*bb)(Real);

	////////////////////////////////////////////////////////////////////////

	Matrix _A1(Real implicitTime) {
		AA = &LinearBDFThree::_A2;
		return this->A1(implicitTime);
	}

	Vector _b1(Real implicitTime) {
		bb = &LinearBDFThree::_b2;
		return this->b1(implicitTime);
	}

	Matrix _A2(Real implicitTime) {
		AA = &LinearBDFThree::A3;
		return this->A2(implicitTime);
	}

	Vector _b2(Real implicitTime) {
		bb = &LinearBDFThree::b3;
		return this->b2(implicitTime);
	}

public:

	template <typename I, typename D, typename L>
	LinearBDFThree(I &iteration, D &domain, L &op) noexcept
			: LinearBDFBase<Dimension, Controls, 3, Forward>(
			iteration, domain, op) {
	}

	virtual void clear() {
		AA = &LinearBDFThree::_A1;
		bb = &LinearBDFThree::_b1;
	}

	virtual Matrix A(Real implicitTime) {
		return (this->*AA)(implicitTime);
	}

	virtual Vector b(Real implicitTime) {
		return (this->*bb)(implicitTime);
	}

};

template <size_t Controls = 0, bool Forward = false>
using LinearBDFThree1 = LinearBDFThree<1, Controls, Forward>;

template <size_t Controls = 0, bool Forward = false>
using LinearBDFThree2 = LinearBDFThree<2, Controls, Forward>;

template <size_t Controls = 0, bool Forward = false>
using LinearBDFThree3 = LinearBDFThree<3, Controls, Forward>;

////////////////////////////////////////////////////////////////////////////////

template <size_t Dimension, size_t Controls = 0, bool Forward = false>
class LinearBDFFour : public LinearBDFBase<Dimension, Controls, 4, Forward> {

	const Domain<Dimension> *domain;
	const LinearOperator<Dimension, Controls> *op;

	Matrix (LinearBDFFour<Dimension, Controls, Forward>::*AA)(Real);
	Vector (LinearBDFFour<Dimension, Controls, Forward>::*bb)(Real);

	////////////////////////////////////////////////////////////////////////

	Matrix _A1(Real implicitTime) {
		AA = &LinearBDFFour::_A2;
		return this->A1(implicitTime);
	}

	Vector _b1(Real implicitTime) {
		bb = &LinearBDFFour::_b2;
		return this->b1(implicitTime);
	}

	Matrix _A2(Real implicitTime) {
		AA = &LinearBDFFour::A3;
		return this->A2(implicitTime);
	}

	Vector _b2(Real implicitTime) {
		bb = &LinearBDFFour::b3;
		return this->b2(implicitTime);
	}

	Matrix _A3(Real implicitTime) {
		AA = &LinearBDFFour::A4;
		return this->A3(implicitTime);
	}

	Vector _b3(Real implicitTime) {
		bb = &LinearBDFFour::b4;
		return this->b3(implicitTime);
	}

public:

	template <typename I, typename D, typename L>
	LinearBDFFour(I &iteration, D &domain, L &op) noexcept
			: LinearBDFBase<Dimension, Controls, 4, Forward>(
			iteration, domain, op) {
	}

	virtual void clear() {
		AA = &LinearBDFFour::_A1;
		bb = &LinearBDFFour::_b1;
	}

	virtual Matrix A(Real implicitTime) {
		return (this->*AA)(implicitTime);
	}

	virtual Vector b(Real implicitTime) {
		return (this->*bb)(implicitTime);
	}

};

template <size_t Controls = 0, bool Forward = false>
using LinearBDFFour1 = LinearBDFFour<1, Controls, Forward>;

template <size_t Controls = 0, bool Forward = false>
using LinearBDFFour2 = LinearBDFFour<2, Controls, Forward>;

template <size_t Controls = 0, bool Forward = false>
using LinearBDFFour3 = LinearBDFFour<3, Controls, Forward>;

////////////////////////////////////////////////////////////////////////////////

template <size_t Dimension, size_t Controls = 0, bool Forward = false>
class LinearBDFFive : public LinearBDFBase<Dimension, Controls, 5, Forward> {

	const Domain<Dimension> *domain;
	const LinearOperator<Dimension, Controls> *op;

	Matrix (LinearBDFFive<Dimension, Controls, Forward>::*AA)(Real);
	Vector (LinearBDFFive<Dimension, Controls, Forward>::*bb)(Real);

	////////////////////////////////////////////////////////////////////////

	Matrix _A1(Real implicitTime) {
		AA = &LinearBDFFive::_A2;
		return this->A1(implicitTime);
	}

	Vector _b1(Real implicitTime) {
		bb = &LinearBDFFive::_b2;
		return this->b1(implicitTime);
	}

	Matrix _A2(Real implicitTime) {
		AA = &LinearBDFFive::A3;
		return this->A2(implicitTime);
	}

	Vector _b2(Real implicitTime) {
		bb = &LinearBDFFive::b3;
		return this->b2(implicitTime);
	}

	Matrix _A3(Real implicitTime) {
		AA = &LinearBDFFive::_A4;
		return this->A3(implicitTime);
	}

	Vector _b3(Real implicitTime) {
		bb = &LinearBDFFive::_b4;
		return this->b3(implicitTime);
	}

	Matrix _A4(Real implicitTime) {
		AA = &LinearBDFFive::A5;
		return this->A4(implicitTime);
	}

	Vector _b4(Real implicitTime) {
		bb = &LinearBDFFive::b5;
		return this->b4(implicitTime);
	}

public:

	template <typename I, typename D, typename L>
	LinearBDFFive(I &iteration, D &domain, L &op) noexcept
			: LinearBDFBase<Dimension, Controls, 5, Forward>(
			iteration, domain, op) {
	}

	virtual void clear() {
		AA = &LinearBDFFive::_A1;
		bb = &LinearBDFFive::_b1;
	}

	virtual Matrix A(Real implicitTime) {
		return (this->*AA)(implicitTime);
	}

	virtual Vector b(Real implicitTime) {
		return (this->*bb)(implicitTime);
	}

};

template <size_t Controls = 0, bool Forward = false>
using LinearBDFFive1 = LinearBDFFive<1, Controls, Forward>;

template <size_t Controls = 0, bool Forward = false>
using LinearBDFFive2 = LinearBDFFive<2, Controls, Forward>;

template <size_t Controls = 0, bool Forward = false>
using LinearBDFFive3 = LinearBDFFive<3, Controls, Forward>;

////////////////////////////////////////////////////////////////////////////////

template <size_t Dimension, size_t Controls = 0, bool Forward = false>
class LinearBDFSix : public LinearBDFBase<Dimension, Controls, 6, Forward> {

	const Domain<Dimension> *domain;
	const LinearOperator<Dimension, Controls> *op;

	Matrix (LinearBDFSix<Dimension, Controls, Forward>::*AA)(Real);
	Vector (LinearBDFSix<Dimension, Controls, Forward>::*bb)(Real);

	////////////////////////////////////////////////////////////////////////

	Matrix _A1(Real implicitTime) {
		AA = &LinearBDFSix::_A2;
		return this->A1(implicitTime);
	}

	Vector _b1(Real implicitTime) {
		bb = &LinearBDFSix::_b2;
		return this->b1(implicitTime);
	}

	Matrix _A2(Real implicitTime) {
		AA = &LinearBDFSix::A3;
		return this->A2(implicitTime);
	}

	Vector _b2(Real implicitTime) {
		bb = &LinearBDFSix::b3;
		return this->b2(implicitTime);
	}

	Matrix _A3(Real implicitTime) {
		AA = &LinearBDFSix::_A4;
		return this->A3(implicitTime);
	}

	Vector _b3(Real implicitTime) {
		bb = &LinearBDFSix::_b4;
		return this->b3(implicitTime);
	}

	Matrix _A4(Real implicitTime) {
		AA = &LinearBDFSix::_A5;
		return this->A4(implicitTime);
	}

	Vector _b4(Real implicitTime) {
		bb = &LinearBDFSix::_b5;
		return this->b4(implicitTime);
	}

	Matrix _A5(Real implicitTime) {
		AA = &LinearBDFSix::A6;
		return this->A5(implicitTime);
	}

	Vector _b5(Real implicitTime) {
		bb = &LinearBDFSix::b6;
		return this->b6(implicitTime);
	}

public:

	template <typename I, typename D, typename L>
	LinearBDFSix(I &iteration, D &domain, L &op) noexcept
			: LinearBDFBase<Dimension, Controls, 6, Forward>(
			iteration, domain, op) {
	}

	virtual void clear() {
		AA = &LinearBDFSix::_A1;
		bb = &LinearBDFSix::_b1;
	}

	virtual Matrix A(Real implicitTime) {
		return (this->*AA)(implicitTime);
	}

	virtual Vector b(Real implicitTime) {
		return (this->*bb)(implicitTime);
	}

};

template <size_t Controls = 0, bool Forward = false>
using LinearBDFSix1 = LinearBDFSix<1, Controls, Forward>;

template <size_t Controls = 0, bool Forward = false>
using LinearBDFSix2 = LinearBDFSix<2, Controls, Forward>;

template <size_t Controls = 0, bool Forward = false>
using LinearBDFSix3 = LinearBDFSix<3, Controls, Forward>;


} // QuantPDE

#endif

