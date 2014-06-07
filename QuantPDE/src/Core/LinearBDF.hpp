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
		return std::get<1>(this->iterands()[0]);
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

		return 4. / 3. * v_n1 - 1. / 3. * v_n0;
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

} // QuantPDE

#endif

