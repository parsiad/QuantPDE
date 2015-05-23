#ifndef QUANT_PDE_CORE_CRANK_NICOLSON_HPP
#define QUANT_PDE_CORE_CRANK_NICOLSON_HPP

#include <limits> // std::numeric_limits

/**
 * QuantPDE namespace.
 */
namespace QuantPDE {

/**
 * The Crank-Nicolson method.
 *
 * Let \f$\Delta t\equiv t^1 - t^0\f$ where \f$t^1\f$ is the current time and
 * \f$t^0\f$ is the previous time. This creates the linear system
 * \f[ \left[I + A(t^1)\Delta t/2\right]\mathbf{x}^1 = \left[I - A(t^0)\Delta
 * t/2\right]\mathbf{x}^0 + \left[b\left(t^1\right) + b\left(t^0\right)\right]/2
 * \f]
 *
 * @tparam ThetaInverse 1 for implicit method, 2 for Crank-Nicolson, and
 *                      std::numeric_limits<Real>::infinity() for infinity.
 *
 * @see QuantPDE::TimeIteration
 */
template <bool Forward, int ThetaInverse = 2>
class CrankNicolson : public IterationNode {

	static_assert(ThetaInverse > 0, "ThetaInverse must be positive.");

	// GCC has trouble with this if NDEBUG is off
	#ifndef NDEBUG
	const Real theta = 1. / ((Real) ThetaInverse);
	#else
	static constexpr Real theta = 1. / ((Real) ThetaInverse);
	#endif

	const DomainBase &domain;
	LinearSystem &system;

	inline Real dt() const {
		const Real
			t1 = this->nextTime(),
			t0 = this->time(0)
		;

		const Real dt = Forward ? t1 - t0 : t0 - t1;
		assert(dt > epsilon);

		return dt;
	}

	virtual Matrix A(Real t1) {
		// Explicit
		if(theta < QuantPDE::epsilon) {
			return this->domain.identity();
		}

		// Not explicit
		return
			this->domain.identity()
			+ system.A(t1) * theta * dt()
		;
	}

	virtual Vector b(Real t1) {
		const Real    t0 = this->time(0);
		const Vector &v0 = this->iterand(0);

		// TODO: Optimize for explicit method

		return (
			this->domain.identity()
			- system.A(t0) * (1-theta) * dt()
		) * v0 + ( theta * system.b(t1) + (1-theta) * system.b(t0) );
	}

public:

	virtual bool isATheSame() const {
		// Explicit
		if(theta < QuantPDE::epsilon) {
			return true;
		}

		// Not explicit
		return this->isTimestepTheSame() && system.isATheSame();
	}

	/**
	 * @param domain The spatial domain to solve the problem on.
	 * @param system Creates
	 *               \f$A\left(t\right)\f$ and \f$b\left(t\right)\f$.
	 */
	template <typename D>
	CrankNicolson(
		D &domain,
		LinearSystem &system
	) noexcept :
		domain(domain),
		system(system)
	{
	}

};

typedef CrankNicolson<false> ReverseCrankNicolson;
typedef CrankNicolson<true>  ForwardCrankNicolson;

typedef CrankNicolson<false, std::numeric_limits<int>::max()>
		ReverseExplicitMethod;
typedef CrankNicolson<true,  std::numeric_limits<int>::max()>
		ForwardExplicitMethod;

} // QuantPDE

#endif
