#ifndef QUANT_PDE_CORE_CRANK_NICOLSON_HPP
#define QUANT_PDE_CORE_CRANK_NICOLSON_HPP

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
 * @see QuantPDE::TimeIteration
 */
template <bool Forward>
class CrankNicolson : public IterationNode {

	const DomainBase *domain;
	LinearSystem *system;

	inline Real dt() {
		const Real
			t1 = nextTime(),
			t0 = time(0)
		;

		return Forward ? t1 - t0 : t0 - t1;
	}

	virtual Matrix A(Real t1) {
		return
			domain->identity()
			+ system->A(t1) * dt() / 2.
		;
	}

	virtual Vector b(Real t1) {
		const Real    t0 = time(0);
		const Vector &v0 = iterand(0);

		return (
			domain->identity()
			- system->A(t0) * dt() / 2.
		) * v0 + ( system->b(t1) + system->b(t0) ) / 2.;
	}

public:

	virtual bool isATheSame() const {
		return isTimestepTheSame() && system->isATheSame();
	}

	/**
	 * @param domain The spatial domain to solve the problem on.
	 * @param system Creates
	 *               \f$A\left(t\right)\f$ and \f$b\left(t\right)\f$.
	 */
	template <typename D>
	CrankNicolson(D &domain, LinearSystem &system) noexcept :
			domain(&domain), system(&system) {
	}

};

typedef CrankNicolson<false> ReverseCrankNicolson;
typedef CrankNicolson<true > ForwardCrankNicolson;

} // QuantPDE

#endif

