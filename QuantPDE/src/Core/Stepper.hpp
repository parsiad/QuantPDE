#ifndef QUANT_PDE_CORE_STEPPER_HPP
#define QUANT_PDE_CORE_STEPPER_HPP

namespace QuantPDE {

/**
 * Steps through time at constant intervals.
 */
template <bool Forward>
class ConstantStepper final : public TimeIteration<Forward> {

	const Real dt;

public:

	virtual Real timestep() {
		return dt;
	}

	/**
	 * Constructor.
	 * @param startTime The initial time.
	 * @param endTime The expiry time (must be larger than the initial).
	 * @param dt The size of the timestep.
	 * @param lookback The number of iterands to keep track of.
	 */
	ConstantStepper(
		Real startTime,
		Real endTime,
		Real dt
	) noexcept :
		TimeIteration<Forward>(startTime, endTime),
		dt(dt)
	{
		assert(dt > 0);
	}

};

typedef ConstantStepper<false> ReverseConstantStepper;
typedef ConstantStepper<true > ForwardConstantStepper;

////////////////////////////////////////////////////////////////////////////////

/**
 * Steps through time at (self-adjusting) variable intervals.
 *
 * The initial interval size \f$\Delta t^0\f$ is given, and subsequent interval
 * sizes are computed by
 * \f$
 * \Delta t^n \equiv \frac{ \text{target} }{
 * 	\max_i \frac{ \left| v_i^{n+1} - v_i^n \right| }{
 * 		\max \left( \text{scale}, \max_j \left(
 * 			\left| v_j^{n+1} \right|
 * 		\right) \right)
 * 	} + \text{epsilon}
 * } \Delta t^{n+1}
 * \f$
 */
template <bool Forward>
class VariableStepper final : public TimeIteration<Forward> {

	Real dt, target, scale;
	Real (VariableStepper::*_step)();

	Real _step0() {
		_step = &VariableStepper::_step1;
		return dt;
	}

	Real _step1() {
		const Vector
			&v1 = this->iterand(0),
			&v0 = this->iterand(1)
		;

		dt *= target / (
			relativeError( v1, v0, scale )
			+ QuantPDE::epsilon
		);

		return dt;
	}

	virtual void clear() {
		_step = &VariableStepper::_step0;
	}

	virtual int minimumLookback() const {
		return 2;
	}

public:

	virtual Real timestep() {
		return (this->*_step)();
	}

	/**
	 * Constructor.
	 * @param startTime The initial time.
	 * @param endTime The expiry time (must be larger than the initial).
	 * @param dt The initial time step size.
	 * @param target The target relative error.
	 * @param scale The scale of the error (e.g. 1 for dollars).
	 * @param lookback The number of iterands to keep track of.
	 */
	VariableStepper(
		Real startTime,
		Real endTime,
		Real dt,
		Real target,
		Real scale = QuantPDE::scale
	) noexcept :
		TimeIteration<Forward>(startTime, endTime),
		dt(dt),
		target(target),
		scale(scale),
		_step(nullptr)
	{
		assert(dt > 0);
		assert(target > 0);
		assert(scale > 0);
	}

};

typedef VariableStepper<false> ReverseVariableStepper;
typedef VariableStepper<true > ForwardVariableStepper;

} // QuantPDE

#endif
