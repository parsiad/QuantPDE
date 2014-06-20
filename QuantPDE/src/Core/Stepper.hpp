#ifndef QUANT_PDE_CORE_STEPPER
#define QUANT_PDE_CORE_STEPPER

namespace QuantPDE {

/**
 * Steps through time at constant intervals.
 */
template <bool>
class ConstantStepper final : public Iteration {

	Real startTime, endTime, dt;
	unsigned n, steps;

	virtual void clear() {
		n = 0;
	}

	virtual Real initialTime(Real) const;

	virtual Real timestep();

	virtual bool done() const {
		return n >= steps;
	}

public:

	/**
	 * Constructor.
	 * @param startTime The initial time.
	 * @param endTime The expiry time (must be larger than the initial).
	 * @param steps The number of steps to take.
	 * @param lookback The number of iterands to keep track of.
	 */
	ConstantStepper(
		Real startTime,
		Real endTime,
		unsigned steps,
		int lookback = DEFAULT_LOOKBACK
	) noexcept :
		Iteration(lookback),
		startTime(startTime),
		endTime(endTime),
		dt( (endTime - startTime) / steps ),
		steps(steps)
	{
		assert(startTime >= 0.);
		assert(startTime < endTime);
		assert(steps > 0);
	}

};

typedef ConstantStepper<false> ReverseConstantStepper;
typedef ConstantStepper<true > ForwardConstantStepper;

template <>
Real ReverseConstantStepper::initialTime(Real) const {
	return endTime;
}

template <>
Real ForwardConstantStepper::initialTime(Real) const {
	return startTime;
}

template <>
Real ReverseConstantStepper::timestep() {
	n++;
	return -dt;
}

template <>
Real ForwardConstantStepper::timestep() {
	n++;
	return dt;
}

////////////////////////////////////////////////////////////////////////////////

/**
 * Steps through time at (self-adjusting) variable intervals.
 */
template <bool>
class VariableStepper final : public Iteration {

	Real startTime, endTime, dt, target, scale, time;
	Real (VariableStepper::*_step)();

	inline void __step1() {
		const Vector
			&v1 = iterand(0),
			&v0 = iterand(1)
		;

		// Tested 2014-06-08
		dt *= target / ( v1 - v0 ).cwiseAbs().cwiseQuotient(
			( scale * Vector::Ones( v1.size() ) ).cwiseMax(
				v1.cwiseAbs().cwiseMax(
					v0.cwiseAbs()
				)
			)
		).maxCoeff();
	}

	Real _step0();

	Real _step1();

	virtual void clear() {
		_step = &VariableStepper::_step0;
		time = initialTime(0.);
	}

	virtual Real initialTime(Real) const;

	virtual Real timestep() {
		return (this->*_step)();
	}

	virtual bool done() const;

public:

	/**
	 * Constructor.
	 * @param startTime The initial time.
	 * @param endTime The expiry time (must be larger than the initial).
	 * @param dt The initial time step size.
	 * @param target The target time step size.
	 * @param scale The scale of the error (e.g. 1 for dollars).
	 * @param lookback The number of iterands to keep track of.
	 */
	VariableStepper(
		Real startTime,
		Real endTime,
		Real dt,
		Real target,
		Real scale = 1,
		int lookback = DEFAULT_LOOKBACK
	) noexcept :
		Iteration(lookback),
		startTime(startTime),
		endTime(endTime),
		dt(dt),
		target(target),
		scale(scale)
	{
		assert(startTime >= 0.);
		assert(startTime < endTime);
		assert(dt > 0);
		assert(target > 0);
		assert(scale > 0);
		assert(lookback >= 2); // Need at least two iterands to
		                       // variable step
	}

};

typedef VariableStepper<false> ReverseVariableStepper;
typedef VariableStepper<true > ForwardVariableStepper;

template <>
Real ReverseVariableStepper::_step1() {
	__step1();

	if(time - dt < startTime) {
		dt = time - startTime;
		time = startTime;
	} else {
		time -= dt;
	}

	return -dt;
}

template <>
Real ForwardVariableStepper::_step1() {
	__step1();

	if(time + dt > endTime) {
		dt = endTime - time;
		time = endTime;
	} else {
		time += dt;
	}

	return dt;
}

template <>
Real ReverseVariableStepper::_step0() {
	_step = &VariableStepper::_step1;
	return -dt;
}

template <>
Real ForwardVariableStepper::_step0() {
	_step = &VariableStepper::_step1;
	return dt;
}

template <>
Real ReverseVariableStepper::initialTime(Real) const {
	return endTime;
}

template <>
Real ForwardVariableStepper::initialTime(Real) const {
	return startTime;
}

template <>
bool ReverseVariableStepper::done() const {
	return time <= startTime;
}

template <>
bool ForwardVariableStepper::done() const {
	return time >= endTime;
}

} // QuantPDE

#endif

