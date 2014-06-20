#ifndef QUANT_PDE_CORE_STEPPER
#define QUANT_PDE_CORE_STEPPER

#include <functional> // std::function
#include <memory>     // std::unique_ptr
#include <queue>      // std::priority_queue
#include <tuple>      // std::tuple
#include <vector>     // std::vector

namespace QuantPDE {

namespace Metafunctions {

namespace NaryFunctionSignatureHelpers {

template <typename R, typename ...Ts>
using TransformTarget = R (const Function<sizeof...(Ts)> &, Real, Ts...);

template <unsigned N>
using TransformSignature = Type<TransformTarget, Matrix, N, Vector>;

} // NaryFunctionSignatureHelpers

} // Metafunctions

template <unsigned N>
using Transform = std::function< Metafunctions::NaryFunctionSignatureHelpers
		::TransformSignature<N> >;

typedef Transform<1> Transform1;
typedef Transform<2> Transform2;
typedef Transform<3> Transform3;

/**
 * Transforms a vector.
 */
class EventBase {

	virtual Vector doEvent(const Vector &vector) const = 0;
	virtual Vector doEvent(Vector &&vector) const = 0;

public:

	template <typename V>
	Vector operator()(V &&vector) const {
		return doEvent( std::forward<V>(vector) );
	}

};

/**
 * An event describes a manner in which a solution on a domain is transformed.
 * In particular, an event creates an interpolant from the solution data,
 * transforms the interpolant using a lambda function, and maps the result back
 * onto a vector (whose components correspond to the domain nodes).
 *
 * For example, a Bermudan put event has the associated lambda function:
 * \code{.cpp}
 * [K] (const Function1 &V, Real S) {
 * 	// V is the solution prior to the event
 * 	// S is the coordinate (asset price)
 * 	max( V(S), K - S );
 * }
 * \endcode
 */
template <Index Dimension>
class Event : public EventBase {

	typedef std::unique_ptr<InterpolantFactory<Dimension>> In;
	typedef std::unique_ptr<Interpolant<Dimension>> I;
	typedef std::unique_ptr<Map<Dimension>> Out;

	const Domain<Dimension> *domain;

	Transform<Dimension> transform;

	In in;
	Out out;

	template <int ...Indices>
	static inline Real packAndCall(
		const Transform<Dimension> &transform,
		const Function<Dimension> &solution,
		const Real *array,
		Metafunctions::GenerateSequenceHelpers::Sequence<
				Indices...>
	) {
		return transform( solution, array[Indices]... );
	}

	template <int N>
	static inline Real packAndCall(
		const Transform<Dimension> &transform,
		const Function<Dimension> &solution,
		const Real *array
	) {
		return Event::packAndCall(
			transform,
			solution,
			array,
			GenerateSequence<N>()
		);
	}

	template <typename V>
	Vector _doEvent(V &&vector) const {
		I i = in->make(std::forward<V>(vector));

		Vector v = domain->vector();
		for(auto node : domain->accessor(v)) {
			*node = packAndCall<Dimension>(
				transform,
				in,
				(&node).data()
			);
		}

		return (*out)(v);
	}

	virtual Vector doEvent(const Vector &vector) const {
		return _doEvent(vector);
	}

	virtual Vector doEvent(Vector &&vector) const {
		return _doEvent(std::move(vector));
	}

public:

	// TODO: Make constructor for general interpolation/mapping

	/**
	 * Default constructor for rectilinear grids.
	 */
	template <typename T, typename G>
	Event(
		T &&transform,
		G &grid
	) noexcept :
		transform( std::forward<T>( transform ) ),
		in(In(new typename PiecewiseLinear<Dimension>::Factory(grid))),
		out( Out(new PointwiseMap<Dimension>(grid)) ) {
	}

};

////////////////////////////////////////////////////////////////////////////////

/**
 * Creates time steppers.
 */
template <bool>
class StepperFactory {

	typedef std::unique_ptr<Iteration> I;

public:

	/**
	 * Destructor.
	 */
	virtual ~StepperFactory() {
	}

	/**
	 * Creates a time stepper.
	 * @param startTime The initial time.
	 * @param endTime The end time.
	 * @return A time stepper.
	 */
	virtual I make(Real startTime, Real endTime) const = 0;

};

typedef StepperFactory<false> ReverseStepperFactory;
typedef StepperFactory<true > ForwardStepperFactory;

////////////////////////////////////////////////////////////////////////////////

/**
 * This pure virtual class should not be directly extended.
 * @see QuantPDE::EventIteration
 */
template <Index Dimension, bool Forward>
class EventIterationBase : public Iteration {

protected:

	typedef typename std::conditional< Forward, std::less<Real>,
			std::greater<Real> >::type Order;

	typedef std::tuple<int, Real, Event<Dimension>> T;

	struct TimeOrder {
		// Returns true if a goes before b in the ordering
		bool operator()(const T &a, const T &b) const {
			return
				Order()( std::get<1>(a), std::get<1>(b) )
				|| (
					std::get<1>(a) == std::get<1>(b)
					&& Order()( std::get<0>(a),
							std::get<0><(b) )
				)
			;
		}
	};

	// TODO: What happens if events are set to the same time?

	typedef std::unique_ptr<Iteration> I;

	////////////////////////////////////////////////////////////////////////

	int id;
	std::priority_queue<T, std::vector<T>, TimeOrder> events;

	Real startTime, endTime, time;

	const StepperFactory<Forward> *factory;
	I iteration;

private:

	virtual Real timestep() {
		return iteration->timestep();
	}

public:

	/**
	 * Constructor.
	 */
	template <typename F>
	EventIterationBase(Real startTime, Real endTime, F &factory) noexcept
			: id(0), startTime(startTime), endTime(endTime),
			factory(&factory) {
		assert(startTime >= 0);
		assert(startTime < endTime);
	}

	/**
	 * Adds an event to be processed.
	 * @param time The time at which the event occurs.
	 * @param event The event.
	 */
	template <typename E>
	void add(Real time, E &&event) {
		assert(time >= startTime);
		assert(time <= endTime);
		events.push( std::make_tuple(id++, time,
				std::forward<E>(event)) );
	}

};

/**
 * Handles timestepping with interleaved events.
 * @see QuantPDE::Event
 */
template <Index Dimension, bool Forward>
class EventIteration final : public EventIterationBase<Dimension, Forward> {

	virtual Vector transformIterand(const Vector &iterand) {
		Vector it = iterand;

		if(this->iteration == nullptr || this->iteration->done()) {

			// Perform events
			while( std::get<1>(this->events.top()) == time ) {
				it = ( std::get<2>(this->events.pop()) )(it);
			}

			this->clearHistory();

			Real t = this->events.empty() ? this->startTime
					: std::get<0>(this->events.top());
			this->iteration = this->factory->make(t, time);

		}

		return it;
	}

	virtual void clear() {
		time = this->endTime;
		this->iteration = nullptr;
	}

	virtual Real initialTime(Real) const {
		return this->endTime;
	}

	virtual bool done() const {
		return time <= this->startTime;
	}

public:

	/**
	 * Constructor.
	 */
	template <typename F>
	EventIteration(Real startTime, Real endTime, F &factory) noexcept
			: EventIterationBase<Dimension, Forward>(startTime,
			endTime, std::forward<F>(factory)) {
	}

};

template <Index Dimension>
class EventIteration<Dimension, true> final
		: public EventIterationBase<Dimension, true> {

	virtual Vector transformIterand(const Vector &iterand) {
		Vector it = iterand;

		if(this->iteration == nullptr || this->iteration->done()) {

			// Perform events
			while( std::get<1>(this->events.top()) == time ) {
				it = ( std::get<2>(this->events.pop()) )(it);
			}

			this->clearHistory();

			Real t = this->events.empty() ? this->endTime
					: std::get<0>(this->events.top());
			this->iteration = this->factory->make(time, t);

		}

		return it;
	}

	virtual void clear() {
		time = this->startTime;
		this->iteration = nullptr;
	}

	virtual Real initialTime(Real) const {
		return this->startTime;
	}

	virtual bool done() const {
		return time >= this->endTime;
	}

public:

	/**
	 * Constructor.
	 */
	template <typename F>
	EventIteration(Real startTime, Real endTime, F &factory) noexcept
			: EventIterationBase<Dimension, true>(startTime,
			endTime, std::forward<F>(factory)) {
	}

};

template <Index Dimension>
using ReverseEventIteration = EventIteration<Dimension, false>;

template <Index Dimension>
using ForwardEventIteration = EventIteration<Dimension, true>;

typedef ReverseEventIteration<1> ReverseEventIteration1;
typedef ReverseEventIteration<2> ReverseEventIteration2;
typedef ReverseEventIteration<3> ReverseEventIteration3;

typedef ForwardEventIteration<1> ForwardEventIteration1;
typedef ForwardEventIteration<2> ForwardEventIteration2;
typedef ForwardEventIteration<3> ForwardEventIteration3;

////////////////////////////////////////////////////////////////////////////////

/**
 * Steps through time at constant intervals.
 */
template <bool Forward>
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

	class Factory : public StepperFactory<Forward> {

		typedef std::unique_ptr<Iteration> I;

		unsigned steps;
		int lookback;

	public:

		/**
		 * Constructor.
		 */
		Factory(unsigned steps, int lookback = DEFAULT_LOOKBACK)
				noexcept : steps(steps), lookback(lookback) {
		}

		/**
		 * Copy constructor.
		 */
		Factory(const Factory &that) noexcept : steps(that.steps),
				lookback(that.lookback) {
		}

		/**
		 * Assignment operator.
		 */
		Factory &operator=(const Factory &that) & noexcept {
			steps = that.steps;
			lookback = that.lookback;
		}

		virtual I make(Real startTime, Real endTime) {
			return I(new ConstantStepper(startTime, endTime,
					steps, lookback));
		}

	};

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
template <bool Forward>
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

	class Factory : public StepperFactory<Forward> {

		typedef std::unique_ptr<Iteration> I;

		Real dt;
		Real target;
		Real scale;
		int lookback;

	public:

		/**
		 * Constructor.
		 */
		Factory(Real dt, Real target, Real scale,
				int lookback = DEFAULT_LOOKBACK) noexcept
				: dt(dt), target(target), scale(scale),
				lookback(lookback) {
		}

		/**
		 * Copy constructor.
		 */
		Factory(const Factory &that) noexcept : dt(that.dt),
				target(that.target), scale(that.scale),
				lookback(that.lookback) {
		}

		/**
		 * Assignment operator.
		 */
		Factory &operator=(const Factory &that) & noexcept {
			dt = that.dt;
			target = that.target;
			scale = that.scale;
			lookback = that.lookback;
		}

		virtual I make(Real startTime, Real endTime) {
			return I(new VariableStepper(startTime, endTime,
					dt, target, scale, lookback));
		}

	};

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

