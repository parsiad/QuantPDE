#ifndef QUANT_PDE_CORE_STEPPER
#define QUANT_PDE_CORE_STEPPER

#include <algorithm>  // std::min, std::max
#include <cstdlib>    // std::abs
#include <functional> // std::function
#include <memory>     // std::unique_ptr
#include <queue>      // std::priority_queue
#include <tuple>      // std::tuple
#include <vector>     // std::vector

namespace QuantPDE {

namespace Metafunctions {

namespace NaryFunctionSignatureHelpers {

template <typename R, typename ...Ts>
using TransformTarget = R (const Interpolant<sizeof...(Ts)> &, Ts...);

template <unsigned N>
using TransformSignature = Type<TransformTarget, Real, N, Real>;

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

	typedef typename InterpolantFactory<Dimension>::Wrapper In;
	typedef typename Map<Dimension>::Wrapper Out;

	Transform<Dimension> transform;

	In in;
	Out out;

	////////////////////////////////////////////////////////////////////////

	// TODO: Use C++1y move-capture in the future; return a bound function

	template <Index N, typename T, typename ...Ts>
	struct Driver : public Driver<N - 1, T, T, Ts...> {
	};

	template <typename T, typename ...Ts>
	struct Driver<0, T, Ts...> {
		static_assert( sizeof...(Ts) == Dimension,
				"The number of binding arguments in an event "
				"should be equal to the dimension");

		static inline Vector function(
			const Transform<Dimension> &transform,
			const Interpolant<Dimension> &interpolant,
			const Map<Dimension> &out
		) {
			return out([&] (Ts ...coordinates) {
				return transform(interpolant, coordinates...);
			});
		}
	};

	////////////////////////////////////////////////////////////////////////

	template <typename V>
	Vector _doEvent(V &&vector) const {
		// 1. Make interpolated function from vector
		// 2. Create a new the function by binding the first argument of
		//    the transform function to the interpolated function
		// 3. Map the function to the domain (producing a vector)

		return Driver<Dimension, Real>::function(
			transform,
			in.make( std::forward<V>(vector) ),
			out
		);
	}

	virtual Vector doEvent(const Vector &vector) const {
		return _doEvent(vector);
	}

	virtual Vector doEvent(Vector &&vector) const {
		return _doEvent(std::move(vector));
	}

public:

	/**
	 * Constructor.
	 * @param transform A function that transforms the solution.
	 * @param in An interpolant factory that handles how the solution is
	 *           interpolated off the domain nodes.
	 * @param out A map that handles how to transfer the transformed
	 *            solution back onto the domain.
	 */
	template <typename T>
	Event(
		T &&transform,
		In in,
		Out out
	) noexcept :
		transform(std::forward<T>(transform)),
		in(std::move(in)),
		out(std::move(out)) {
	}

	/**
	 * Default constructor for rectilinear grids.
	 * @param transform A function that transforms the solution.
	 * @param grid A rectilinear grid.
	 */
	template <typename T>
	Event(
		T &&transform,
		const RectilinearGrid<Dimension> &grid
	) noexcept :
		transform(std::forward<T>(transform)),
		in(grid.defaultInterpolantFactory()),
		out( std::unique_ptr<Map<Dimension>>(
				new PointwiseMap<Dimension>(grid)) )
	{
	}

	// Disable copy constructor and copy assignment operator.
	Event(const Event &) = delete;
	Event &operator=(const Event &) = delete;

	/**
	 * Move constructor.
	 */
	Event(Event &&that) noexcept :
		transform(std::move(that.transform)),
		in(std::move(that.in)),
		out(std::move(that.out))
	{
	}

	/**
	 * Move assignment operator.
	 */
	Event &operator=(Event &&that) & noexcept {
		transform = std::move(that.transform);
		in = std::move(that.in);
		out = std::move(that.out);
		return *this;
	}

};

typedef Event<1> Event1;
typedef Event<2> Event2;
typedef Event<3> Event3;

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

	typedef typename std::conditional<
		Forward,
		std::greater<Real>,
		std::less<Real>
	>::type Order;

	typedef std::tuple<
		int,
		Real,
		std::unique_ptr<EventBase>
	> T;

	struct TimeOrder {
		// Returns true if a goes before b in the ordering
		bool operator()(const T &a, const T &b) const {
			return
				Order()( std::get<1>(a), std::get<1>(b) )
				|| (
					std::get<1>(a) == std::get<1>(b)
					&& Order()( std::get<0>(a),
							std::get<0>(b) )
				)
			;
		}
	};

	// If events are set to the same time, ties are broken depending on the
	// order they were added. Events added later are assumed to occur later
	// in time (e.g. handled earlier if Forward == true; later otherwise).

	typedef std::unique_ptr<Iteration> I;

	////////////////////////////////////////////////////////////////////////

	int id;
	std::priority_queue<T, std::vector<T>, TimeOrder> events;

	Real startTime, endTime, time;

	const StepperFactory<Forward> *factory;
	I iteration;

private:

	virtual Real timestep() {
		Real step = iteration->timestep();
		time += step;
		return step;
	}

public:

	/**
	 * Constructor.
	 */
	template <typename F>
	EventIterationBase(Real startTime, Real endTime, F &factory,
			int lookback = DEFAULT_LOOKBACK) noexcept
			: Iteration(lookback), id(0), startTime(startTime),
			endTime(endTime), factory(&factory) {
		assert(startTime >= 0);
		assert(startTime < endTime);
	}

	/**
	 * Adds an event to be processed.
	 * @param time The time at which the event occurs.
	 * @param event The event.
	 */
	void add(Real time, std::unique_ptr<EventBase> event) {
		assert(time >= startTime);
		assert(time <= endTime);
		assert(time != initialTime(0.));

		events.emplace( id++, time, std::move(event) );
	}

	/**
	 * Adds an event to be processed.
	 * @param time The time at which the event occurs.
	 * @param args Arguments passed to Event constructor.
	 * @see QuantPDE::Event
	 */
	template <typename ...Ts>
	void add(Real time, Ts &&...args) {
		// TODO: Does not work with clang-503.0.40 with -std=c++1y.
		//       Tries to use Event's copy constructor.
		//       Is this a Clang bug?

		assert(time >= startTime);
		assert(time <= endTime);
		assert(time != initialTime(0.));

		events.emplace(
			id++,
			time,
			std::unique_ptr<EventBase>(
				new Event<Dimension>(
					std::forward<Ts>(args)...
				)
			)
		);
	}

};

#define QUANT_PDE_TMP(TRANSFORMED)                                             \
		do { if(this->iteration && !this->events.empty()) {            \
			Real t = std::get<1>( this->events.top() );            \
			do {                                                   \
				TRANSFORMED = true;                            \
				it = ( *std::get<2>(this->events.top()) )(it); \
				this->events.pop();                            \
				if(this->events.empty()) {                     \
					break;                                 \
				}                                              \
			} while( std::get<1>(this->events.top()) == t );       \
		} } while(0)

/**
 * Handles timestepping with interleaved events.
 * @see QuantPDE::Event
 */
template <Index Dimension, bool Forward>
class EventIteration final : public EventIterationBase<Dimension, Forward> {

	virtual bool transformIterand(Vector &it) {
		bool transformed = false;

		if(!this->iteration || this->iteration->done()) {
			QUANT_PDE_TMP(transformed);

			// Get next terminal time
			Real t = this->events.empty()
				? this->startTime
				: std::get<1>(this->events.top());

			if(!done()) {
				this->iteration = this->factory->make(
						t, this->time);
				this->iteration->clear();
			}
		}

		return transformed;
	}

	virtual void clear() {
		this->time = this->endTime;
		this->iteration = nullptr;
	}

	virtual Real initialTime(Real) const {
		return this->endTime;
	}

	virtual bool done() const {
		return this->time <= this->startTime;
	}

public:

	/**
	 * Constructor.
	 */
	template <typename F>
	EventIteration(
		Real startTime,
		Real endTime,
		F &factory,
		int lookback = DEFAULT_LOOKBACK
	) noexcept : EventIterationBase<Dimension, Forward>(
		startTime,
		endTime,
		factory,
		lookback
	) {
	}

};

template <Index Dimension>
class EventIteration<Dimension, true> final
		: public EventIterationBase<Dimension, true> {

	virtual bool transformIterand(Vector &it) {
		bool transformed = false;

		if(!this->iteration || this->iteration->done()) {
			QUANT_PDE_TMP(transformed);

			// Get next terminal time
			Real t = this->events.empty()
				? this->endTime
				: std::get<1>(this->events.top());

			if(!done()) {
				this->iteration = this->factory->make(
						this->time, t);
				this->iteration->clear();
			}
		}

		return transformed;
	}

	virtual void clear() {
		this->time = this->startTime;
		this->iteration = nullptr;
	}

	virtual Real initialTime(Real) const {
		return this->startTime;
	}

	virtual bool done() const {
		return this->time >= this->endTime;
	}

public:

	/**
	 * Constructor.
	 */
	template <typename F>
	EventIteration(
		Real startTime,
		Real endTime,
		F &factory,
		int lookback = DEFAULT_LOOKBACK
	) noexcept : EventIterationBase<Dimension, true>(
		startTime,
		endTime,
		factory,
		lookback
	) {
	}

};

#undef QUANT_PDE_TMP

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

		virtual I make(Real startTime, Real endTime) const {
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

	static_assert(std::numeric_limits<Real>::is_iec559,
			"IEEE 754 required");

	Real startTime, endTime, dt, target, scale, time;
	Real (VariableStepper::*_step)();

	inline void __step1() {
		const Real epsilon = 1e-6;

		const Vector
			&v1 = iterand(0),
			&v0 = iterand(1)
		;

		// Tested 2014-06-08
		const Real quotient = ( v1 - v0 ).cwiseAbs().cwiseQuotient(
			( scale * Vector::Ones( v1.size() ) ).cwiseMax(
				v1.cwiseAbs().cwiseMax(v0.cwiseAbs())
			)
		).maxCoeff();
		dt *= target / (quotient + epsilon);

		//std::cout << v0 << std::endl << std::endl << v1 << std::endl
		//		<< std::endl << quotient << " | " << dtOld
		//		<< " -> " << dt << std::endl << std::endl;
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

		Real dt, target, dtMin, dtMax, scale;
		int lookback;

	public:

		/**
		 * Constructor.
		 */
		Factory(
			Real dt,
			Real target,
			Real scale = 1.,
			int lookback = DEFAULT_LOOKBACK
		) noexcept :
			dt(dt),
			target(target),
			scale(scale),
			lookback(lookback)
		{
		}

		/**
		 * Copy constructor.
		 */
		Factory(const Factory &that) noexcept :
			dt(that.dt),
			target(that.target),
			scale(that.scale),
			lookback(that.lookback)
		{
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
			return I(
				new VariableStepper(
					startTime,
					endTime,
					dt,
					target,
					scale,
					lookback
				)
			);
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
		Real scale = 1.,
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

