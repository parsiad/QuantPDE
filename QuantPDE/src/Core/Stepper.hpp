#ifndef QUANT_PDE_CORE_STEPPER
#define QUANT_PDE_CORE_STEPPER

#include <functional>    // std::function
#include <string>        // std::string
#include <unordered_map> // std::unordered_map

namespace QuantPDE {

/**
 * Used as a wrapper around solution routines to allow for caching.
 */
template <typename S>
class Cache {

	bool solved;
	S s;

	virtual S solve() = 0;

public:

	Cache() : solved(false) {
	}

	Cache(const Cache &) = delete;
	Cache &operator=(const Cache &) & = delete;

	const S &solution() {
		if(solved) {
			return s;
		}

		s = solve();
		solved = true;
		return s;
	}

};

// TODO: Document

////////////////////////////////////////////////////////////////////////////////

namespace Metafunctions {

namespace NaryFunctionSignatureHelpers {

template <typename R, typename... Ts>
using RoutineTarget = R (const Constraint &, Real, Ts...);

template<unsigned N>
using RoutineSignature = Type<RoutineTarget, Matrix, N, Vector>;

} // NaryFunctionSignatureHelpers

} // Metafunctions

template <Index N>
using Routine = std::function< Metafunctions::NaryFunctionSignatureHelpers
		::RoutineSignature<N> >;

typedef Routine<0> Routine0;
typedef Routine<1> Routine1;
typedef Routine<2> Routine2;
typedef Routine<3> Routine3;

/**
 * Delegates to each constraint a routine to handle it.
 * @see QuantPDE::Constraint
 */
template <Index Controls = 0>
class Discretizer {

	static_assert(Controls >= 0, "Number of controls must be nonnegative");

	std::unordered_map<std::string, Routine<Controls> > routines;

protected:

	/**
	 * Associates a routine with a type of constraint.
	 * @param identifier The identifier associated with the type of
	 *                   constraint.
	 * @param routine The routine.
	 */
	template <typename F>
	void registerRoutine(const std::string &identifier, F &&routine) {
		routines[identifier] = routine;
	}

public:

	/**
	 * Constructor.
	 */
	Discretizer() noexcept {
	}

	/**
	 * Destructor.
	 */
	virtual ~Discretizer() {
	}

	// Disable copy constructor and assignment operator
	Discretizer(const Discretizer &that) = delete;
	Discretizer &operator=(const Discretizer &that) & = delete;

	template <typename ...Args>
	Matrix handle(const Constraint &constraint, Real time, Args ...controls)
			{
		return routines[constraint.identifier()](constraint, time,
				controls...);
	}

};

typedef Discretizer<0> Discretizer0;
typedef Discretizer<1> Discretizer1;
typedef Discretizer<2> Discretizer2;
typedef Discretizer<3> Discretizer3;

#define QUANT_PDE_REGISTER_ROUTINE(class)                      \
	registerRoutine(                                       \
		QUANT_PDE_IDENTIFIER(class),                   \
		[this] (const Constraint &a, Real b, Real c) { \
			_##class(a, b, c);                     \
		}                                              \
	)

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

/*
template <typename E>
class ConstantStepper : public Stepper<E> {

	unsigned N, n;
	double start, dt;

	virtual E next() {
		E e = next( start + dt * n, start + dt * (n + 1) );
		n++;
		return e;
	}

	virtual E next(double previousTime, double nextTime) = 0;

public:

	// TODO: Document
	ConstantStepper(const Grid &grid, Vector &initial, double start,
			double end, unsigned steps)
			: Stepper<E>(grid, initial, start), N(steps), n(0),
			start(start), dt((end - start) / steps) {
		assert(N > 0);
	}

	// TODO: Copy constructor

	virtual bool done() const {
		return n >= N;
	}

};

// TODO: Document
template <typename E>
class VariableStepper : public Stepper<E> {

	double end, dt, dnorm, scale;

	virtual E next() {
		double nextTime = this->currentTime() + dt;
		if(nextTime > end) {
			nextTime = end;
		}
		return next(this->currentTime(), nextTime);
	}

	virtual void postProcess(double previousTime, double nextTime,
			const Vector &previousSolution,
			const Vector &nextSolution) {
		Vector divisor = (scale * this->grid.ones()).cwiseMax(
				nextSolution.cwiseAbs().cwiseMax(
				previousSolution.cwiseAbs()));

		Vector relative = (previousSolution - nextSolution)
				.cwiseAbs().cwiseQuotient(divisor);

		dt *= dnorm / relative.maxCoeff();

		// TODO: Not sure if this currently works. Fix it.
	}


	virtual E next(double previousTime, double nextTime) = 0;

public:

	// TODO: Document
	VariableStepper(const Grid &grid, Vector &initial, double start,
			double end, double dt, double dnorm, double scale = 1)
			: Stepper<E>(grid, initial, start),
			end(end), dt(dt), dnorm(dnorm), scale(scale) {
	}

	// TODO: Copy constructor

	virtual bool done() const {
		return this->currentTime() >= end;
	}

};
*/

}

#endif

