#ifndef QUANT_PDE_STEPPER_HPP
#define QUANT_PDE_STEPPER_HPP

namespace QuantPDE {

/**
 * An abstract interface to generate events that advance the solution until
 * some condition has been met.
 * @see QuantPDE::Event
 */
template <typename E>
class Stepper {

	Vector s;
	double t;
	unsigned m;

	virtual E next() = 0;

	virtual void postProcess(const E &e, const Vector &previousSolution,
			const Vector &nextSolution) {
	}

protected:

	/**
	 * The underlying grid.
	 */
	const Grid &grid;

public:

	/**
	 * Constructor.
	 * @param grid The underlying grid.
	 * @param initial The initial solution.
	 * @param start The time at which the initial solution is defined at.
	 */
	Stepper(const Grid &grid, Vector initial, double start = 0)
			: s(initial), t(start), m(0), grid(grid) {
	}

	/**
	 * Copy constructor.
	 */
	Stepper(const Stepper &that) : s(that.s), t(that.t), m(that.m),
			grid(that.grid) {
	}

	/**
	 * @return True if and only if all events have been popped.
	 */
	virtual bool done() const = 0;

	/**
	 * Processes the next event.
	 */
	void process() {
		assert(!done());

		E e = next();

		Vector newSolution = e.advance(grid, s);
		postProcess(e, s, newSolution);
		s = newSolution;

		t = e.nextTime();

		m++;
	}

	/**
	 * Processes all remaining events.
	 */
	void processAll() {
		while(!done()) {
			process();
		}
	}

	/**
	 * @return The current solution.
	 */
	const Vector &solution() {
		return s;
	}

	/**
	 * @return The current solution time.
	 */
	double currentTime() const {
		return t;
	}

	/**
	 * @return The number of events processed so far.
	 */
	unsigned steps() const {
		return m;
	}

};

// TODO: Document
template <typename E>
class ConstantStepper : public Stepper<E> {

	unsigned N, n;
	double start, dt;

	virtual E next() {
		E e = next( start + dt * n, start + dt * (n + 1) );
		n++;
		return e;
	}

	/**
	 * @param nextTime The time after the event.
	 * @return An event.
	 */
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

	/**
	 * @param nextTime The time after the event.
	 * @return An event.
	 */
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

}

#endif

