#ifndef QUANT_PDE_STEPPER_HPP
#define QUANT_PDE_STEPPER_HPP

namespace QuantPDE {

template <Index dim>
class Stepper {

};

typedef Stepper<1> Stepper1;
typedef Stepper<2> Stepper2;
typedef Stepper<3> Stepper3;

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

