#ifndef QUANT_PDE_EVENT_HPP
#define QUANT_PDE_EVENT_HPP

namespace QuantPDE {

/**
 * Used to advance the solution over time.
 */
class Event {

	const double previous, next;

public:

	/**
	 * Constructor.
	 * @param previousTime The time before the event.
	 * @param nextTime The time after the event.
	 */
	Event(double previousTime, double nextTime) : previous(previousTime),
			next(nextTime) {
	}

	/**
	 * Copy constructor.
	 */
	Event(const Event &that) : previous(that.previous), next(that.next) {
	}

	/**
	 * Advances the solution past this event.
	 * @param grid The grid the solution is on.
	 * @param solution The solution at the current time.
	 * @param previousTime The time before the event.
	 * @return The new solution.
	 */
	virtual Vector advance(const Grid &grid, const Vector &solution)
			const = 0;

	/**
	 * @return The time before the event.
	 */
	double previousTime() const {
		return previous;
	}

	/**
	 * @return The time after the event.
	 */
	double nextTime() const {
		return next;
	}

};

}

#endif

