#ifndef QUANT_PDE_PROBLEM_HPP
#define QUANT_PDE_PROBLEM_HPP

#include <memory>        // std::unique_ptr
#include <string>        // std::string
#include <typeinfo>      // typeid
#include <vector>        // std::vector

namespace QuantPDE {

template <Index dim>
class Event {

	static_assert(dim > 0, "Dimension must be positive");

public:

	Event() {
	}

	virtual ~Event() {
	}

	Event(const Event &) = delete;
	Event &operator=(const Event &) & = delete;

	virtual Function<dim> advance(const Function<dim> &solution) const = 0;

};

class Constraint {

public:

	Constraint() {
	}

	virtual ~Constraint() {
	}

	Constraint(const Constraint &) = delete;
	Constraint &operator=(const Constraint &) & = delete;

	virtual std::string identifier() const = 0;

};

#define QUANT_PDE_IDENTIFIER(class) typeid(class).name()

////////////////////////////////////////////////////////////////////////////////

/**
 * Used to describe an initial value problem abstractly. Such a problem
 * description is not inherently coupled with a method to solve it.
 */
template <Index dim>
class Problem {

	static_assert(dim > 0, "Dimension must be positive");

	////////////////////////////////////////////////////////////////////////

	std::vector< std::tuple< std::unique_ptr<Constraint>, Real, Real > >
			constraints;
	std::vector< std::tuple< std::unique_ptr<Event<dim>>, Real > > events;

	const Function<dim> initial;

public:

	/**
	 * Constructor.
	 * @param initial The initial condition (as a function).
	 */
	template <typename F>
	Problem(F &&initial) noexcept : initial( std::forward<F>(initial) ) {
	}

	/**
	 * Destructor.
	 */
	virtual ~Problem() {
	}

	// Disable copy constructor and assignment operator
	Problem(const Problem &) = delete;
	Problem &operator=(const Problem &) & = delete;

	/**
	 * Adds an event to be processed.
	 * @param event The event.
	 * @param start The time at which the event occurs.
	 */
	void add(std::unique_ptr<Event<dim>> event, const Real &start) {
		events.push_back(std::make_tuple(std::move(event), start));
	}

	/**
	 * Adds a constraint to be processed.
	 * @param constraint The constraint.
	 * @param start The time at which this constraint goes into effect.
	 * @param end The time at which this constraint goes out of effect.
	 */
	void add(std::unique_ptr<Constraint> constraint, const Real &start,
			const Real &end) {
		constraints.push_back(std::make_tuple(std::move(constraint),
				start, end));
	}

	template <Index N> friend class Solver;

};

typedef Problem<1> Problem1;
typedef Problem<2> Problem2;
typedef Problem<3> Problem3;

}

#endif

