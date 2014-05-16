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
#define QUANT_PDE_IDENTIFIER_METHOD(class) virtual std::string identifier() \
		const { return QUANT_PDE_IDENTIFIER(class);  }

////////////////////////////////////////////////////////////////////////////////

/**
 * Used to describe an initial value problem abstractly. Such a problem
 * description is not inherently coupled with a method to solve it.
 */
template <Index dim, typename T = Real>
class Problem {

	static_assert(dim > 0, "Dimension must be positive");

	template <typename U>
	using Dated = std::tuple<U, T>;

	////////////////////////////////////////////////////////////////////////

	std::vector< std::tuple< std::unique_ptr<Constraint>, T, T > >
			constraints;
	std::vector< std::tuple< std::unique_ptr<Event<dim>>, T > > events;

	const Function<dim> initial;

public:

	/**
	 * Constructor.
	 * @param initial The initial condition (as a function).
	 */
	Problem(const Function<dim> &initial) noexcept : initial(initial) {
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
	template <template <Index> class E>
	void add(std::unique_ptr<Event<dim>> event, const T &start) {
		events.push_back(std::make_tuple(std::move(event), start));
	}

	/**
	 * Adds a constraint to be processed.
	 * @param constraint The constraint.
	 * @param start The time at which this constraint goes into effect.
	 * @param end The time at which this constraint goes out of effect.
	 */
	void add(std::unique_ptr<Constraint> constraint, const T &start,
			const T &end) {
		constraints.push_back(std::make_tuple(std::move(constraint),
				start, end));
	}

	/**
	 * @return The initial condition.
	 */
	const Function<dim> &initialCondition() {
		return initial;
	}

};

template <typename T = Real>
using Problem1 = Problem<1, T>;

template <typename T = Real>
using Problem2 = Problem<2, T>;

template <typename T = Real>
using Problem3 = Problem<3, T>;

}

#endif

