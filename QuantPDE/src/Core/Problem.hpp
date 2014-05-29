#ifndef QUANT_PDE_CORE_PROBLEM
#define QUANT_PDE_CORE_PROBLEM

#include <limits>        // std::numeric_limits
#include <memory>        // std::unique_ptr
#include <string>        // std::string
#include <typeinfo>      // typeid
#include <unordered_set> // std::unordered_set
#include <vector>        // std::vector

namespace QuantPDE {

/**
 * Let $\u$ denote a solution. An event is used to describe the relation
 * \f$u\left(\mathbf{x},t^+\right)=f\left(\mathbf{x},t,u\left(\cdot,t^-\right)\right)\f$.
 * Here, the event describes \f$f\f$.
 */
template <Index Dimension>
class Event {

	static_assert(Dimension > 0, "Dimension must be positive");

public:

	/**
	 * Constructor.
	 */
	Event() noexcept {
	}

	/**
	 * Destructor.
	 */
	virtual ~Event() {
	}

	// Disable copy constructor and assignment operator
	Event(const Event &) = delete;
	Event &operator=(const Event &) & = delete;

	/**
	 * Transforms the solution across the event.
	 * @param solution The solution.
	 * @return The new solution.
	 */
	virtual Function<Dimension> transform(const Function<Dimension>
			&solution) const = 0;

};

////////////////////////////////////////////////////////////////////////////////

/**
 * A constraint that holds over an interval (in time).
 */
class Constraint {

public:

	/**
	 * Constructor.
	 */
	Constraint() noexcept {
	}

	/**
	 * Destructor.
	 */
	virtual ~Constraint() {
	}

	// Disable copy constructor and assignment operator.
	Constraint(const Constraint &) = delete;
	Constraint &operator=(const Constraint &) & = delete;

	/**
	 * @return An identifier unique to the type of constraint.
	 */
	virtual std::string identifier() const = 0;

};

typedef std::unordered_set<const Constraint *> ConstraintSet;

#define QUANT_PDE_IDENTIFIER(DERIVED_CLASS) typeid(DERIVED_CLASS).name()

////////////////////////////////////////////////////////////////////////////////

/**
 * Used to describe an initial value problem abstractly. Such a problem
 * description is not inherently coupled with a method to solve it.
 * @param InitialAtLeft Whether the initial condition occurs at the left
 *                      endpoint or the right one.
 */
template <Index Dimension, bool InitialAtLeft = false>
class Problem {

	static_assert(Dimension > 0, "Dimension must be positive");

	////////////////////////////////////////////////////////////////////////

	typedef std::vector<std::tuple<std::unique_ptr<const Constraint>, Real,
			Real>> C;
	typedef  std::vector<std::tuple<std::unique_ptr<const Event<Dimension>>,
			Real>> E;

	C c;
	E e;

	const Function<Dimension> initial;

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
	 * @param time The time at which the event occurs.
	 */
	void add(std::unique_ptr<const Event<Dimension>> event,
			const Real &time) {
		assert(time >= 0);

		e.push_back(std::make_tuple(std::move(event), time));
	}

	/**
	 * Adds a constraint to be processed.
	 * @param constraint The constraint.
	 * @param startTime The time at which this constraint goes into effect.
	 * @param endTime The time at which this constraint goes out of effect.
	 */
	void add(std::unique_ptr<const Constraint> constraint,
			const Real &startTime, const Real &endTime) {
		assert(startTime >= 0);
		assert(endTime > startTime);

		c.push_back(std::make_tuple(std::move(constraint),
				startTime, endTime));
	}

	/**
	 * @return A vector of the constraints in this problem.
	 */
	const C &constraints() {
		return c;
	}

	/**
	 * @return A vector of the events in this problem.
	 */
	const E &events() {
		return e;
	}

	/**
	 * @return The initial condition.
	 */
	const Function<Dimension> &initialCondition() {
		return initial;
	}

};

typedef Problem<1> Problem1;
typedef Problem<2> Problem2;
typedef Problem<3> Problem3;

}

#endif

