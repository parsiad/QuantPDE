#ifndef QUANT_PDE_PROBLEM_HPP
#define QUANT_PDE_PROBLEM_HPP

#include <cassert>       // assert
#include <queue>         // std::priority_queue
#include <string>        // std::string
#include <tuple>         // std::tuple
#include <typeinfo>      // typeid
#include <unordered_map> // std::unordered_map
#include <unordered_set> // std::unordered_set
#include <vector>        // std::vector

namespace QuantPDE {

template <Index dim> // TODO: Check if positive
class Event {

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

////////////////////////////////////////////////////////////////////////////////

/**
 * Used to describe a problem abstractly. Such a problem description is not
 * inherently coupled with a method to solve it.
 */
template <Index dim, typename T = Real> // TODO: Check if positive
class Problem {

	template <typename U>
	using Dated = std::tuple<U, T>;

	////////////////////////////////////////////////////////////////////////

	std::vector< std::tuple<Constraint, T, T> > constraints;
	std::vector< std::tuple<Event<dim>, T   > > events;

	const Function<dim> initialCondition;

public:

	/**
	 * Constructor.
	 * @param initialCondition The initial condition (as a function).
	 */
	Problem(const Function<dim> &initialCondition)
			noexcept : initialCondition(initialCondition) {
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
	void add(E<dim> &&event, const T &start) {
		events.push_back(std::make_tuple( std::forward(event), start));
	}

	/**
	 * Adds a constraint to be processed.
	 * @param constraint The constraint.
	 * @param start The time at which this constraint goes into effect.
	 * @param end The time at which this constraint goes out of effect.
	 */
	template <typename C>
	void add(C &&constraint, const T &start, const T &end) {
		constraints.push_back(std::make_tuple( std::forward(constraint),
				start, end ));
	}

	template <Index N, typename U> friend class Solver;

};

template <Index dim, typename T = Real> // TODO: Check if positive
class Solver {

	template <typename U>
	using Dated = std::tuple<U, T>;

	template <typename U>
	struct Later {
		bool operator()(const U &a, const U &b) const {
			return std::get<1>(a) > std::get<1>(b);
		}
	};

	template <typename U>
	using MinPriorityQueue = std::priority_queue<U, std::vector<U>,
			Later<U>>;

	////////////////////////////////////////////////////////////////

	MinPriorityQueue< Dated<const Constraint *> > startQueue, endQueue;
	MinPriorityQueue< Dated<const Event<dim> *> > eventQueue;

	std::unordered_set<Constraint *> active;

	typedef std::function<void (const Constraint &)> Routine;
	std::unordered_map<std::string, Routine> routines;

	T earliest, now;

	T topmost() {
		T r;

		if(!startQueue.empty()) {
			r = std::get<1>( startQueue.top() );
		}

		if(!endQueue.empty()) {
			T t = std::get<1>( endQueue.top() );
			r = std::min(t, r);
		}

		if(!eventQueue.empty()) {
			T t = std::get<1>( eventQueue.top() );
			r = std::min(t, r);
		}

		return r;
	}

protected:

	Function<dim> solution;

public:

	/**
	 * Constructor.
	 */
	Solver(const Problem<dim, T> &problem) noexcept
			: startQueue(startQueue), endQueue(endQueue),
			eventQueue(eventQueue), earliest(topmost()),
			now(earliest), solution(problem.initialCondition) {

		for(const std::tuple<Constraint, T, T> &tuple
				: problem.constraints) {
			const Constraint *constraint = &std::get<0>(tuple);
			startQueue.push( std::make_tuple( constraint,
					std::get<1>(tuple) ) );
			endQueue.push( std::make_tuple( constraint,
					std::get<2>(tuple) ) );
		}

		for(const std::tuple<Event<dim>, T> &tuple : problem.events) {
			eventQueue.push( std::make_tuple( &std::get<0>(tuple),
				std::get<1>(tuple) ) );
		}

	}

	Solver(Problem<dim, T> &&problem) = delete;

	/**
	 * Destructor.
	 */
	virtual ~Solver() {
	}

	// Disable copy constructor and assignment operator
	Solver(const Solver &that) = delete;
	Solver &operator=(const Solver &that) & = delete;

	/**
	 * @return True if and only if the solution cannot be advanced
	 *         further.
	 */
	bool done() {
		return startQueue.empty() && endQueue.empty()
				&& eventQueue.empty();
	}

	/**
	 * Advance time.
	 */
	void advance() {
		assert(topmost() == now);

		while(!eventQueue.empty() && std::get<1>(
				eventQueue.top()) == now) {
			std::get<0>(eventQueue.pop())->handle(solution);
		}

		while(!endQueue.empty() && std::get<1>(
				endQueue.top()) == now) {
			active.erase( std::get<0>(endQueue.pop()) );
		}

		while(!startQueue.empty() && std::get<1>(
				startQueue.top()) == now) {
			active.insert( std::get<0>(startQueue.pop()) );
		}

		T next = topmost();

		beforeConstraints();
		for(const Constraint *constraint : active) {
			// next - now
			routines[constraint->identifier()](*constraint);
		}
		afterConstraints();
	}

	/**
	 * Advance time until the terminal time.
	 */
	void advanceUntilDone() {
		while(!done()) {
			advance();
		}
	}

	virtual void beforeConstraints() = 0;
	virtual void afterConstraints() = 0;

};

template <typename T = Real>
using Problem1 = Problem<1, T>;

template <typename T = Real>
using Problem2 = Problem<2, T>;

template <typename T = Real>
using Problem3 = Problem<3, T>;

}

#endif

