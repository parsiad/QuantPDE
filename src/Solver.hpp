#ifndef QUANT_PDE_SOLVER_HPP
#define QUANT_PDE_SOLVER_HPP

#include <cassert>       // assert
#include <queue>         // std::priority_queue
#include <tuple>         // std::tuple
#include <unordered_map> // std::unordered_map
#include <unordered_set> // std::unordered_set
#include <vector>        // std::vector

namespace QuantPDE {

template <Index dim, typename T = Real>
class Solver {

	static_assert(dim > 0, "Dimension must be positive");

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

	std::unordered_set<const Constraint *> active;

	typedef std::function<void (const Constraint &, T, T)> Routine;
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

	void advance() {
		assert(!done());
		assert(topmost() == now);

		while(!eventQueue.empty() && std::get<1>(
				eventQueue.top()) == now) {
			solution = std::get<0>(eventQueue.pop())->advance(
					solution);
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
			routines[constraint->identifier()](*constraint, now,
					next);
		}
		afterConstraints();
	}

protected:

	const Grid<dim> &grid;
	Function<dim> solution;

	/**
	 * Run before anything has been done.
	 */
	virtual void onStart() {
	}

	/**
	 * Run before constraints have been processed.
	 */
	virtual void beforeConstraints() {
	}

	/**
	 * Run after constraints have been processed.
	 */
	virtual void afterConstraints() {
	}

	/**
	 * Run when time cannot be advanced further.
	 */
	virtual void onFinish() {
	}

public:

	/**
	 * Constructor.
	 */
	Solver(const Problem<dim, T> &problem, const Grid<dim> &grid) noexcept
			: startQueue(startQueue), endQueue(endQueue),
			eventQueue(eventQueue), earliest(topmost()),
			now(earliest), grid(&grid),
			solution(problem.initialCondition()) {

		for(const std::tuple<Constraint, T, T> &tuple
				: problem.constraints) {
			const Constraint *constraint = std::get<0>(tuple).get();
			startQueue.push( std::make_tuple( constraint,
					std::get<1>(tuple) ) );
			endQueue.push( std::make_tuple( constraint,
					std::get<2>(tuple) ) );
		}

		for(const std::tuple<Event<dim>, T> &tuple : problem.events) {
			eventQueue.push(std::make_tuple(
				std::get<0>(tuple).get(), std::get<1>(tuple) ));
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
	 * Associates a routine with a family of constraints.
	 * @param identifier The identifier associated with the family of
	 *                   constraints.
	 * @param routine The routine.
	 */
	template <typename F>
	void registerRoutine(const std::string &identifier, F &&routine) {
		routines[identifier] = routine;
	}

	/**
	 * @return True if and only if the solution cannot be advanced
	 *         further.
	 */
	bool done() {
		return startQueue.empty() && endQueue.empty()
				&& eventQueue.empty();
	}

	/**
	 * Advance time until the terminal time.
	 */
	void advanceUntilDone() {
		onStart();
		while(!done()) {
			advance();
		}
		onFinish();
	}

	/**
	 * Returns the solution.
	 */
	const Function<dim> &currentSolution() const {
		return solution;
	}

};

template <typename T = Real>
using Solver1 = Solver<1, T>;

template <typename T = Real>
using Solver2 = Solver<2, T>;

template <typename T = Real>
using Solver3 = Solver<3, T>;

}

#endif

