#ifndef QUANT_PDE_SOLVER_HPP
#define QUANT_PDE_SOLVER_HPP

#include <cassert>       // assert
#include <queue>         // std::priority_queue
#include <tuple>         // std::tuple
#include <unordered_map> // std::unordered_map
#include <unordered_set> // std::unordered_set
#include <vector>        // std::vector

namespace QuantPDE {

template <Index dim>
class Solver {

	static_assert(dim > 0, "Dimension must be positive");

	template <typename U>
	using Dated = std::tuple<U, Real>;

	template <typename U>
	struct Later {
		bool operator()(const U &a, const U &b) const {
			return std::get<1>(a) > std::get<1>(b);
		}
	};

	template <typename U>
	using MinPriorityQueue = std::priority_queue<U, std::vector<U>,
			Later<U>>;

	////////////////////////////////////////////////////////////////////////

	MinPriorityQueue< Dated<const Constraint *> > startQueue, endQueue;
	MinPriorityQueue< Dated<const Event<dim> *> > eventQueue;

	std::unordered_set<const Constraint *> active;

	typedef std::function<void (const Constraint &, Real, Real)> Routine;
	std::unordered_map<std::string, Routine> routines;

	Real earliest, now;

	const Stepper<dim> *stepper;

	Real topmost() {
		Real r;

		if(!startQueue.empty()) {
			r = std::get<1>( startQueue.top() );
		}

		if(!endQueue.empty()) {
			Real t = std::get<1>( endQueue.top() );
			r = std::min(t, r);
		}

		if(!eventQueue.empty()) {
			Real t = std::get<1>( eventQueue.top() );
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

		Real next = topmost();

		beforeConstraints();
		//do {
			//Real dt = nextStep();
			for(const Constraint *constraint : active) {
				routines[constraint->identifier()](*constraint,
						now, next);
			}
		//} while(stepping());
		afterConstraints();
	}

protected:

	const Grid<dim> &grid;
	Function<dim> solution;

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
	Solver(const Problem<dim> &problem, const Grid<dim> &grid,
			const Stepper<dim> &stepper) noexcept
			: earliest(topmost()), now(earliest), stepper(&stepper),
			grid(grid), solution(problem.initial) {

		for(const std::tuple< std::unique_ptr<Constraint>, Real,
				Real > &t : problem.constraints) {

			const Constraint *constraint = std::get<0>(t).get();

			startQueue.push(std::make_tuple(
				constraint,
				std::get<1>(t)
			));

			endQueue.push(std::make_tuple(
				constraint,
				std::get<2>(t)
			));

		}

		for(const std::tuple< std::unique_ptr<Event<dim>>, Real> &tuple
				: problem.events) {

			eventQueue.push(std::make_tuple(
				std::get<0>(tuple).get(),
				std::get<1>(tuple)
			));

		}

	}

	// Disable rvalue constructors
	Solver(const Problem<dim> &, const Grid<dim> &, Stepper<dim> &&)
			= delete;
	Solver(const Problem<dim> &, Grid<dim> &&, const Stepper<dim> &)
			= delete;
	Solver(const Problem<dim> &, Grid<dim> &&, Stepper<dim> &&) = delete;
	Solver(Problem<dim> &&, const Grid<dim> &, const Stepper<dim> &)
			= delete;
	Solver(Problem<dim> &&, const Grid<dim> &, Stepper<dim> &&) = delete;
	Solver(Problem<dim> &&, Grid<dim> &&, const Stepper<dim> &) = delete;
	Solver(Problem<dim> &&, Grid<dim> &&, Stepper<dim> &&) = delete;

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

#define QUANT_PDE_REGISTER_ROUTINE(class)                      \
	registerRoutine(                                       \
		QUANT_PDE_IDENTIFIER(class),                   \
		[this] (const Constraint &a, Real b, Real c) { \
			_##class(a, b, c);                     \
		}                                              \
	)

typedef Solver<1> Solver1;
typedef Solver<2> Solver2;
typedef Solver<3> Solver3;

}

#endif

