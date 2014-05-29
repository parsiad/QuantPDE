#ifndef QUANT_PDE_CORE_SOLVER
#define QUANT_PDE_CORE_SOLVER

#include <cassert>       // assert
#include <memory>        // std::unique_ptr
#include <queue>         // std::priority_queue
#include <tuple>         // std::tuple
#include <type_traits>   // std::conditional
#include <vector>        // std::vector

namespace QuantPDE {

/**
 * Used to solve (approximately) a problem description.
 * @see QuantPDE::Problem
 */
template <Index Dimension, bool InitialAtLeft = false>
class Solver {

	static_assert(Dimension > 0, "Dimension must be positive");

	template <typename U>
	using Dated = std::tuple<U, Real>;

	template <typename U>
	struct Earlier {
		bool operator()(const U &a, const U &b) const {
			return std::get<1>(a) < std::get<1>(b);
		}
	};

	template <typename U>
	struct Later {
		bool operator()(const U &a, const U &b) const {
			return std::get<1>(a) > std::get<1>(b);
		}
	};

	template <typename U>
	using Order = typename std::conditional<InitialAtLeft, Later<U>,
			Earlier<U>>::type;

	template <typename U>
	using MinPriorityQueue = std::priority_queue<U, std::vector<U>,
			Order<U>>;

	////////////////////////////////////////////////////////////////////////

	MinPriorityQueue< Dated<const Constraint *> > startQueue, endQueue;
	MinPriorityQueue< Dated<const Event<Dimension> *> > eventQueue;

	ConstraintSet active;

	//std::queue< std::unique_ptr<Stepper<Dimension, InitialAtLeft>> >
	//		steppers;

	Real now;

	////////////////////////////////////////////////////////////////////////

	Function<Dimension> solution;

	////////////////////////////////////////////////////////////////////////

	Real topmost() {
		Real r;

		if(!startQueue.empty()) {
			r = std::get<1>( startQueue.top() );
		}

		if(!endQueue.empty()) {
			Real t = std::get<1>( endQueue.top() );
			r = std::max(t, r);
		}

		if(!eventQueue.empty()) {
			Real t = std::get<1>( eventQueue.top() );
			r = std::max(t, r);
		}

		return r;
	}

	bool done() {
		return startQueue.empty() && endQueue.empty()
				&& eventQueue.empty();
	}

protected:

	/**
	 * Pushes a stepper onto the queue. Each time a set of active
	 * constraints is enconutered, a stepper is popped from the queue to
	 * handle the set and advance the solution.
	 * @param stepper The stepper.
	 */
	//void add(std::unique_ptr< Stepper<Dimension, InitialAtLeft> > stepper) {
	//	steppers.push( std::move(stepper) );
	//}

public:

	/**
	 * Constructor.
	 */
	Solver(const Problem<Dimension, InitialAtLeft> &problem) noexcept
			: now( topmost() ), solution(problem.initialCondition())
			{

		for(auto t : problem.constraints()) {

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

		for(auto t : problem.events()) {

			eventQueue.push(std::make_tuple(
				std::get<0>(t).get(),
				std::get<1>(t)
			));

		}

	}

	// Disable rvalue constructors
	Solver(Problem<Dimension, InitialAtLeft> &&problem) = delete;

	/**
	 * Destructor.
	 */
	virtual ~Solver() {
	}

	// Disable copy constructor and assignment operator
	Solver(const Solver &that) = delete;
	Solver &operator=(const Solver &that) & = delete;

	/**
	 * Advance time.
	 */
	void advance() {
		assert(!done());
		assert(topmost() == now);

		while(!endQueue.empty() && std::get<1>(
				endQueue.top()) == now) {
			if(InitialAtLeft) { // Unroll conditional
				active.erase ( std::get<0>(endQueue.pop()) );
			} else {
				active.insert( std::get<0>(endQueue.pop()) );
			}
		}

		while(!startQueue.empty() && std::get<1>(
				startQueue.top()) == now) {
			if(InitialAtLeft) { // Unroll conditional
				active.insert( std::get<0>(startQueue.pop()) );
			} else {
				active.erase ( std::get<0>(startQueue.pop()) );
			}
		}

		// If the active set is empty, this is either an isolated event
		// or we have advanced too many times
		assert(!active.empty());

		Real next = topmost();
		//solution = steppers.pop()->transform(solution, active, now,
		//		next);
		now = next;

		while(!eventQueue.empty() && std::get<1>(
				eventQueue.top()) == now) {
			solution = std::get<0>(eventQueue.pop())->transform(
					solution);
		}
	}

	/**
	 * Advance time until the terminal time.
	 */
	void advanceUntilDone() {
		while(!done()) {
			advance();
		}
	}

	/**
	 * @return The current solution. Calling this immediately after
	 *         initializing the solver is guaranteed to return the initial
	 *         condition associated with the problem.
	 */
	const Function<Dimension> &currentSolution() const {
		return solution;
	}

};

}

#endif

