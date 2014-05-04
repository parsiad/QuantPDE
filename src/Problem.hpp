#ifndef QUANT_PDE_FEATURE_HPP
#define QUANT_PDE_FEATURE_HPP

// TODO: This file has not been tested at all

#include <cassert>       // assert
#include <queue>         // std::priority_queue
#include <tuple>         // std::tuple
#include <unordered_set> // std::unordered_set
#include <vector>        // std::vector

namespace QuantPDE {

class Problem final {

	template <typename T>
	using Dated = std::tuple<T, DateTime>;

	template <typename T>
	class Greater {
		bool operator()(const Dated<T> &a, const Dated<T> &b) {
			return std::get<1>(a) > std::get<1>(b);
		}
	};

	typedef Dated<Constraint> A;
	typedef Dated<ExplicitConstraint> B;

	////////////////////////////////////////////////////////////////////////

	DateTime dateTime;

	std::priority_queue<A, std::vector<A>, Greater<A>> startQueue, endQueue;
	std::priority_queue<B, std::vector<B>, Greater<B>> explicitQueue;

	std::unordered_set<const Constraint *> active;

	const Grid *grid;
	const Function *initial;

	#ifndef NDEBUG
	bool started = false;
	#endif

public:

	/**
	 * Constructor.
	 * @param grid The underlying grid.
	 * @param initial The initial solution.
	 */
	Problem(const Grid &grid, const Function &initial) : grid(&grid),
			initial(&initial) {
	}

	// Disable copy constructor and assignment operator
	Problem(const Problem &) = delete;
	Problem &operator=(const Problem &) = delete;

	void add(const ExplicitConstraint &constraint,
			const DateTime &dateTime) {
		assert(!started);
		explicitQueue.push( std::make_tuple(constraint, dateTime) );
	}

	void add(const Constraint &constraint, const DateTime &start,
			const DateTime &end) {
		assert(!started);
		startQueue.push( std::make_tuple(constraint, start) );
		endQueue.push( std::make_tuple(constraint, end) );
	}

	bool done() const {
		return start.empty() && end.empty();
	}

	void process() {
		#ifndef NDEBUG
		started = true;
		#endif

		// Check invariants
		assert(startQueue.empty() ||
				std::get<1>(startQueue.top()) >= dateTime);
		assert(endQueue.empty() ||
				std::get<1>(endQueue.top()) >= dateTime);
		// At least one of the above should hold with equality

		// Remove elements from explicit queue
		while(!explicitQueue.empty()
				&& std::get<1>(explicitQueue.top())
				== dateTime) {
			// TODO: Apply explicit constraint
			explicitQueue.pop();
		}

		// Remove elements from end queue
		while(!endQueue.empty()
				&& std::get<1>(endQueue.top()) == dateTime) {
			active.erase( &std::get<0>(endQueue.pop()) );
		}

		// Remove elements from start queue
		while(!startQueue.empty()
				&& std::get<1>(startQueue.top()) == dateTime) {
			active.insert( &std::get<0>(startQueue.pop()) );
		}

		Date nextDateTime;
		if(startQueue.empty()) {
			if(!endQueue.empty()) {
				nextDateTime = std::get<1>(endQueue.top());
			} else {
				// No need to do anything; date is already at
				// the expiry
			}
		} else if(endQueue.empty()) {
			nextDateTime = std::get<1>(startQueue.top());
		} else {
			// Both empty
			nextDateTime = std::min( std::get<1>(endQueue.top()),
					std::get<1>(startQueue.top()) );
		}

		// TODO: Act on active set on dateTime -> nextDateTime
		for(const Constraint *constraint : active) {
			// Take sum of constraints and apply them
		}

		// Update date
		dateTime = nextDateTime;
	}

	void processAll() {
		while(!done()) {
			process();
		}
	}

};

}

#endif

