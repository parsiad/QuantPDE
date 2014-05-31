#ifndef QUANT_PDE_ITERATIVE_METHOD
#define QUANT_PDE_ITERATIVE_METHOD

#include <deque>   // std::deque
#include <list>    // std::list
#include <memory>  // std::unique_ptr
#include <utility> // std::forward, std::move

namespace QuantPDE {

/**
 * A pure virtual class representing a solver for equations of type \f$Ax=b\f$.
 */
class LinearSolver {

public:

	/**
	 * Constructor.
	 */
	LinearSolver() noexcept {
	}

	/**
	 * Destructor.
	 */
	virtual ~LinearSolver() {
	}

	// Disable copy constructor and assignment operator.
	LinearSolver(const LinearSolver &) = delete;
	LinearSolver &operator=(const LinearSolver &) & = delete;

	/**
	 * Initializes the linear solver with a matrix. If solving a linear
	 * system with a constant left-hand-side multiple times, this call
	 * should occur only once.
	 * @param A The left-hand-side matrix.
	 */
	virtual void initialize(const Matrix &A) = 0;

	/**
	 * Solves the linear system. This should only be called after a call
	 * to initialize.
	 * @param b The right-hand-side.
	 * @param guess An initial guess (ignored for noniterative methods).
	 * @return The solution.
	 * @see QuantPDE::LinearSolver::initialize
	 */
	virtual Vector solve(const Vector &b, const Vector &guess) const = 0;

};

/**
 * Solves \f$Ax=b\f$ with BiCGSTAB.
 */
class BiCGSTABSolver : public LinearSolver {

	BiCGSTAB solver;

public:

	virtual void initialize(const Matrix &A) {
		solver.compute(A);
	}

	BiCGSTABSolver() noexcept : LinearSolver() {
	}

	virtual Vector solve(const Vector &b, const Vector &guess) const {
		return solver.solveWithGuess(b, guess);
	}

};

////////////////////////////////////////////////////////////////////////////////

/**
 * Used to generate the left and right-hand sides of the linear system at each
 * iteration.
 */

class IterativeMethod;

class Linearizer {

	// Nonownership
	const IterativeMethod *method;

protected:

	const std::deque<Vector> &iterands() const;

public:

	template <typename M>
	Linearizer(M &method);

	virtual ~Linearizer() {
	}

	// Disable copy constructor and assignment operator.
	Linearizer(const Linearizer &) = delete;
	Linearizer &operator=(const Linearizer &) & = delete;

	virtual void onIterationStart() = 0;

	virtual bool doesAChange() const = 0;

	virtual Matrix A() const = 0;

	virtual Vector b() const = 0;

};

class IterativeMethod {

	unsigned lookback;
	std::deque<Vector> iterands;

	// Nonownership
	std::list<Linearizer *> linearizers;

	// TODO: Change to ownership
	IterativeMethod *child;

	template <typename V>
	Vector iterateUntilDone(V &&initialIterand, const Linearizer &root,
			LinearSolver &solver, bool solverInitialized = false) {

		// Initialize to the start of the iteration
		initialize();

		// Keep track of the initial iterand
		iterands.push_back( std::forward<V>(initialIterand) );

		// Iterate until done
		while(!done()) {
			for(auto linearizer : linearizers) {
				linearizer->onIterationStart();
			}

			if(child) {
				// Add a new iterand
				iterands.push_back(
					child->iterateUntilDone(
						iterands.back(),
						root,
						solver,
						solverInitialized
					)
				);
			} else {
				// Base case

				// Only compute LHS if necessary
				if(!solverInitialized || root.doesAChange()) {
					solver.initialize( root.A() );
				}

				iterands.push_back(
					solver.solve(
						root.b(),
						iterands.back()
					)
				);
			}

			// Remove the oldest iterand
			// It is important that this happens only after we
			// execute the child's iterateUntilDone routine!
			if(iterands.size() >= lookback) {
				iterands.pop_front();
			}

			// Must be initialized after the first iteration
			solverInitialized = true;
		}

		return iterands.back();

	}

	virtual void initialize() = 0;

public:

	IterativeMethod(unsigned lookback) : lookback(lookback) {
		assert(lookback > 0);
	}

	virtual ~IterativeMethod() {
	}

	// Disable copy constructor and assignment operator.
	IterativeMethod(const IterativeMethod &) = delete;
	IterativeMethod &operator=(const IterativeMethod &) & = delete;

	template <typename V>
	Vector iterateUntilDone(V &&initialIterand, const Linearizer &root,
			LinearSolver &solver) {
		iterateUntilDone( std::forward<V>( initialIterand ), root,
				solver, false );
	}

	virtual bool done() const = 0;

	friend class Linearizer;

};

const std::deque<Vector> &Linearizer::iterands() const {
	return method->iterands;
}

template <typename M>
Linearizer::Linearizer(M &method) : method(&method) {
	method.linearizers.push_back(this);
}

}

#endif

