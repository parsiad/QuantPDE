#ifndef QUANT_PDE_NONLINEAR_ITERATION
#define QUANT_PDE_NONLINEAR_ITERATION

#include <list>    // std::list
#include <memory>  // std::unique_ptr

namespace QuantPDE {

/**
 * Used to solve linear systems.
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
	 * @param x0 An initial guess (ignored for noniterative methods).
	 * @return The solution.
	 * @see QuantPDE::LinearSolver::initialize
	 */
	virtual Vector solve(const Vector &b, const Vector &x0) const = 0;

};

////////////////////////////////////////////////////////////////////////////////

enum class LinearizerStatus { WORKING, PASSTHROUGH, DONE };

class Linearizer {

public:

	virtual LinearizerStatus status() const = 0;

	virtual void beforeBody() = 0;

	virtual void afterBody() = 0;

	virtual bool doesLhsChange() const = 0;

	virtual Matrix lhs() const = 0;

	virtual Vector rhs() const = 0;

};

////////////////////////////////////////////////////////////////////////////////

class IterativeMethod final {

	Vector iterand;
	std::unique_ptr<IterativeMethod> child;
	std::list< std::unique_ptr<Linearizer> > linearizers;

	// TODO: Change std::list to singly-linked list
	template <typename T>
	void run(T &&linearSolver,
			std::list<Linearizer *> &&pass,
			bool initialize = false) {

		std::list<Linearizer *> active;
		for(auto it = linearizers.cbegin(); it != linearizers.cend();
				it++) {
			active.push_back( it->get() );
		}

		while(true) {
			// Run beforeBody() methods and maintain active list
			bool done = true;
			for(auto it = active.cbegin(); it != active.cend();
					it++) {
				switch( (*it)->status() ) {
					case LinearizerStatus::WORKING:
						done = false;
					case LinearizerStatus::PASSTHROUGH:
						(*it)->beforeBody();
						break;
					case LinearizerStatus::DONE:
						active.erase(it);
						break;
				}
			}

			// If the list is empty, we are done this level of
			// iteration
			if(done) {
				return;
			}

			if(child.get() == nullptr) {
				// Base case

				// Create right-hand-side vector and check to
				// see if the left-hand-side will change w.r.t.
				// previous iteration
				Vector b;
				for(auto linearizer : active) {
					initialize = initialize || linearizer
							->doesLhsChange();
					b += linearizer->rhs();
				}

				// Only compute A if necessary
				// TODO: Force compute if has never been
				//       initialized
				if(initialize) {
					Matrix A;
					for(auto linearizer : active) {
						A += linearizer->lhs();
					}

					linearSolver.compute(A);
				}

				iterand = linearSolver.solve(b, iterand);
			} else {
				std::list<Linearizer *> tmp(pass);
				tmp.insert(tmp.cend(), active.cbegin(),
						active.cend());

				child->run(
					std::forward(linearSolver),
					std::move(tmp)
				);
			}

			// Run afterBody() methods (active list is guaranteed
			// not to change here)
			for(auto it = active.cbegin(); it != active.cend();
					it++) {
				(*it)->afterBody();
			}
		}
	}

public:

	template <typename T>
	void run(T &&linearSolver) {
		run(std::forward(linearSolver),
				std::list<Linearizer *>());
	}

};

}

#endif

