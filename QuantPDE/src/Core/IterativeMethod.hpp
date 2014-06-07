#ifndef QUANT_PDE_ITERATIVE_METHOD
#define QUANT_PDE_ITERATIVE_METHOD

#include <deque>   // std::deque
#include <list>    // std::list
#include <memory>  // std::unique_ptr
#include <tuple>   // std::tuple
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

// Forward declaration
template <size_t> class Iteration;
template <size_t> class Iterands;

/**
 * Used to generate the left and right-hand sides of the linear system at each
 * iteration. This class should not be extended directly.
 * @see QuantPDE::Linearizer
 */
class LinearizerBase {

public:

	virtual ~LinearizerBase() {
	}

	// Disable copy constructor and assignment operator.
	LinearizerBase(const LinearizerBase &) = delete;
	LinearizerBase &operator=(const LinearizerBase &) = delete;

	virtual void onIterationStart() {
		// Default: do nothing
	}

	virtual bool doesAChange() const {
		// Default: assume A changes
		return true;
	}

	virtual Matrix A(Real implicitTime) const = 0;

	virtual Vector b() const = 0;

};

/**
 * Used to generate the left and right-hand sides of the linear system at each
 * iteration.
 */
template <size_t Lookback>
class Linearizer : public LinearizerBase {

	// Nonownership
	const Iteration<Lookback> *iteration;

protected:

	const Iterands<Lookback> &iterands() const;

public:

	template <typename I>
	Linearizer(I &iteration);

};

/**
 * Executes an iterative method. This class should notbe extended directly.
 * @see QuantPDE::Iteration
 */
class IterationBase {

	IterationBase *child;

	virtual Real nextTime(Real time) {
		return time;
	}

	virtual Real initialTime(Real time) const {
		return time;
	}

	virtual Vector iterateUntilDone(
		const Vector &initialIterand,
		const LinearizerBase &root,
		LinearSolver &solver,
		Real time,
		bool initialized
	) = 0;

public:

	template <typename I>
	IterationBase(I &childIteration) noexcept : child(&childIteration) {
	}

	virtual ~IterationBase() {
	}

	// Disable copy constructor and assignment operator.
	IterationBase(const IterationBase &) = delete;
	IterationBase &operator=(const IterationBase &) = delete;

	virtual bool done() const = 0;

	template <size_t> friend class Iteration;

};

/**
 * Executes an iterative method.
 */
template <size_t Lookback>
class Iteration : public IterationBase {

	typedef std::tuple<Real, Vector> T;

	class IterandHistory final {

		T data[Lookback];
		size_t tail;

		#ifndef NDEBUG
		size_t size;
		#endif

	public:

		IterandHistory() noexcept : tail(0) {
			#ifndef NDEBUG
			size = 0;
			#endif
		}

		IterandHistory(const IterandHistory &) = delete;
		IterandHistory &operator=(const IterandHistory &) = delete;

		void clear() {
			tail = 0;

			#ifndef NDEBUG
			size = 0;
			#endif
		}

		void push(T &&element) {
			data[tail] = std::forward<T>(element);
			tail = (tail + 1) % Lookback;

			#ifndef NDEBUG
			if(size < Lookback) {
				size++;
			}
			#endif
		}

		const T &operator[](size_t index) {
			size_t position = (tail - 1 + index) % Lookback;
			assert(position < size);
			return data[position];
		}

	};

	// Nonownership
	std::list<Linearizer<Lookback> *> linearizers;

	IterandHistory history;

	virtual Vector iterateUntilDone(
		const Vector &initialIterand,
		const LinearizerBase &root,
		LinearSolver &solver,
		Real time,
		bool initialized
	) {

		assert(child != nullptr);

		time = initialTime(time);

		// Keep track of the initial iterand
		history.clear();
		history.push( std::make_tuple(
			time,
			initialIterand
		) );

		// Iterate until done
		Real implicitTime = time;
		do {
			// Iterating at least once allows us to make assumptions
			// about how many iterands we have available to us for
			// processing when implementing the done() method

			implicitTime = nextTime(implicitTime);

			for(auto linearizer : linearizers) {
				linearizer->onIterationStart();
			}

			if(child) {
				// Inductive case (perform inner iteration)

				// Add a new iterand
				const T &tuple = history[0];
				history.push( std::make_tuple(
					implicitTime,
					child->iterateUntilDone(
						std::get<1>( tuple ),
						root,
						solver,
						std::get<0>( tuple ),
						initialized
					)
				) );
			} else {
				// Base case (solve linear system)

				// Only compute A if necessary
				if(!initialized || root.doesAChange()) {
					solver.initialize(root.A(implicitTime));
				}

				history.push( std::make_tuple(
					implicitTime,
					solver.solve(
						root.b(),

						// Initial guess
						std::get<1>( history[0] )
					)
				) );
			}

			// Must be initialized after the first iteration
			initialized = true;

		} while( !done() );

		return std::get<1>( history[0] );

	}

protected:

	const Iterands<Lookback> &iterands() const;

public:

	template <typename I>
	Iteration(I &childIteration) noexcept : IterationBase(childIteration) {
	}

	Vector solve(
		const Vector &initialIterand,
		const LinearizerBase &root,
		LinearSolver &solver
	) {
		return solve(
			initialIterand,
			root,
			solver,
			0.,
			false
		);
	}

	template <size_t> friend class Linearizer;
	template <size_t> friend class Iterands;

};

template <size_t Lookback>
class Iterands final {

	typedef typename Iteration<Lookback>::IterandHistory H;
	typedef std::tuple<Real, Vector> T;

	const H *history;

public:

	Iterands(const H &history) noexcept
			: history(&history) {
	}

	Iterands(const Iterands &that) noexcept : history(that.history) {
	}

	Iterands &operator=(const Iterands &that) & noexcept {
		history = that.history;
	}

	const T &operator[](size_t index) {
		return (*history)[index];
	}

};

template <size_t Lookback>
const Iterands<Lookback> &Linearizer<Lookback>::iterands() const {
	return Iterands<Lookback>(iteration->history);
}

template <size_t Lookback>
const Iterands<Lookback> &Iteration<Lookback>::iterands() const {
	return Iterands<Lookback>(history);
}

////////////////////////////////////////////////////////////////////////////////

template <size_t Lookback>
class BackwardsConstantStepper : public Iteration<Lookback> {

	Real endTime, dt;
	unsigned n, steps;

public:

	template <typename I>
	BackwardsConstantStepper(
		I &childIteration,
		Real startTime,
		Real endTime,
		unsigned steps
	) noexcept :
		Iteration<Lookback>(childIteration),
		endTime(endTime),
		dt( (endTime - startTime) / steps ),
		n(0),
		steps(steps)
	{
		assert(startTime >= 0.);
		assert(startTime < endTime);
		assert(steps > 0);
	}

	virtual bool done() const {
		return n >= steps;
	}

	virtual Real nextTime(Real) {
		return endTime - (++n) * dt;
	}

	virtual Real initialTime(Real) const {
		return endTime;
	}

};

////////////////////////////////////////////////////////////////////////////////

template <size_t Lookback>
class ToleranceIteration : public Iteration<Lookback> {

	Real toleranceSquared;

public:

	template <typename I>
	ToleranceIteration(
		I &childIteration,
		Real tolerance = 1e-6
	) noexcept :
		Iteration<Lookback>(childIteration),
		toleranceSquared(tolerance * tolerance)
	{
	}

	virtual bool done() const {
		auto x = this->iterands();
		return ( std::get<1>( x[-1] ) - std::get<1>( x[0] ) )
				.squaredNorm() < toleranceSquared;
	}

};

} // QuantPDE

#endif

