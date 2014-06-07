#ifndef QUANT_PDE_CORE_ITERATIVE_METHOD
#define QUANT_PDE_CORE_ITERATIVE_METHOD

#include <deque>        // std::deque
#include <forward_list> // std::list
#include <memory>       // std::unique_ptr
#include <tuple>        // std::tuple
#include <type_traits>  // std::conditional
#include <utility>      // std::forward, std::move

namespace QuantPDE {

/**
 * A pure virtual class representing a solver for equations of type \f$Ax=b\f$.
 */
class LinearSolver {

protected:

	Matrix A;

	virtual void initialize() = 0;

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
	template <typename M>
	void initialize(M &&A) {
		this->A = std::forward<M>(A);
		initialize();
	}

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

	virtual void initialize() {
		solver.compute(A);
		assert( solver.info() == Eigen::Success );
	}

public:

	/**
	 * Constructor.
	 */
	BiCGSTABSolver() noexcept : LinearSolver() {
	}

	virtual Vector solve(const Vector &b, const Vector &guess) const {
		Vector v = solver.solveWithGuess(b, guess);
		assert( solver.info() == Eigen::Success );
		return v;
	}

};

////////////////////////////////////////////////////////////////////////////////

// Forward declaration
template <size_t> class Iteration;
template <size_t> class Iterands;

// TODO: Document
template <size_t Lookback>
class IterandHistory final {

	static_assert(Lookback > 0, "The number of iterands to store must be "
			"positive");

	typedef std::tuple<Real, Vector> T;

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

	const T &operator[](size_t index) const {
		size_t position = (tail - 1 + index) % Lookback;
		assert(position < size);
		return data[position];
	}

};

/**
 * Used to generate the left and right-hand sides of the linear system at each
 * iteration. This class should not be extended directly.
 * @see QuantPDE::Linearizer
 */
class LinearizerBase {

public:

	/**
	 * Constructor.
	 */
	LinearizerBase() {
	}

	/**
	 * Destructor.
	 */
	virtual ~LinearizerBase() {
	}

	// Disable copy constructor and assignment operator.
	LinearizerBase(const LinearizerBase &) = delete;
	LinearizerBase &operator=(const LinearizerBase &) = delete;

	/**
	 * Method called before iteration begins.
	 */
	virtual void clear() {
		// Default: do nothing
	}

	/**
	 * Method called on the start of an iteration.
	 */
	virtual void onIterationStart() {
		// Default: do nothing
	}

	/**
	 * @return True if and only if the left-hand-side matrix (A) changes
	 *         from iteration to iteration.
	 */
	virtual bool doesAChange() const {
		// Default: assume A changes
		return true;
	}

	/**
	 * @return The left-hand-side matrix (A).
	 */
	virtual Matrix A(Real dt) = 0;

	/**
	 * @return The right-hand-side vector (b).
	 */
	virtual Vector b(Real dt) = 0;

};

/**
 * Used to generate the left and right-hand sides of the linear system at each
 * iteration.
 */
template <size_t Lookback>
class Linearizer : public LinearizerBase {

	static_assert(Lookback > 0, "The number of iterands to store must be "
			"positive");

	// Nonownership
	const Iteration<Lookback> *iteration;

protected:

	/**
	 * @return The most recent iterands.
	 */
	const IterandHistory<Lookback> &iterands() const;

public:

	/**
	 * Constructor.
	 */
	template <typename I>
	Linearizer(I &iteration);

};

/**
 * Executes an iterative method. This class should notbe extended directly.
 * @see QuantPDE::Iteration
 */
class IterationBase {

	IterationBase *child;

	virtual Real step() {
		return 0.;
	}

	virtual Real initialTime(Real time) const {
		return time;
	}

	virtual Vector iterateUntilDone(
		const Vector &initialIterand,
		LinearizerBase &root,
		LinearSolver &solver,
		Real time,
		bool initialized
	) = 0;

public:

	/**
	 * Constructor.
	 */
	IterationBase() noexcept : child(nullptr) {
	}

	// TODO: Document
	template <typename I>
	void setChildIteration(I &childIteration) {
		child = &childIteration;
	}

	virtual ~IterationBase() {
	}

	// Disable copy constructor and assignment operator.
	IterationBase(const IterationBase &) = delete;
	IterationBase &operator=(const IterationBase &) = delete;

	/**
	 * @return True if and only if this iteration has finished.
	 */
	virtual bool done() const = 0;

	template <size_t> friend class Iteration;

};

/**
 * Executes an iterative method.
 */
template <size_t Lookback>
class Iteration : public IterationBase {

	static_assert(Lookback > 0, "The number of iterands to store must be "
			"positive");

	typedef std::tuple<Real, Vector> T;

	// Nonownership
	std::forward_list<Linearizer<Lookback> *> linearizers;

	IterandHistory<Lookback> history;

	virtual Vector iterateUntilDone(
		const Vector &initialIterand,
		LinearizerBase &root,
		LinearSolver &solver,
		Real time,
		bool initialized
	) {

		for(auto linearizer : linearizers) {
			linearizer->clear();
		}

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

			implicitTime += step();

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
						root.b(implicitTime),
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

	/**
	 * @return The most recent iterands.
	 */
	const IterandHistory<Lookback> &iterands() const;

public:

	/**
	 * Constructor.
	 */
	Iteration() noexcept {
	}

	// TODO: Document
	template <typename L>
	Vector iterateUntilDone(
		const Vector &initialIterand,
		L &root,
		LinearSolver &solver
	) {
		return iterateUntilDone(
			initialIterand,
			root,
			solver,
			0.,
			false
		);
	}

	template <size_t> friend class Linearizer;

};

template <size_t Lookback> template <typename I>
Linearizer<Lookback>::Linearizer(I &iteration) {
	this->iteration = &iteration;
	iteration.linearizers.push_front(this);
}

template <size_t Lookback>
const IterandHistory<Lookback> &Linearizer<Lookback>::iterands() const {
	return iteration->history;
}

template <size_t Lookback>
const IterandHistory<Lookback> &Iteration<Lookback>::iterands() const {
	return history;
}

////////////////////////////////////////////////////////////////////////////////

template <size_t Dimension, size_t Controls>
class LinearOperator {

public:

	LinearOperator() noexcept {
	}

	virtual ~LinearOperator() {
	}

	// Disable copy constructor and assignment operator.
	LinearOperator(const LinearOperator &) = delete;
	LinearOperator &operator=(const LinearOperator &) & = delete;

	virtual Matrix discretize(Real time) const = 0;

};

////////////////////////////////////////////////////////////////////////////////

template <size_t Lookback = 1, bool Forward = false>
class ConstantStepper : public Iteration<Lookback> {

	Real startTime, endTime, dt;
	unsigned n, steps;

public:

	ConstantStepper(
		Real startTime,
		Real endTime,
		unsigned steps
	) noexcept :
		startTime(startTime),
		endTime(endTime),
		dt( (endTime - startTime) / steps ),
		n(0),
		steps(steps)
	{
		assert(startTime >= 0.);
		assert(startTime < endTime);
		assert(steps > 0);
	}

	virtual Real initialTime(Real) const {
		return Forward ? startTime : endTime; // TODO: Optimize
	}

	virtual Real step() {
		n++;
		return Forward ? dt : -dt; // TODO: Optimize
	}

	virtual bool done() const {
		return n >= steps;
	}

};

////////////////////////////////////////////////////////////////////////////////

template <size_t Lookback = 1, bool Forward = false>
class VariableStepper : public Iteration<Lookback> {

	Real startTime, endTime, dt, dnorm, scale, time;

public:

	VariableStepper(
		Real startTime,
		Real endTime,
		Real dt,
		Real dnorm,
		Real scale = 1
	) noexcept :
		startTime(startTime),
		endTime(endTime),
		dt(dt),
		dnorm(dnorm),
		scale(scale),
		time(endTime)
	{
	}

	virtual Real initialTime(Real) const {
		return endTime;
	}

	virtual Real step() {
		const Vector
			&v_n0 = std::get<1>( this->iterands()[-1] ),
			&v_n1 = std::get<1>( this->iterands()[ 0] )
		;

		dt *= dnorm / ( v_n1 - v_n0 ).cwiseAbs().cwiseQuotient(
			( scale * Vector::Ones( v_n1.size() ) ).cwiseMax(
				v_n1.cwiseAbs().cwiseMax(
					v_n0.cwiseAbs()
				)
			)
		).maxCoeff;

		// TODO: Optimize
		if(Forward) {
			if(time + dt > endTime) {
				dt = endTime - time;
			}

			return dt;
		} else {
			if(time - dt < startTime) {
				dt = time - startTime;
			}

			return -dt;
		}
	}

	virtual bool done() {
		return time <= startTime;
	}

};

////////////////////////////////////////////////////////////////////////////////

template <size_t Lookback = 1>
class ToleranceIteration : public Iteration<Lookback> {

	Real toleranceSquared;

public:

	ToleranceIteration(
		Real tolerance = 1e-6
	) noexcept :
		toleranceSquared(tolerance * tolerance)
	{
	}

	virtual bool done() const {
		const Vector
			&v_n0 = std::get<1>( this->iterands()[-1] ),
			&v_n1 = std::get<1>( this->iterands()[ 0] )
		;

		return ( v_n1 - v_n0 ).squaredNorm() < toleranceSquared;
	}

};

} // QuantPDE

#endif

