#ifndef QUANT_PDE_CORE_ITERATIVE_METHOD
#define QUANT_PDE_CORE_ITERATIVE_METHOD

// TODO: Remove <Lookback> parameter

#include <deque>        // std::deque
#include <forward_list> // std::list
#include <functional>   // std::function
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

	const T &operator[](int index) const {
		// Lookback is size_t; need to convert to something signed
		// to make the following comparisons
		assert(index >  -(int) Lookback);
		assert(index <=  (int) Lookback);

		size_t position = (tail - 1 + Lookback + index) % Lookback;

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

protected:

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
	virtual Matrix A() = 0;

	/**
	 * @return The right-hand-side vector (b).
	 */
	virtual Vector b() = 0;

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

	template <size_t> friend class Iteration;

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
	Iteration<Lookback> *iteration;

protected:

	Linearizer() noexcept : LinearizerBase(), iteration(nullptr) {
	}

	/**
	 * @return The most recent iterands.
	 */
	const IterandHistory<Lookback> &iterands() const;

	/**
	 * @return The time with which the next solution is associated with.
	 */
	Real nextTime() const;

public:

	// TODO: Document
	void setIteration(Iteration<Lookback> &iteration);

};

/**
 * Executes an iterative method. This class should notbe extended directly.
 * @see QuantPDE::Iteration
 */
class IterationBase {

	IterationBase *child;

	/**
	 * Method called before iteration begins.
	 */
	virtual void clear() {
		// Default: do nothing
	}

	/**
	 * @return A step taken in time (positive or negative). Zero if none is
	 *         taken.
	 */
	virtual Real step() {
		return 0.;
	}

	/**
	 * @return The initial time associated with this iteration.
	 */
	virtual Real initialTime(Real time) const {
		return time;
	}

	/**
	 * @return True if and only if this iteration is done.
	 */
	virtual bool done() const = 0;

	// TODO: Document
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

	virtual ~IterationBase() {
	}

	// TODO: Document
	void setChildIteration(IterationBase &childIteration) {
		child = &childIteration;
	}

	// Disable copy constructor and assignment operator.
	IterationBase(const IterationBase &) = delete;
	IterationBase &operator=(const IterationBase &) = delete;

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
	Real implicitTime;

	size_t its;

	virtual Vector iterateUntilDone(
		const Vector &initialIterand,
		LinearizerBase &root,
		LinearSolver &solver,
		Real time,
		bool initialized
	) {

		its = 0;
		clear();

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
		implicitTime = time;
		do {
			// Iterating at least once allows us to make assumptions
			// about how many iterands we have available to us for
			// processing when implementing the done() method (at
			// least two)

			// Importing that this occurs before onIterationStart()
			// calls
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
					solver.initialize(root.A());
				}

				history.push( std::make_tuple(
					implicitTime,
					solver.solve(
						root.b(),
						std::get<1>( history[0] )
					)
				) );
			}

			// Must be initialized after the first iteration
			initialized = true;

			its++;

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
	Vector iterateUntilDone(
		const Vector &initialIterand,
		LinearizerBase &root,
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

	/**
	 * @return Number of iterations performed.
	 */
	size_t iterations() const {
		return its;
	}

	template <size_t> friend class Linearizer;

};

template <size_t Lookback>
void Linearizer<Lookback>::setIteration(Iteration<Lookback> &iteration) {
	if(this->iteration) {
		this->iteration->linearizers.remove(this);
	}

	this->iteration = &iteration;
	iteration.linearizers.push_front(this);
}

template <size_t Lookback>
const IterandHistory<Lookback> &Linearizer<Lookback>::iterands() const {
	return iteration->history;
}

template <size_t Lookback>
Real Linearizer<Lookback>::nextTime() const {
	return iteration->implicitTime;
}

template <size_t Lookback>
const IterandHistory<Lookback> &Iteration<Lookback>::iterands() const {
	return history;
}

////////////////////////////////////////////////////////////////////////////////

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

template <size_t Controls>
class ControlledLinearOperator : public LinearOperator {
	// TODO
};

////////////////////////////////////////////////////////////////////////////////

template <size_t Lookback = 1, bool Forward = false>
class ConstantStepper : public Iteration<Lookback> {

	Real startTime, endTime, dt;
	unsigned n, steps;

	virtual void clear() {
		n = 0;
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

public:

	ConstantStepper(
		Real startTime,
		Real endTime,
		unsigned steps
	) noexcept :
		startTime(startTime),
		endTime(endTime),
		dt( (endTime - startTime) / steps ),
		steps(steps)
	{
		assert(startTime >= 0.);
		assert(startTime < endTime);
		assert(steps > 0);
	}

};

////////////////////////////////////////////////////////////////////////////////

template <size_t Lookback = 1, bool Forward = false>
class VariableStepper : public Iteration<Lookback> {

	Real startTime, endTime, dt, target, scale, time;
	Real (VariableStepper<Lookback, Forward>::*stepMethod)();

	Real _step0() {
		stepMethod = &VariableStepper::_step1;
		return Forward ? dt : -dt; // TODO: Optimize
	}

	Real _step1() {
		const Vector
			&v1 = std::get<1>( this->iterands()[ 0] ),
			&v0 = std::get<1>( this->iterands()[-1] )
		;

		dt *= target / ( v1 - v0 ).cwiseAbs().cwiseQuotient(
			( scale * Vector::Ones( v1.size() ) ).cwiseMax(
				v1.cwiseAbs().cwiseMax(
					v0.cwiseAbs()
				)
			)
		).maxCoeff();

		// TODO: Optimize
		if(Forward) {
			if(time + dt > endTime) {
				dt = endTime - time;
				time = endTime;
			} else {
				time += dt;
			}

			return dt;
		} else {
			if(time - dt < startTime) {
				dt = time - startTime;
				time = startTime;
			} else {
				time -= dt;
			}

			return -dt;
		}
	}

	virtual void clear() {
		stepMethod = &VariableStepper::_step0;
		time = Forward ? startTime : endTime; // TODO: Optimize
	}

	virtual Real initialTime(Real) const {
		return endTime;
	}

	virtual Real step() {
		return (this->*stepMethod)();
	}

	virtual bool done() const {
		return time <= startTime;
	}

public:

	VariableStepper(
		Real startTime,
		Real endTime,
		Real dt,
		Real target,
		Real scale = 1
	) noexcept :
		startTime(startTime),
		endTime(endTime),
		dt(dt),
		target(target),
		scale(scale)
	{
	}

};

////////////////////////////////////////////////////////////////////////////////

template <size_t Lookback = 1>
class ToleranceIteration : public Iteration<Lookback> {

	Real tolerance;
	Real scale;

	virtual bool done() const {
		const Vector
			&v1 = std::get<1>( this->iterands()[ 0] ),
			&v0 = std::get<1>( this->iterands()[-1] )
		;

		// 2014-06-07: Tested this; it works
		return
			( v1 - v0 ).cwiseAbs().cwiseQuotient(
				( scale * Vector::Ones(v1.size()) ).cwiseMax(
					v1.cwiseAbs()
				)
			).maxCoeff() < tolerance;
	}

public:

	ToleranceIteration(Real tolerance = 1e-6, Real scale = 1) noexcept
			: tolerance(tolerance), scale(scale) {
		assert(tolerance > 0);
		assert(scale > 0);
	}

};

////////////////////////////////////////////////////////////////////////////////

/*
template <Index Dimension, size_t Lookback = 1>
class PenaltyMethod : public Linearizer<Lookback> {

	typedef std::function< bool (
		const Domain<Dimension> &, Index,
		Real,
		const Vector &
	)> F;

	const Domain<Dimension> *domain;
	LinearizerBase *left, *right;

	F predicate;
	Real large;
	Matrix P;

public:

	template <typename D, typename F1>
	PenaltyMethod(D &domain, LinearizerBase &left, LinearizerBase &right,
			F1 &&predicate, Real tolerance = 1e-6) noexcept
			: domain(&domain), left(&left), right(&right),
			predicate( std::forward<F1>(predicate) ),
			large( 1. / tolerance ) {
	}

	virtual void onIterationStart() {
		P = Matrix(domain->size(), domain->size());

		for(Index i = 0; i < domain->size(); i++) {

			bool penalize = predicate(
				*domain, i,
				this->nextTime(),
				std::get<1>( this->iterands()[0] )
			);

			if(penalize) {
				P.insert(i, i) = large;
			}
		}

		P.makeCompressed();
	}

	virtual Matrix A() {
		return left->A() + P * right->A();
	}

	virtual Vector b() {
		return left->b() + P * right->b();
	}

};


template <size_t Lookback = 1>
using PenaltyMethod1 = PenaltyMethod<1, Lookback>;

template <size_t Lookback = 1>
using PenaltyMethod2 = PenaltyMethod<2, Lookback>;

template <size_t Lookback = 1>
using PenaltyMethod3 = PenaltyMethod<3, Lookback>;
*/

} // QuantPDE

#endif

