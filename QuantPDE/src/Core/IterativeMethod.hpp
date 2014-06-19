#ifndef QUANT_PDE_CORE_ITERATIVE_METHOD
#define QUANT_PDE_CORE_ITERATIVE_METHOD

#include <array>        // std::array
#include <cstdlib>      // size_t
#include <forward_list> // std::forward_list
#include <memory>       // std::unique_ptr
#include <tuple>        // std::tuple
#include <utility>      // std::forward, std::move
#include <vector>       // std::vector

namespace QuantPDE {

const int DEFAULT_LOOKBACK = 6;

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
	 * should occur only once so that the matrix is factored only once.
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

/**
 * Circular buffer. When pushing occurs at capacity, the oldest element is
 * removed.
 */
template <typename T>
class CircularBuffer {

protected:

	T *data;

	int tail, n;

	#ifndef NDEBUG
	size_t size;
	#endif

public:

	/**
	 * Constructor.
	 * @param lookback How many iterands to keep track of.
	 */
	CircularBuffer(int lookback) noexcept : n(lookback) {
		assert(lookback > 0);
		data = new T[lookback];
		clear();
	}

	/**
	 * Destructor.
	 */
	virtual ~CircularBuffer() {
		delete [] data;
	}

	// Disable copy constructor and assignment operator.
	CircularBuffer(const CircularBuffer &) = delete;
	CircularBuffer &operator=(const CircularBuffer &) = delete;

	/**
	 * Removes everything from the data structure.
	 */
	void clear() {
		tail = 0;

		#ifndef NDEBUG
		size = 0;
		#endif
	}

	/**
	 * Pushes an element into the buffer.
	 * @param element The element.
	 */
	template <typename E>
	void push(E &&element) {
		data[tail] = std::forward<E>(element);
		tail = (tail + 1) % n;

		#ifndef NDEBUG
		if(size < n) {
			size++;
		}
		#endif
	}

	/**
	 * Retrieves an element from the buffer given an index. The index 0
	 * corresponds to the most recently pushed element (1 corresponds to the
	 * element pushed in after that one, etc.).
	 * @param index The index.
	 * @return The element.
	 */
	const T &operator[](int index) const {
		assert(index >= 0);
		assert(index < n);

		int position = (tail - 1 + n - index) % n;

		assert(position < size);

		return data[position];
	}

	/**
	 * @return The maximum number of iterands one can store in this buffer.
	 */
	int lookback() {
		return n;
	}

};

class Iteration;

class LinearSystem {

public:

	/**
	 * Constructor.
	 */
	LinearSystem() noexcept {
	}

	/**
	 * Destructor.
	 */
	virtual ~LinearSystem() {
	}

	// Disable copy constructor and assignment operator.
	LinearSystem(const LinearSystem &) = delete;
	LinearSystem &operator=(const LinearSystem &) = delete;

	/**
	 * @return False if and only if the left-hand-side matrix (A) has
	 *         changed since the previous call to A.
	 * @see QuantPDE::LinearSystem::A
	 */
	virtual bool isATheSame() const {
		// Default: assume A changes
		return false;
	}

	/**
	 * @return The left-hand-side matrix (A).
	 */
	virtual Matrix A(Real time) = 0;

	/**
	 * @return The right-hand-side vector (b).
	 */
	virtual Vector b(Real time) = 0;

};

////////////////////////////////////////////////////////////////////////////////

/**
 * A convenience class that takes as input (during construction) any one of a
 * constant, a function, or an interpolant factory. It can be used to build
 * robust discrete operators.
 *
 * Many Black-Scholes models assume constant coefficients. For example,
 * \code{.cpp}
 * RectilinearGrid1 grid(Axis { 0., 50., 100., 150., 200. } );
 *
 * BlackScholes blackScholes(
 * 	grid,
 * 	0.04, // Interest rate
 * 	0.2,  // Volatility
 * 	0.,   // Dividend rate
 * );
 * \endcode
 *
 * However, more exotic models are often used. For example,
 * \code{.cpp}
 * const Real alpha = 0.4;
 *
 * BlackScholes blackScholes(
 * 	grid,
 *
 * 	// Controllable interest rate
 * 	Control1::make(grid),
 *
 * 	// Local volatility
 * 	[alpha] (Real t, Real S) { return alpha * t/S; },
 *
 * 	// Constant dividend rate
 * 	0.   // Dividend rate
 * );
 * \endcode
 * The above operator has a controllable interest rate and a local volatility
 * model.
 *
 * This flexibility is made possible by the WrapperFunction class, a wrapper for
 * constants, functions of space and time, functions of space, and controls.
 *
 * @tparam Dimension The dimension of the associated spatial domain (not the
 *                   control domain)
 */
template <Index Dimension>
class WrapperFunction final {

	static_assert(Dimension > 0, "Dimension must be positive");

	class Base {

	public:

		Base() noexcept {
		}

		virtual ~Base() {
		}

		// Disable copy constructor and assignment operator.
		Base(const Base &) = delete;
		Base &operator=(const Base &) = delete;

		virtual Real value(std::array<Real, Dimension + 1> coordinates)
				const = 0;

		virtual bool isConstantInTime() const {
			// Default: not constant w.r.t. time
			return false;
		}

		virtual bool isControllable() const {
			// Default: not controllable
			return false;
		}

		virtual void setInput(const Vector &input) {
			// Default: do nothing
		}

		virtual void setInput(Vector &&input) {
			// Default: do nothing
		}

		virtual std::unique_ptr<Base> clone() const = 0;

	};

	class Constant final : public Base {

		Real constant;

	public:

		Constant(Real constant) noexcept : Base(), constant(constant) {
		}

		virtual Real value(std::array<Real, Dimension + 1> coordinates)
				const {
			return constant;
		}

		virtual bool isConstantInTime() const {
			return true;
		}

		virtual std::unique_ptr<Base> clone() const {
			return std::unique_ptr<Base>(new Constant(constant));
		}

	};

	class FunctionOfSpaceAndTime final : public Base {

		Function<Dimension + 1> function;

	public:

		template <typename F>
		FunctionOfSpaceAndTime(F &&function) noexcept : Base(),
				function(function) {
		}

		virtual Real value(std::array<Real, Dimension + 1> coordinates)
				const {
			return packAndCall<Dimension + 1>(function,
					coordinates.data());
		}

		virtual std::unique_ptr<Base> clone() const {
			return std::unique_ptr<Base>(
					new FunctionOfSpaceAndTime(function));
		}

	};

	class FunctionOfSpace final : public Base {

		Function<Dimension> function;

	public:

		template <typename F>
		FunctionOfSpace(F &&function) noexcept : Base(),
				function(function) {
		}

		virtual Real value(std::array<Real, Dimension + 1> coordinates)
				const {
			return packAndCall<Dimension>(function,
					coordinates.data() + 1);
		}

		virtual bool isConstantInTime() const {
			return true;
		}

		virtual std::unique_ptr<Base> clone() const {
			return std::unique_ptr<Base>(
					new FunctionOfSpace(function));
		}

	};

	class Control final : public Base {

		typedef std::unique_ptr<InterpolantFactoryBase<Dimension>> F;
		typedef std::unique_ptr<Interpolant<Dimension>> I;

		F factory;
		I interpolant;

	public:

		Control(F factory, I interpolant = nullptr) noexcept
				: factory( std::move(factory) ),
				interpolant( std::move(interpolant) ) {
		}

		virtual Real value(std::array<Real, Dimension + 1> coordinates)
				const {
			assert(interpolant != nullptr);

			// Pack and call
			NaryMethodConst<Real, Interpolant<Dimension>,
					Dimension, Real> tmp =
					&Interpolant<Dimension>::operator();
			return packAndCall<Dimension>(*interpolant, tmp,
					coordinates.data() + 1);
		}

		virtual bool isControllable() const {
			return true;
		}

		virtual void setInput(const Vector &input) {
			interpolant = factory->interpolant(input);
		}

		virtual void setInput(Vector &&input) {
			interpolant = factory->interpolant( std::move(input) );
		}

		virtual std::unique_ptr<Base> clone() const {
			return std::unique_ptr<Base>(
				new Control(
					factory->clone(),
					interpolant->clone()
				)
			);
		}

	};

	std::unique_ptr<Base> base;

public:

	/**
	 * Constructor for a constant.
	 */
	WrapperFunction(Real constant) noexcept {
		base = std::unique_ptr<Base>(new Constant(constant));
	}

	/**
	 * Constructor for a function of space and time.
	 */
	WrapperFunction(const Function<Dimension + 1> &function) noexcept {
		base = std::unique_ptr<Base>(
			new FunctionOfSpaceAndTime(
				function
			)
		);
	}

	/**
	 * Move constructor for a function of space and time.
	 */
	WrapperFunction(Function<Dimension + 1> &&function) noexcept {
		base = std::unique_ptr<Base>(
			new FunctionOfSpaceAndTime(
				std::move(function)
			)
		);
	}

	/**
	 * Constructor for a function of space.
	 */
	WrapperFunction(const Function<Dimension> &function) noexcept {
		base = std::unique_ptr<Base>(
			new FunctionOfSpace(
				function
			)
		);
	}

	/**
	 * Move constructor for a function of space.
	 */
	WrapperFunction(Function<Dimension> &&function) noexcept {
		base = std::unique_ptr<Base>(
			new FunctionOfSpace(
				std::move(function)
			)
		);
	}

	/**
	 * Constructor for a control.
	 * @param factory The interpolant factory used to create interpolants
	 *                for this control on the spatial domain.
	 * @see QuantPDE::InterpolantFactoryBase
	 */
	WrapperFunction(std::unique_ptr<InterpolantFactoryBase<Dimension>>
			factory) noexcept {
		base = std::unique_ptr<Base>(
			new Control(
				std::move(factory)
			)
		);
	}

	/**
	 * Copy constructor.
	 */
	WrapperFunction(const WrapperFunction &that) noexcept
			: base(that.base->clone()) {
	}

	/**
	 * Move constructor.
	 */
	WrapperFunction(WrapperFunction &&that) noexcept
			: base( std::move(that.base) ) {
	}

	/**
	 * Copy assignment operator.
	 */
	WrapperFunction &operator=(const WrapperFunction &that) & noexcept {
		base = that.base->clone();
	}

	/**
	 * Move assignment operator.
	 */
	WrapperFunction &operator=(WrapperFunction &&that) & noexcept {
		base = std::move(that.base);
	}

	/**
	 * Query the value of this function at the specified time and
	 * (spatial) coordinates.
	 * @param time The time.
	 * @param coordinates The (spatial) coordinates.
	 * @return The function's value.
	 */
	template <typename ...Ts>
	Real operator()(Real time, Ts ...coordinates) const {
		return base->value( {{time, coordinates...}} );
	}

	/**
	 * @return True if and only if this is not a function of time.
	 */
	bool isConstantInTime() const {
		return base->isConstantInTime();
	}

	/**
	 * @return True if and only if this is a (wrapper for a) control.
	 */
	bool isControllable() const {
		return base->isControllable();
	}

	/**
	 * Sets the value of the control. If this is not a control, nothing is
	 * done.
	 * @param input The value to take on.
	 */
	void setInput(const Vector &input) {
		base->setInput(input);
	}

	/**
	 * Creates a control of particular dimension.
	 * @tparam T The interpolant to use.
	 * @param domain The domain to interpolate on.
	 */
	template <template <Index> class T = PiecewiseLinear, typename D>
	static WrapperFunction make(D &domain) {
		return WrapperFunction(
			std::unique_ptr<InterpolantFactoryBase<Dimension>>(
				new typename T<Dimension>::Factory(domain)
			)
		);
	}

};

typedef WrapperFunction<1> WrapperFunction1;
typedef WrapperFunction<2> WrapperFunction2;
typedef WrapperFunction<3> WrapperFunction3;

template <Index Dimension>
using Control = WrapperFunction<Dimension>;

typedef Control<1> Control1;
typedef Control<2> Control2;
typedef Control<3> Control3;

////////////////////////////////////////////////////////////////////////////////

/**
 * A controllable linear system.
 */
class ControlledLinearSystem : public LinearSystem {

	/**
	 * Controls the system.
	 * @param inputs An array of inputs. It can be assumed that the inputs
	 *               will not be used again and hence invoking move
	 *               semantics on each input is safe.
	 */
	virtual void setInputs(Vector *inputs) = 0;

	/**
	 * @return The number of controls.
	 */
	virtual size_t controlDimension() const = 0;

public:

	/**
	 * Constructor.
	 */
	ControlledLinearSystem() noexcept : LinearSystem() {
	}

	/**
	 * Controls the system.
	 * @param inputs The inputs.
	 */
	template <typename ...Ts>
	void setInputs(Ts &&...inputs) {
		assert(controlDimension() == sizeof...(Ts));

		Vector in[] { std::forward<Ts>(inputs)... };
		setInputs(in);
	}

};

/**
 * A controllable linear system using wrappers as the controls.
 * @see QuantPDE::WrapperFunction
 */
template <Index Dimension>
class SimpleControlledLinearSystem : public ControlledLinearSystem {

	std::vector<WrapperFunction<Dimension> *> controls;

	virtual void setInputs(Vector *inputs) {
		for(auto control : controls) {
			control->setInput( std::move(*(inputs++)) );
		}
	}

	virtual size_t controlDimension() const {
		return controls.size();
	}

protected:

	/**
	 * Must be called to register a control. Once registered, a control can
	 * be controlled using setInputs.
	 * @param wrapper The control.
	 * @see QuantPDE::ControlledLinearSystem::setInputs
	 */
	void registerControl(WrapperFunction<Dimension> &wrapper) {
		if(wrapper.isControllable()) {
			controls.push_back(&wrapper);
		}
	}

public:

	/**
	 * Constructor.
	 */
	SimpleControlledLinearSystem() noexcept : ControlledLinearSystem() {
	}

};

typedef SimpleControlledLinearSystem<1> SimpleControlledLinearSystem1;
typedef SimpleControlledLinearSystem<2> SimpleControlledLinearSystem2;
typedef SimpleControlledLinearSystem<3> SimpleControlledLinearSystem3;

////////////////////////////////////////////////////////////////////////////////

/**
 * Used to generate the left and right-hand sides of the linear system at each
 * iteration.
 */
class LinearSystemIteration : public LinearSystem {

	Iteration *iteration;

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
	 * Method called at the end of an iteration.
	 */
	virtual void onIterationEnd() {
		// Default: do nothing
	}

	/**
	 * @return The left-hand-side matrix (A).
	 */
	virtual Matrix A() = 0;

	/**
	 * @return The right-hand-side vector (b).
	 */
	virtual Vector b() = 0;

	/**
	 * @return The minimum number of previous iterands required to function
	 *         properly.
	 * @see QuantPDE::CircularBuffer
	 */
	virtual int minimumLookback() const = 0;

protected:

	/**
	 * @param index See CircularBuffer for indexing information.
	 * @return Previously encountered time.
	 * @see QuantPDE::CircularBuffer
	 */
	Real time(int index) const;

	/**
	 * @param index See CircularBuffer for indexing information.
	 * @return Previously encountered iterand.
	 * @see QuantPDE::CircularBuffer
	 */
	const Vector &iterand(int index) const;

	/**
	 * @return The time with which the next solution is associated with.
	 */
	Real nextTime() const;

	/**
	 * @return False if and only if the timestep size was different on the
	 *         previous iteration.
	 */
	bool isTimestepTheSame() const {
		Real t2 = nextTime();
		Real t1 = time(0);
		Real t0 = time(1);
		return (t2 - t1) == (t1 - t0);
	}

public:

	/**
	 * Constructor.
	 */
	LinearSystemIteration() noexcept : LinearSystem(), iteration(nullptr) {
	}

	virtual Matrix A(Real) {
		// TODO: Save A
		Matrix _A = A();
		return _A;
	}

	virtual Vector b(Real) {
		// TODO: Save b
		Vector _b = b();
		return _b;
	}

	/**
	 * Associates with this linear system an iterative method.
	 * @param iteration The iterative method.
	 */
	void setIteration(Iteration &iteration);

	friend Iteration;

};

/**
 * Represents an iterative method.
 */
class Iteration {

	Iteration *child;

	// Nonownership
	std::forward_list<LinearSystemIteration *> systems;

	CircularBuffer< std::tuple<Real, Vector> > history;
	Real implicitTime;

	std::vector<size_t> its;

	void clearIterations() {
		its.clear();
		if(child) {
			child->clearIterations();
		}
	}

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
	virtual Real timestep() {
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

	Vector iterateUntilDone(
		const Vector &initialIterand,
		LinearSystemIteration &root,
		LinearSolver &solver,
		Real time,
		bool initialized
	) {
		its.push_back(0);

		clear();

		for(auto system : systems) {
			system->clear();
		}

		time = initialTime(time);

		// Keep track of the initial iterand
		history.clear();
		history.push( std::make_tuple(time, initialIterand) );

		// Iterate until done
		implicitTime = time;

		// Use macros to prevent branching inside iteration loop

		// Important that time is advanced before onIterationStart()
		// calls
		#define QUANT_PDE_TMP_HEAD                                     \
				do {                                           \
					implicitTime += timestep();            \
					for(auto system : systems) {           \
						system->onIterationStart();    \
					}                                      \
				} while(0)

		// Must be initialized after the first iteration
		#define QUANT_PDE_TMP_TAIL                                     \
				do {                                           \
					initialized = true;                    \
					its.back()++;                          \
					for(auto system : systems) {           \
						system->onIterationEnd();      \
					}                                      \
				} while(0)

		// Iterating at least once allows us to make assumptions
		// about how many iterands we have available to us for
		// processing when implementing the done() method

		if(child) {

			// Inductive case (perform inner iteration)

			do {
				QUANT_PDE_TMP_HEAD;

				// Add a new iterand
				history.push( std::make_tuple(
					implicitTime,
					child->iterateUntilDone(
						iterand(0),
						root,
						solver,
						this->time(0),
						initialized
					)
				) );

				QUANT_PDE_TMP_TAIL;
			} while( !done() );

		} else {

			// Base case (solve linear system)

			do {
				QUANT_PDE_TMP_HEAD;

				// Only compute A if necessary
				if(!initialized || !root.isATheSame()) {
					solver.initialize(root.A());
				}

				history.push( std::make_tuple(
					implicitTime,
					solver.solve(
						root.b(),
						iterand(0)
					)
				) );

				QUANT_PDE_TMP_TAIL;
			} while( !done() );

		}

		#undef QUANT_PDE_TMP_HEAD
		#undef QUANT_PDE_TMP_TAIL

		return iterand(0);
	}

protected:

	/**
	 * @param index See CircularBuffer for indexing information.
	 * @return Previously encountered time.
	 * @see QuantPDE::CircularBuffer
	 */
	Real time(int index) const {
		return std::get<0>(history[index]);
	}

	/**
	* @param index See CircularBuffer for indexing information.
	* @return Previously encountered iterand.
	* @see QuantPDE::CircularBuffer
	*/
	const Vector &iterand(int index) const {
		return std::get<1>(history[index]);
	}

public:

	/**
	 * Constructor.
	 */
	Iteration(int lookback) noexcept : child(nullptr), history(lookback) {
		assert(lookback >= 2); // Guarantee at least two iterands are
		                       // saved
	}

	virtual ~Iteration() {
	}

	/**
	 * Associates with this iterative method an inner iterative method,
	 * called on each iteration.
	 * @param innerIteration The inner iterative method.
	 */
	void setInnerIteration(Iteration &innerIteration) {
		child = &innerIteration;
	}

	// Disable copy constructor and assignment operator.
	Iteration(const Iteration &) = delete;
	Iteration &operator=(const Iteration &) = delete;

	/**
	 * Iterates until completion and returns the final iterand.
	 * @param initialIterand The initial iterand.
	 * @param root The linear system used to generate the left and right
	 *             hand sides of the linear equation to be solved at each
	 *             iteration.
	 * @param solver A linear solver.
	 * @return The final iterand.
	 */
	Vector iterateUntilDone(
		const Vector &initialIterand,
		LinearSystemIteration &root,
		LinearSolver &solver
	) {
		clearIterations();

		return iterateUntilDone(
			initialIterand,
			root,
			solver,
			0.,
			false
		);
	}

	/**
	 * @return Vector with number of iterations.
	 */
	const std::vector<size_t> &iterations() const {
		return its;
	}

	/**
	 * @return The average number of iterations.
	 */
	Real meanIterations() const {
		Real mean = its[0];
		for(size_t i = 1; i < its.size(); i++) {
			mean = (i * mean + its[i]) / (i + 1);
		}
		return mean;
	}

	friend LinearSystemIteration;
};

void LinearSystemIteration::setIteration(Iteration &iteration) {
	if(this->iteration) {
		this->iteration->systems.remove(this);
	}

	this->iteration = &iteration;
	iteration.systems.push_front(this);

	assert(minimumLookback() <= this->iteration->history.lookback());
}

Real LinearSystemIteration::time(int index) const {
	return iteration->time(index);
}

const Vector &LinearSystemIteration::iterand(int index) const {
	return iteration->iterand(index);
}

Real LinearSystemIteration::nextTime() const {
	return iteration->implicitTime;
}

////////////////////////////////////////////////////////////////////////////////

/**
 * Iterative method whose stopping condition is
 * \f$ \max_{i}\frac{\left|x_{i}^{k+1}-x_{i}^{k}\right|}{\max\left(\text{scale},\left|x_{i}^{k+1}\right|\right)}<\text{tolerance} \f$
 *
 */
class ToleranceIteration : public Iteration {

	Real tolerance;
	Real scale;

	virtual bool done() const {
		const Vector
			&v1 = iterand(0),
			&v0 = iterand(1)
		;

		// Tested 2014-06-07
		return
			( v1 - v0 ).cwiseAbs().cwiseQuotient(
				( scale * Vector::Ones(v1.size()) ).cwiseMax(
					v1.cwiseAbs()
				)
			).maxCoeff() < tolerance;
	}

public:

	/**
	 * Constructor.
	 * @param tolerance The stopping tolerance.
	 * @param scale The scale parameter.
	 * @param lookback The number of iterands to keep track of.
	 */
	ToleranceIteration(
		Real tolerance = 1e-6,
		Real scale = 1,
		int lookback = DEFAULT_LOOKBACK
	) noexcept :
		Iteration(lookback),
		tolerance(tolerance),
		scale(scale)
	{
		assert(tolerance > 0);
		assert(scale > 0);
	}

};

} // QuantPDE

#endif

