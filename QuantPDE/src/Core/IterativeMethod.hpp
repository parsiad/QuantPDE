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

		//#ifndef NDEBUG
		//std::cout << index << " " << size << std::endl;
		//#endif

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

	class Base; typedef std::unique_ptr<Base> B;
	typedef std::unique_ptr<Interpolant<Dimension>> I;

	class Base {

	public:

		Base() noexcept {
		}

		virtual ~Base() {
		}

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

		virtual B clone() const = 0;

	};

	class Constant final : public Base {

		Real constant;

	public:

		Constant(Real constant) noexcept : constant(constant) {
		}

		Constant(const Constant &that) noexcept
				: constant(that.constant) {
		}

		Constant &operator=(const Constant &that) & noexcept {
			constant = that.constant;
		}

		virtual Real value(std::array<Real, Dimension + 1> coordinates)
				const {
			return constant;
		}

		virtual bool isConstantInTime() const {
			return true;
		}

		virtual B clone() const {
			return B(new Constant(*this));
		}

	};

	class FunctionST final : public Base {

		Function<Dimension + 1> function;

	public:

		template <typename F>
		FunctionST(F &&function) noexcept : function(function) {
		}

		FunctionST(const FunctionST &that) noexcept
				: function(that.function) {
		}

		FunctionST(FunctionST &&that) noexcept
				: function( std::move(that.function) ) {
		}

		FunctionST &operator=(const FunctionST &that) & noexcept {
			function = that.function;
		}

		FunctionST &operator=(FunctionST &&that) & noexcept {
			function = std::move(that.function);
		}

		virtual Real value(std::array<Real, Dimension + 1> coordinates)
				const {
			return packAndCall<Dimension + 1>(function,
					coordinates.data());
		}

		virtual B clone() const {
			return B(new FunctionST(*this));
		}

	};

	class FunctionS final : public Base {

		Function<Dimension> function;

	public:

		template <typename F>
		FunctionS(F &&function) noexcept : function(function) {
		}

		FunctionS(const FunctionS &that) noexcept
				: function(that.function) {
		}

		FunctionS(FunctionS &&that) noexcept
				: function( std::move(that.function) ) {
		}

		FunctionS &operator=(const FunctionS &that) & noexcept {
			function = that.function;
		}

		FunctionS &operator=(FunctionS &&that) & noexcept {
			function = std::move(that.function);
		}

		virtual Real value(std::array<Real, Dimension + 1> coordinates)
				const {
			return packAndCall<Dimension>(function,
					coordinates.data() + 1);
		}

		virtual bool isConstantInTime() const {
			return true;
		}

		virtual B clone() const {
			return B(new FunctionS(*this));
		}

	};

	B base;

public:

	class Control final : public Base {

		const InterpolantFactory<Dimension> *factory;
		I interpolant;

		template <typename F>
		Control(F &factory, I interpolant) noexcept
				: factory(&factory),
				interpolant(std::move(interpolant)) {
		}

	public:

		template <typename F>
		Control(F &factory) noexcept
				: factory(&factory),
				interpolant(nullptr) {
		}

		Control(const Control &that) noexcept : factory(that.factory),
				interpolant( that.interpolant->clone() ) {
		}

		Control(Control &&that) noexcept : factory(that.factory),
				interpolant( std::move(that.interpolant) ) {
		}

		Control &operator=(const Control &that) & noexcept {
			factory = that.factory;
			interpolant = that.interpolant->clone();
		}

		Control &operator=(Control &&that) & noexcept {
			factory = that.factory;
			interpolant = std::move(that.interpolant);
		}

		virtual Real value(std::array<Real, Dimension + 1> coordinates)
				const {
			assert(interpolant);

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
			interpolant = factory->make(input);
		}

		virtual void setInput(Vector &&input) {
			interpolant = factory->make( std::move(input) );
		}

		virtual B clone() const {
			return B(new Control(*this));
		}

	};

	/**
	 * Constructor for a constant.
	 */
	WrapperFunction(Real constant) noexcept {
		base = B(new Constant(constant));
	}

	/**
	 * Constructor for a function of space and time.
	 */
	WrapperFunction(const Function<Dimension + 1> &function) noexcept {
		base = B(new FunctionST(function));
	}

	/**
	 * Move constructor for a function of space and time.
	 */
	WrapperFunction(Function<Dimension + 1> &&function) noexcept {
		base = B(new FunctionST(std::move(function)));
	}

	/**
	 * Constructor for a function of space.
	 */
	WrapperFunction(const Function<Dimension> &function) noexcept {
		base = B(new FunctionS(function));
	}

	/**
	 * Move constructor for a function of space.
	 */
	WrapperFunction(Function<Dimension> &&function) noexcept {
		base = B(new FunctionS(std::move(function)));
	}

	/**
	 * Move constructor for a control.
	 */
	WrapperFunction(Control &&control) noexcept {
		base = B(new Control( std::move(control) ));
	}
	WrapperFunction(const Control &control) = delete;

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

};

typedef WrapperFunction<1> WrapperFunction1;
typedef WrapperFunction<2> WrapperFunction2;
typedef WrapperFunction<3> WrapperFunction3;

// Expose Control
template <Index Dimension>
using Control = typename WrapperFunction<Dimension>::Control;

typedef Control<1> Control1;
typedef Control<2> Control2;
typedef Control<3> Control3;

////////////////////////////////////////////////////////////////////////////////

/**
 * A controllable linear system.
 */
class ControlledLinearSystemBase : public LinearSystem {

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
	ControlledLinearSystemBase() noexcept : LinearSystem() {
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
class ControlledLinearSystem : public ControlledLinearSystemBase {

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
	ControlledLinearSystem() noexcept : ControlledLinearSystemBase() {
	}

};

typedef ControlledLinearSystem<1> ControlledLinearSystem1;
typedef ControlledLinearSystem<2> ControlledLinearSystem2;
typedef ControlledLinearSystem<3> ControlledLinearSystem3;

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

	void _clear() {
		for(auto system : systems) {
			system->clear();
		}
		history.clear();
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
	 * Called at the beginning of each iteration.
	 * @param iterand The iterand.
	 * @return True if and only if the iterand was changed.
	 */
	virtual bool transformIterand(Vector &iterand) {
		return false;
	}

	/**
	 * @return True if and only if this iteration is done.
	 */
	virtual bool done() const = 0;

	Vector iterateUntilDone(
		Vector iterand,
		LinearSystemIteration &root,
		LinearSolver &solver,
		Real time,
		bool initialized
	) {
		its.push_back(0);

		clear();
		_clear();

		time = initialTime(time);

		// Iterate until done
		implicitTime = time;
		history.push( std::make_tuple(implicitTime, iterand) );

		// Use macros to prevent branching inside iteration loop

		// Important that time is advanced before onIterationStart()
		// calls
		#define QUANT_PDE_TMP_HEAD                                     \
				do {                                           \
					if(transformIterand(iterand)) {        \
						initialized = false;           \
						_clear();                      \
						history.push( std::make_tuple( \
							implicitTime,          \
							iterand                \
						) );                           \
					}                                      \
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

				iterand = child->iterateUntilDone(
					iterand,
					root,
					solver,
					this->time(0),
					initialized
				);

				history.push( std::make_tuple(implicitTime,
						iterand) );

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

				iterand = solver.solve(
					root.b(),
					iterand
				);

				history.push( std::make_tuple(implicitTime,
						iterand) );

				QUANT_PDE_TMP_TAIL;
			} while( !done() );

		}

		#undef QUANT_PDE_TMP_HEAD
		#undef QUANT_PDE_TMP_TAIL

		// Events occuring at the end
		transformIterand(iterand);

		return iterand;
	}

	template <Index Dimension>
	class InterpolantWrapper {

		typedef std::unique_ptr<Interpolant<Dimension>> I;
		I i;

	public:

		InterpolantWrapper(I i) noexcept : i(std::move(i)) {
		}

		InterpolantWrapper(const InterpolantWrapper &that) noexcept
				: i(that.i->clone()) {
		}

		InterpolantWrapper &operator=(const InterpolantWrapper &that)
				noexcept {
			i = that.i->clone();
			return *this;
		}

		InterpolantWrapper(InterpolantWrapper &&that) noexcept
				: i( std::move(that.i) ) {
		}

		InterpolantWrapper &operator=(InterpolantWrapper &&that)
				noexcept {
			i = std::move(that.i);
			return *this;
		}

		template <typename ...Ts>
		Real operator()(Ts ...coordinates) {
			return (*i)(coordinates...);
		}

	};

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
	 * @param initialCondition The initial condition (as a lambda function).
	 * @param map A mapping from a function to a domain.
	 * @param root The linear system used to generate the left and right
	 *             hand sides of the linear equation to be solved at each
	 *             iteration.
	 * @param solver A linear solver.
	 * @return The solution.
	 * @see QuantPDE::Map
	 */
	template <typename F, Index Dimension>
	InterpolantWrapper<Dimension> solve(
		const Map<Dimension> &map,
		const InterpolantFactory<Dimension> &factory,
		F &&initialCondition,
		LinearSystemIteration &root,
		LinearSolver &solver
	) {
		clearIterations();

		return InterpolantWrapper<Dimension>(
			factory.make(
				iterateUntilDone(
					map(std::forward<F>(
							initialCondition)),
					root,
					solver,
					0.,
					false
				)
			)
		);
	}

	/**
	 * @param initialCondition The initial condition (as a lambda function).
	 * @param domain The spatial domain this problem is defined on.
	 * @param root The linear system used to generate the left and right
	 *             hand sides of the linear equation to be solved at each
	 *             iteration.
	 * @param solver A linear solver.
	 * @return The solution.
	 * @see QuantPDE::Map
	 */
	template <typename F, Index Dimension>
	InterpolantWrapper<Dimension> solve(
		const Domain<Dimension> &domain,
		F &&initialCondition,
		LinearSystemIteration &root,
		LinearSolver &solver
	) {
		PointwiseMap<Dimension> map(domain);
		auto factory = domain.defaultInterpolantFactory();

		return solve(
			map,
			*factory,
			std::forward<F>(initialCondition),
			root,
			solver
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

	// TODO: Remove friendship
	template <Index, bool> friend class EventIterationBase;
	template <Index, bool> friend class EventIteration;
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

