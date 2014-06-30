#ifndef QUANT_PDE_CORE_ITERATIVE_METHOD
#define QUANT_PDE_CORE_ITERATIVE_METHOD

#include <array>        // std::array
#include <cstdlib>      // std:abs, size_t
#include <forward_list> // std::forward_list
#include <memory>       // std::shared_ptr, std::unique_ptr
#include <queue>        // std::priority_queue
#include <tuple>        // std::tuple
#include <utility>      // std::forward, std::move
#include <vector>       // std::vector

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
 * 	Control1(grid),
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
 * This flexibility is made possible by the MultiFunction class, a wrapper for
 * constants, functions of space and time, functions of space, and controls.
 *
 * @tparam Dimension The dimension of the associated spatial domain (not the
 *                   control domain)
 */
template <Index Dimension>
class MultiFunction final {

	static_assert(Dimension > 0, "Dimension must be positive");

	class Base;
	typedef std::unique_ptr<Base> B;

	typedef InterpolantWrapper<Dimension> WI;
	typedef InterpolantFactoryWrapper<Dimension> WF;

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

		WF factory;
		WI interpolant;

	public:

		Control(WF factory) noexcept
				: factory(std::move(factory)),
				interpolant(nullptr) {
		}

		// explicit Control(const RectilinearGrid<Dimension> &grid)
		//
		// Clang seems to get tripped up at trying to disambiguate the
		// (WF factory) constructor from the one above, hence the
		// following workaround with templates:
		template <typename G>
		Control(G &grid) noexcept
				: factory(grid.defaultInterpolantFactory()),
				interpolant(nullptr) {
		}

		Control(const Control &that) noexcept : factory(that.factory),
				interpolant(that.interpolant) {
		}

		Control(Control &&that) noexcept
				: factory(std::move(that.factory)),
				interpolant(std::move(that.interpolant)) {
		}

		Control &operator=(const Control &that) & noexcept {
			factory = that.factory;
			interpolant = that.interpolant;
		}

		Control &operator=(Control &&that) & noexcept {
			factory = std::move(that.factory);
			interpolant = std::move(that.interpolant);
		}

		virtual Real value(std::array<Real, Dimension + 1> coordinates)
				const {
			assert(interpolant);

			// Pack and call
			NaryMethodConst<Real, Interpolant<Dimension>,
					Dimension, Real> tmp =
					&Interpolant<Dimension>::operator();
			return packAndCall<Dimension>(interpolant, tmp,
					coordinates.data() + 1);
		}

		virtual bool isControllable() const {
			return true;
		}

		virtual void setInput(const Vector &input) {
			interpolant = factory.make(input);
		}

		virtual void setInput(Vector &&input) {
			interpolant = factory.make( std::move(input) );
		}

		virtual B clone() const {
			return B(new Control(*this));
		}

	};

	/**
	 * Constructor for a constant.
	 */
	MultiFunction(Real constant) noexcept {
		base = B(new Constant(constant));
	}

	/**
	 * Constructor for a function of space and time.
	 */
	MultiFunction(const Function<Dimension + 1> &function) noexcept {
		base = B(new FunctionST(function));
	}

	/**
	 * Move constructor for a function of space and time.
	 */
	MultiFunction(Function<Dimension + 1> &&function) noexcept {
		base = B(new FunctionST(std::move(function)));
	}

	/**
	 * Constructor for a function of space.
	 */
	MultiFunction(const Function<Dimension> &function) noexcept {
		base = B(new FunctionS(function));
	}

	/**
	 * Move constructor for a function of space.
	 */
	MultiFunction(Function<Dimension> &&function) noexcept {
		base = B(new FunctionS(std::move(function)));
	}

	/**
	 * Move constructor for a control.
	 */
	MultiFunction(Control &&control) noexcept {
		base = B(new Control(std::move(control)));
	}
	MultiFunction(const Control &control) = delete;

	/**
	 * Copy constructor.
	 */
	MultiFunction(const MultiFunction &that) noexcept
			: base(that.base->clone()) {
	}

	/**
	 * Move constructor.
	 */
	MultiFunction(MultiFunction &&that) noexcept
			: base( std::move(that.base) ) {
	}

	/**
	 * Copy assignment operator.
	 */
	MultiFunction &operator=(const MultiFunction &that) & noexcept {
		base = that.base->clone();
	}

	/**
	 * Move assignment operator.
	 */
	MultiFunction &operator=(MultiFunction &&that) & noexcept {
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

typedef MultiFunction<1> MultiFunction1;
typedef MultiFunction<2> MultiFunction2;
typedef MultiFunction<3> MultiFunction3;

// Expose Control
template <Index Dimension>
using Control = typename MultiFunction<Dimension>::Control;

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
 * @see QuantPDE::MultiFunction
 */
template <Index Dimension>
class ControlledLinearSystem : public ControlledLinearSystemBase {

	std::vector<MultiFunction<Dimension> *> controls;

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
	void registerControl(MultiFunction<Dimension> &wrapper) {
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
class IterationNode : public LinearSystem {

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
	 * @return The minimum number of previous iterands required to function
	 *         properly.
	 * @see QuantPDE::CircularBuffer
	 */
	virtual int minimumLookback() const {
		// Default is 1
		return 1;
	}

protected:

	/**
	 * @param index See CircularBuffer for indexing information.
	 * @return Previously encountered time.
	 * @see QuantPDE::CircularBuffer
	 */
	inline Real time(int index) const;

	/**
	 * @param index See CircularBuffer for indexing information.
	 * @return Previously encountered iterand.
	 * @see QuantPDE::CircularBuffer
	 */
	inline const Vector &iterand(int index) const;

	/**
	 * @return The time with which the next solution is associated with.
	 */
	inline Real nextTime() const;

	/**
	 * @return False if and only if the timestep size was different on the
	 *         previous iteration.
	 */
	inline bool isTimestepTheSame() const;

public:

	/**
	 * Constructor.
	 */
	IterationNode() noexcept : LinearSystem(), iteration(nullptr) {
	}

	/**
	 * Associates with this linear system an iterative method.
	 * @param iteration The iterative method.
	 */
	void setIteration(Iteration &iteration);

	friend Iteration;

};

////////////////////////////////////////////////////////////////////////////////

// TODO: Make solution object so that Iteration can seem stateless

/**
 * An iterative method.
 */
class Iteration {

	typedef CircularBuffer<std::tuple<Real, Vector>> CB;

	virtual Vector iterateUntilDone(
		Vector iterand,
		IterationNode &root,
		LinearSolver &solver,
		Real time,
		bool initialized
	) = 0;

	inline void clearNodes() {
		for(auto node : nodes) {
			node->clear();
		}
	}

	inline void startNodes() {
		for(auto node : nodes) {
			node->onIterationStart();
		}
	}

	inline void endNodes() {
		for(auto node : nodes) {
			node->onIterationEnd();
		}
	}

	inline void solveLinearSystemAndSaveResult(
		IterationNode &root,
		LinearSolver &solver,
		bool initialized
	) {
		assert(implicitTime >= 0.);
		if(!initialized || !root.isATheSame()) {
			solver.initialize( root.A(implicitTime) );
		}

		this->history->push(std::make_tuple(
			this->implicitTime,
			solver.solve(
				root.b(implicitTime),
				this->iterand(0)
			)
		));
	}

	virtual bool isTimestepTheSame() const = 0;

	virtual int minimumLookback() const {
		// Default is 1
		return 1;
	}

protected:

	Iteration *child;
	std::forward_list<IterationNode *> nodes;
	CB *history;
	Real implicitTime;
	std::vector<size_t> its;

	/**
	 * @param index See CircularBuffer for indexing information.
	 * @return Previously processed time, in order of most-to-least recent.
	 * @see QuantPDE::CircularBuffer
	 */
	Real time(int index) const {
		return std::get<0>( (*history)[index] );
	}

	/**
	* @param index See CircularBuffer for indexing information.
	* @return Previously processed iterand, in order of most-to-least
	*         recent.
	* @see QuantPDE::CircularBuffer
	*/
	const Vector &iterand(int index) const {
		return std::get<1>( (*history)[index] );
	}

public:

	/**
	 * Constructor.
	 */
	Iteration() noexcept : child(nullptr), history(nullptr) {
	}

	/**
	 * Destructor.
	 */
	virtual ~Iteration() {
		delete history;
	}

	/**
	 * Sets the inner iterative method.
	 * @param innerIteration An iterative method.
	 */
	void setInnerIteration(Iteration &innerIteration) {
		child = &innerIteration;
	}

	// Disable copy constructor and assignment operator.
	Iteration(const Iteration &) = delete;
	Iteration &operator=(const Iteration &) = delete;

	/**
	 * @return The time at which this iterative method is currently
	 *         computing a solution for.
	 */
	Real nextTime() const {
		return implicitTime;
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

	/**
	 * @param map Maps the initial condition to the domain nodes.
	 * @param factory Interpolates the solution on the domain nodes to the
	 *                entire domain. The resulting interpolant is returned.
	 * @param initialCondition The initial condition.
	 * @param root The node used to generate the left and right hand sides
	 *             of the linear equation to be solved at each iteration.
	 * @param solver A linear solver.
	 * @return The solution.
	 * @see QuantPDE::Map
	 */
	template <typename F, Index Dimension>
	InterpolantWrapper<Dimension> solve(
		const Map<Dimension> &map,
		const InterpolantFactory<Dimension> &factory,
		F &&initialCondition,
		IterationNode &root,
		LinearSolver &solver
	) {
		// Clear iteration count
		Iteration *current = this;
		do {
			// Initialize history
			int lookback = current->minimumLookback();
			for(auto node : current->nodes) {
				int tmp = node->minimumLookback();
				if(tmp > lookback) {
					lookback = tmp;
				}
				delete current->history;
				current->history = new CB(lookback);
			}

			// Clear number of iterations
			current->its.clear();

			// Go to the next iterative method
			current = current->child;
		} while(current);

		// Map, iterate, interpolate
		return InterpolantWrapper<Dimension>(
			factory.make(
				iterateUntilDone(
					map(std::forward<F>(initialCondition)),
					root,
					solver,
					-1., // Use a time value that is bogus
					false
				)
			)
		);
	}

	/**
	 * @param domain The spatial domain this problem is defined on.
	 * @param initialCondition The initial condition.
	 * @param root The node used to generate the left and right hand sides
	 *             of the linear equation to be solved at each iteration.
	 * @param solver A linear solver.
	 * @return The solution.
	 * @see QuantPDE::Map
	 */
	template <typename F, Index Dimension>
	InterpolantWrapper<Dimension> solve(
		const Domain<Dimension> &domain,
		F &&initialCondition,
		IterationNode &root,
		LinearSolver &solver
	) {
		return solve(
			PointwiseMap<Dimension>(domain),
			domain.defaultInterpolantFactory(),
			std::forward<F>(initialCondition),
			root,
			solver
		);
	}

	friend IterationNode;
	friend class ToleranceIteration;
	template <bool> friend class TimeIteration;

};

////////////////////////////////////////////////////////////////////////////////

void IterationNode::setIteration(Iteration &iteration) {
	if(this->iteration) {
		this->iteration->nodes.remove(this);
	}
	this->iteration = &iteration;

	iteration.nodes.push_front(this);
}

Real IterationNode::time(int index) const {
	return iteration->time(index);
}

const Vector &IterationNode::iterand(int index) const {
	return iteration->iterand(index);
}

Real IterationNode::nextTime() const {
	return iteration->implicitTime;
}

bool IterationNode::isTimestepTheSame() const {
	return iteration->isTimestepTheSame();
}

////////////////////////////////////////////////////////////////////////////////

#define QUANT_PDE_TMP_OUTER_HEAD
#define QUANT_PDE_TMP_OUTER_TAIL
#define QUANT_PDE_TMP_TIMESTEP

#define QUANT_PDE_TMP_SET_TIME \
	do { \
		this->implicitTime = time; \
	} while(0)

#define QUANT_PDE_TMP_HEAD \
	do { \
		this->startNodes(); \
	} while(0)

#define QUANT_PDE_TMP_TAIL \
	do { \
		initialized = true; \
		this->its.back()++; \
		this->endNodes(); \
	} while(0)

#define QUANT_PDE_TMP_NOT_DONE \
	( ( this->iterand(0) - this->iterand(1) ).cwiseAbs().cwiseQuotient( \
		( scale * Vector::Ones(this->iterand(0).size()) ).cwiseMax( \
			this->iterand(0).cwiseAbs() \
		) \
	).maxCoeff() > tolerance )

#define QUANT_PDE_TMP_ITERATE_UNTIL_DONE \
	do { \
		QUANT_PDE_TMP_SET_TIME; \
		this->its.push_back(0); \
		this->clearNodes(); \
		this->history->clear(); \
		this->history->push(std::make_tuple( \
			this->implicitTime, \
			std::move(iterand) \
		)); \
		if(this->child) { \
			QUANT_PDE_TMP_OUTER_HEAD; \
			do { \
				QUANT_PDE_TMP_TIMESTEP; \
				QUANT_PDE_TMP_HEAD; \
				this->history->push(std::make_tuple( \
					this->implicitTime, \
					this->child->iterateUntilDone( \
						this->iterand(0), \
						root, \
						solver, \
						this->implicitTime, \
						initialized \
					) \
				)); \
				QUANT_PDE_TMP_TAIL; \
			} while(QUANT_PDE_TMP_NOT_DONE); \
			QUANT_PDE_TMP_OUTER_TAIL; \
		} else { \
			QUANT_PDE_TMP_OUTER_HEAD; \
			do { \
				QUANT_PDE_TMP_TIMESTEP; \
				QUANT_PDE_TMP_HEAD; \
				this->solveLinearSystemAndSaveResult(root, \
						solver, initialized); \
				QUANT_PDE_TMP_TAIL; \
			} while(QUANT_PDE_TMP_NOT_DONE); \
			QUANT_PDE_TMP_OUTER_TAIL; \
		} \
		return this->iterand(0); \
	} while(0)

/**
 * An iterative method that terminates when adjacent iterands are within a
 * certain error tolerance. Specifically, the stopping condition is
 * \f$ \max_{i}\frac{\left|x_{i}^{k+1}-x_{i}^{k}\right|}{\max\left(\text{scale},\left|x_{i}^{k+1}\right|\right)}<\text{tolerance} \f$
 */
class ToleranceIteration final : public Iteration {

	virtual Vector iterateUntilDone(
		Vector iterand,
		IterationNode &root,
		LinearSolver &solver,
		Real time,
		bool initialized
	) {
		QUANT_PDE_TMP_ITERATE_UNTIL_DONE;
	}

	virtual bool isTimestepTheSame() const {
		return true;
	}

	Real tolerance, scale;

	virtual int minimumLookback() const {
		return 2;
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
		Real scale = 1
	) noexcept :
		tolerance(tolerance),
		scale(scale)
	{
		assert(tolerance > 0);
		assert(scale > 0);
	}

};

////////////////////////////////////////////////////////////////////////////////

#undef QUANT_PDE_TMP_SET_TIME
#define QUANT_PDE_TMP_SET_TIME \
	do { \
		this->implicitTime = initialTime(); \
		dt = -1.; \
	} while(0)

// Copy the priority queue so this iteration can be used again in the future
// and add the NullEvent so as to make sure that iteration terminates
#undef QUANT_PDE_TMP_OUTER_HEAD
#define QUANT_PDE_TMP_OUTER_HEAD \
	PriorityQueue events(this->events); \
	events.emplace( \
		std::numeric_limits<unsigned>::max(), \
		terminalTime(), \
		std::shared_ptr<EventBase>(new NullEvent{}) \
	); \
	do { \
		const Real nextEventTime = std::get<1>(events.top());

#undef QUANT_PDE_TMP_OUTER_TAIL
#define QUANT_PDE_TMP_OUTER_TAIL \
		this->implicitTime = nextEventTime; \
		const Vector *current = &this->iterand(0); \
		Vector transformed; \
		while(!events.empty() && std::get<1>(events.top()) \
				== this->implicitTime) { \
			transformed = (*std::get<2>(events.top()))(*current); \
			events.pop(); \
			current = &transformed; \
		} \
		this->clearNodes(); \
		this->history->clear(); \
		this->history->push(std::make_tuple( \
			this->implicitTime, \
			std::move(transformed) \
		)); \
	} while( Order()(terminalTime(), this->implicitTime) ) \

// TODO: Optimize
#undef QUANT_PDE_TMP_TIMESTEP
#define QUANT_PDE_TMP_TIMESTEP \
	do { \
		dtPrevious = dt; \
		dt = timestep(); \
		Real tmp = this->implicitTime + direction * dt; \
		if( std::abs(tmp - nextEventTime) < epsilon ) { \
			tmp = nextEventTime; \
		} else if( Order()(tmp, nextEventTime) ) { \
			tmp = nextEventTime; \
			dt = direction * (nextEventTime - this->implicitTime); \
		} \
		this->implicitTime = tmp; \
	} while(0)

#undef QUANT_PDE_TMP_NOT_DONE
#define QUANT_PDE_TMP_NOT_DONE \
	( Order()(nextEventTime, this->implicitTime + direction * epsilon) )

/**
 * An iterative method that terminates when a specified terminal time is
 * reached.
 * @tparam Forward True if and only if time moves forward in the relevant
 *                 initial value problem.
 */
template <bool Forward>
class TimeIteration : public Iteration {

	typedef typename std::conditional<
		Forward,
		std::greater<Real>,
		std::less<Real>
	>::type Order;

	typedef std::tuple<
		unsigned,
		Real,
		std::shared_ptr<EventBase>
	> T;

	struct TimeOrder {
		// Returns true if a goes before b in the ordering
		bool operator()(const T &a, const T &b) const {
			return
				Order()( std::get<1>(a), std::get<1>(b) )
				|| (
					std::get<1>(a) == std::get<1>(b)
					&& Order()( std::get<0>(a),
							std::get<0>(b) )
				)
			;
		}
	};

	// If events are set to the same time, ties are broken depending on the
	// order they were added. Events added later are assumed to occur later
	// in time (e.g. handled earlier if Forward == true; later otherwise).

	typedef std::priority_queue<T, std::vector<T>, TimeOrder> PriorityQueue;

	////////////////////////////////////////////////////////////////////////

	static constexpr Real direction = Forward ? 1. : -1.;
	static constexpr Real epsilon = 1e-12;

	////////////////////////////////////////////////////////////////////////

	virtual Vector iterateUntilDone(
		Vector iterand,
		IterationNode &root,
		LinearSolver &solver,
		Real time,
		bool initialized
	) {
		QUANT_PDE_TMP_ITERATE_UNTIL_DONE;
	}

	virtual bool isTimestepTheSame() const {
		return dt == dtPrevious;
	}

	virtual Real timestep() = 0;

	Real initialTime() const;

	Real terminalTime() const;

	////////////////////////////////////////////////////////////////////////

	unsigned id;
	PriorityQueue events;

	Real startTime, endTime, dt, dtPrevious;

public:

	/**
	 * Constructor.
	 */
	TimeIteration(Real startTime, Real endTime) noexcept : id(0),
			startTime(startTime), endTime(endTime) {
		assert(startTime >= 0.);
		assert(startTime < endTime);
	}

	/**
	 * Adds an event to be processed.
	 * @param time The time at which the event occurs.
	 * @param event The event.
	 */
	void add(Real time, std::unique_ptr<EventBase> event) {
		assert(time >= startTime);
		assert(time <= endTime);
		assert(time != initialTime());

		events.emplace( id++, time, std::move(event) );

		// TODO: Test to see if this works
	}

	// TODO: Infer Dimension

	/**
	 * Adds an event to be processed.
	 * @param time The time at which the event occurs.
	 * @param args Arguments passed to Event constructor.
	 * @see QuantPDE::Event
	 */
	template <Index Dimension, typename ...Ts>
	void add(Real time, Ts &&...args) {
		// TODO: Does not work with clang-503.0.40 with -std=c++1y.
		//       Tries to use Event's copy constructor.
		//       Is this a Clang bug?

		assert(time >= startTime);
		assert(time <= endTime);
		assert(time != initialTime());

		events.emplace(
			id++,
			time,
			std::shared_ptr<EventBase>(
				new Event<Dimension>(
					std::forward<Ts>(args)...
				)
			)
		);
	}

};

typedef TimeIteration<false> ReverseTimeIteration;
typedef TimeIteration<true > ForwardTimeIteration;

template <>
Real ReverseTimeIteration::initialTime() const {
	return endTime;
}

template <>
Real ForwardTimeIteration::initialTime() const {
	return startTime;
}

template <>
Real ReverseTimeIteration::terminalTime() const {
	return startTime;
}

template <>
Real ForwardTimeIteration::terminalTime() const {
	return endTime;
}

#undef QUANT_PDE_TMP_SET_TIME
#undef QUANT_PDE_TMP_OUTER_HEAD
#undef QUANT_PDE_TMP_OUTER_TAIL
#undef QUANT_PDE_TMP_TIMESTEP
#undef QUANT_PDE_TMP_HEAD
#undef QUANT_PDE_TMP_TAIL
#undef QUANT_PDE_TMP_NOT_DONE
#undef QUANT_PDE_TMP_ITERATE_UNTIL_DONE

} // QuantPDE

#endif

