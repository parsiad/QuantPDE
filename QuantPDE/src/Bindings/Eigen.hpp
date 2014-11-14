#include "AssertUnbound.hpp"

#define QUANT_PDE_BOUND

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

#include <unsupported/Eigen/FFT>

#if defined(VIENNACL_WITH_OPENMP) || defined(VIENNACL_WITH_OPENCL) \
		|| defined(VIENNACL_WITH_CUDA)

// Must be set prior to any ViennaCL includes if you want to use ViennaCL
// algorithms on Eigen objects
#define VIENNACL_WITH_EIGEN 1

// ViennaCL headers
//#include <viennacl/linalg/ilu.hpp>
#include <viennacl/linalg/bicgstab.hpp>

#endif

namespace QuantPDE {

typedef Eigen::SparseMatrix<Real>::Index Index;

typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> Vector;
typedef Eigen::SparseMatrix<Real, Eigen::RowMajor> Matrix;

typedef Eigen::VectorXi IntegerVector;

// BiCGSTAB with IncompleteLUT preconditioner
typedef Eigen::BiCGSTAB<Matrix, Eigen::IncompleteLUT<Real>> BiCGSTAB;

////////////////////////////////////////////////////////////////////////////////

class Entry : public Eigen::Triplet<Real> {

public:

	Entry(Index i, Index j, Real v) : Triplet(i, j, v) {
	}

	Real &value() {
		return this->m_value;
	}

};

////////////////////////////////////////////////////////////////////////////////

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

		#ifndef VIENNACL_WITH_EIGEN

			solver.compute(A);
			assert( solver.info() == Eigen::Success );

		#endif

	}

public:

	/**
	 * Constructor.
	 */
	BiCGSTABSolver() noexcept : LinearSolver() {
	}

	virtual Vector solve(const Vector &b, const Vector &guess) const {
		#ifndef VIENNACL_WITH_EIGEN

			Vector v = solver.solveWithGuess(b, guess);
			assert( solver.info() == Eigen::Success );

		#else

			// ViennaCL does not allow for direct input of an
			// initial guess. Solve instead Ax = b - Ax_0 and then
			// translate back x_sol = x + x_0.

			Vector c = b - A * guess;

			//viennacl::linalg::ilut_tag ilut_config;
			//viennacl::linalg::ilut_precond<Matrix> vcl_ilut(vcl_matrix, ilut_config);

			Vector v = viennacl::linalg::solve(
				A, c,
				viennacl::linalg::bicgstab_tag()
				/*, vcl_ilut*/
			) + guess;

		#endif

		return v;
	}

};

}

