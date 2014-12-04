#include "AssertUnbound.hpp"

#define QUANT_PDE_BOUND

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

#include <unsupported/Eigen/FFT>

#if defined(VIENNACL_WITH_OPENMP) || defined(VIENNACL_WITH_OPENCL) \
		|| defined(VIENNACL_WITH_CUDA)

// Must be set prior to any ViennaCL includes if you want to use ViennaCL
// algorithms on Eigen objects
#define VIENNACL_WITH_EIGEN

// ViennaCL headers
#include <viennacl/linalg/ilu.hpp>
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
	BiCGSTABSolver(Real droptol = QuantPDE::droptol) noexcept
			: LinearSolver() {
		solver.preconditioner().setDroptol(droptol);
	}

	virtual Vector solve(const Vector &b, const Vector &guess) const {
		#ifndef VIENNACL_WITH_EIGEN

			Vector v = solver.solveWithGuess(b, guess);
			assert( solver.info() == Eigen::Success );

		#else

			// ViennaCL does not allow for direct input of an
			// initial guess. Solve instead Ax = b - Ax_0 and then
			// translate back x_sol = x + x_0.

			Vector rhs = b - A * guess;

			////////////////////////////////////////////////////////
			// Solve directly on Eigen objects (no preconditioner)
			////////////////////////////////////////////////////////

			//Vector v = viennacl::linalg::solve(
			//	A, rhs,
			//	viennacl::linalg::bicgstab_tag()
			//) + guess;

			////////////////////////////////////////////////////////
			// Copy and solve
			////////////////////////////////////////////////////////

			typedef viennacl::compressed_matrix<Real> vcl_sparse;

			// Eigen -> ViennaCL
			int dim = rhs.size();

			viennacl::vector<Real> vcl_rhs(dim);
			viennacl::copy(rhs, vcl_rhs);

			vcl_sparse vcl_A(dim, dim);
			viennacl::copy(A, vcl_A);

			// Preconditioner
			viennacl::linalg::ilut_tag ilut_config;
			viennacl::linalg::ilut_precond<vcl_sparse> vcl_ilut(
					vcl_A, ilut_config);

			// Solve
			viennacl::vector<Real> vcl_result;
			vcl_result = viennacl::linalg::solve(
				vcl_A,
				vcl_rhs,
				viennacl::linalg::bicgstab_tag(),
				vcl_ilut
			);

			// ViennaCL -> Eigen
			Vector v;
			viennacl::copy(vcl_result, v);
			v += guess; // Adjust for initial guess

		#endif

		return v;
	}

};

}

