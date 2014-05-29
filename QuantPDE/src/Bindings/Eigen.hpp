#include "AssertUnbound.hpp"

#define QUANT_PDE_BOUND

#include <eigen3/Eigen/SparseCore>

namespace QuantPDE {

typedef double Real;

typedef Eigen::SparseMatrix<Real>::Index Index;

typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> Vector;
typedef Eigen::SparseMatrix<Real, Eigen::RowMajor> Matrix;

typedef Eigen::VectorXi IntegerVector;

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

// TODO: Find a place for this

/*
class BiCGSTABSolver : public LinearSolver {

	// BiCGSTAB with IncompleteLUT preconditioner
	BiCGSTAB<Matrix, IncompleteLUT<Real>> solver;

	virtual void initialize(const Matrix &A) {
		solver.compute(A);
	}

	virtual Vector solve(const Vector &b, const Vector &x0) const {
		return solver.solveWithGuess(b, x0);
	}

};
*/

}

