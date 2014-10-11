#include "AssertUnbound.hpp"

#define QUANT_PDE_BOUND

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

#include <unsupported/Eigen/FFT>

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

}

