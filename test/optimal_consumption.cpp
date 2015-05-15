////////////////////////////////////////////////////////////////////////////////
// optimal_consumption.cpp
// -----------------------
//
// TODO: Description
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>

////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////

using namespace QuantPDE;

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Solution grid
////////////////////////////////////////////////////////////////////////////////

// Grid
RectilinearGrid2 initialGrid(
	// Hand-picked axis, scaled by w0
	// Used in GMWB paper by Zhuliang and Forsyth
	(w0 / 100.) * Axis {
		0., 5., 10., 15., 20., 25.,
		30., 35., 40., 45.,
		50., 55., 60., 65., 70., 72.5, 75., 77.5, 80., 82., 84.,
		86., 88., 90., 91., 92., 93., 94., 95.,
		96., 97., 98., 99., 100.,
		101., 102., 103., 104., 105., 106.,
		107., 108., 109., 110., 112., 114.,
		116., 118., 120., 123., 126.,
		130., 135., 140., 145., 150., 160., 175., 200., 225.,
		250., 300., 500., 750., 1000.
	},

	(w0 / 100.) * Axis {
		0., 5., 10., 15., 20., 25.,
		30., 35., 40., 45.,
		50., 55., 60., 65., 70., 72.5, 75., 77.5, 80., 82., 84.,
		86., 88., 90., 91., 92., 93., 94., 95.,
		96., 97., 98., 99., 100.,
		101., 102., 103., 104., 105., 106.,
		107., 108., 109., 110., 112., 114.,
		116., 118., 120., 123., 126.,
		130., 135., 140., 145., 150., 160., 175., 200., 225.,
		250., 300., 500., 750., 1000.
	}
);

////////////////////////////////////////////////////////////////////////////////
// Operator to discretize
// i.e. Lu - rho + f, where L is the infinitesimal generator of process
//      (B_t, W_t)
////////////////////////////////////////////////////////////////////////////////

class Discretizee final : public RawControlledLinearSystem2_1 {

	const RectilinearGrid2 &grid;

	const Real r, v, mu, rho, gamma;

	const bool controlled;
	const Vector zero;

	template <typename G>
	Discretizee(
		G1 &grid,
		Real interest,
		Real volatility,
		Real drift,
		Real discount,
		Real aversion,
		bool controlled = true
	) noexcept :
		grid( grid ),
		r( interest ),
		v( volatility ),
		mu( drift ),
		rho( discount ),
		gamma( aversion ),
		controlled( controlled ),
		zero( grid.zero() )
	{
	}

	virtual Matrix A(Real) {
		Matrix M = grid.matrix();
		M.reserve( IntegerVector::Constant(grid.size(), 4) );

		// Axes
		const Axis &W = grid[0];
		const Axis &B = grid[1];

		// Control as vector
		Index k = 0;
		const Vector &raw = controlled ? control(0) : zero;

		for(Index j = 0; j < B.size(); ++j) {
			for(Index i = 0; i < W.size(); ++i) {

				const Real
					dWb = W[i  ] - W[i-1],
					dWc = W[i+1] - W[i-1],
					dWf = W[i+1] - W[i  ],
					dBb = B[j  ] - B[j-1],
					dBc = B[j+1] - B[j-1],
					dBf = B[j+1] - B[j  ]
				;

				const Real tmp1 = v * v * W[i] * W[i];
				const Real tmp2 = mu * W[i];

				const Real alpha_common = tmp1 / dWb / dWc;
				const Real  beta_common = tmp1 / dWf / dWc;

				// Central
				Real alpha_i = alpha_common - tmp2 / dWc;
				Real  beta_i =  beta_common + tmp2 / dWc;

				if(alpha_i < 0.) {
					// Forward
					alpha_i = alpha_common;
					 beta_i =  beta_common + tmp2 / dWf;
				} else if(beta_i < 0.) {
					// Backward
					alpha_i = alpha_common - tmp2 / dWb;
					 beta_i =  beta_common;
				}

				if(j == 0) {
					if(i < W.size()) {
						// (B=0, 0<=W<W_max)
					} else {
						// (B=0, W=W_max)
					}
				} else if(j < B.size()) {
					if(i < W.size()) {
						// (0<B<B_max, 0<=W<W_max)
					} else {
						// (0<B<B_max, W=W_max)
					}
				} else {
					if(i < W.size()) {
						// (B=B_max, 0<=W<W_max)
					} else {
						// (B=B_max, W=W_max)
					}
				}

				++k;

			}
		}

		M.makeCompressed();
		return M;
	}

	virtual Vector b(Real) {
		Vector b = grid.vector();

		const Axis &W = grid[0];
		const Axis &B = grid[1];

		// Control as a vector
		Index k = 0;
		const Vector &raw = controlled ? control(0) : zero;

		// B = 0 (no consumption)
		for(Index i = 0; i < W.size(); ++i) {
			b(k) = 0.;
			++k;
		}

		// B > 0
		for(Index j = 1; j < B.size(); ++j) {
			// 0 <= W <= W_max
			for(Index i = 0; i < W.size(); ++i) {
				b(k) = raw(k);
				++k;
			}
		}

		return b;
	}

	virtual bool isATheSame() const {
		return !controlled;
	}

};

int main() {
	return 0;
}
