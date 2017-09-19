////////////////////////////////////////////////////////////////////////////////
// glwb.cpp
// --------
//
// Computes the price of a guaranteed lifelong withdrawal benefits contract as
// specified in [1]
//
// [1] Azimzadeh, Parsiad, and Peter A. Forsyth. "The existence of optimal
// bang-bang controls for GMxB contracts." SIAM Journal on Financial Mathematics
// 6.1 (2015): 117-139.
//
// Author: Parsiad Azimzadeh
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>
#include <QuantPDE/Modules/Operators>
#include <QuantPDE/Modules/Utilities>

using namespace QuantPDE;
using namespace QuantPDE::Modules;

///////////////////////////////////////////////////////////////////////////////

#include <iostream> // cerr
#include <cmath> // abs

using namespace std;

///////////////////////////////////////////////////////////////////////////////

Real *R, *P;

class GLWBOperator : public BlackScholes1 {

public:

	GLWBOperator(
		RectilinearGrid1 &grid,
		Real r, Real vol, Real management_rate
	) noexcept : BlackScholes1(grid, r, vol, management_rate) {}

	virtual Vector b(Real t) {
		Vector S_vec = G.image([] (Real S) { return S; });
		int n = floor(t);
		return -(R[n+1] - R[n]) * S_vec;
	}
};

///////////////////////////////////////////////////////////////////////////////

Real r, vol, S_0, management_rate, withdrawal_rate, bonus_rate;
int N, T, kn, k0;
RectilinearGrid1 *grid, *mortality, *penalty_rates;

// Used to store the optimal control
Vector optimal_control;

ResultsTuple1 run(int k) {

	// Refinement
	int factor = 1;
	for(int i = 0; i < k; ++i) {
		factor *= 2;
	}
	const int N2k = N * factor;
	auto refined_grid = grid->refined(k);

	SparseLUSolver solver;
	typedef ReverseBDFTwo Discretization;
    //typedef ReverseBDFOne Discretization;

	////////////////////////////////////////////////////////////////////////

	GLWBOperator op(refined_grid, r, vol, management_rate);

	ReverseConstantStepper stepper(0., (Real) T, (Real) T / N2k);

	Discretization discretization(refined_grid, op);
	discretization.setIteration(stepper);

	////////////////////////////////////////////////////////////////////////

	const Real Q = S_0;
	for(int n = 1; n <= T; ++n) {
		auto event = [Q, n, k] (const Interpolant1 &V, Real S) {
			const Real GQ = withdrawal_rate * Q;

			// Nonwithdrawal
			Real best = (1. + bonus_rate) * V(S / (1. + bonus_rate));
			int optimal_control = 0;

			// Withdraw
			{
				const Real S_plus = max(S - GQ, 0.);
				const Real tmp = V(S_plus) + R[n] * GQ;
				if(tmp > best) {
					best = tmp;
					optimal_control = 1;
				}
			}

			// Surrender
			if(S > GQ) {
				const Real tmp = R[n] * (
					GQ
					+ ( 1. - (*penalty_rates)[0][n-1] ) * (S - GQ)
				);
				if(tmp > best) {
					best = tmp;
					optimal_control = 2;
				}
			}

			// Print optimal control
			/*if(k == kn && n == 1) {
				cerr << S << "\t" << optimal_control << endl;
			}*/

			return best;
		};
		stepper.add((Real) n, event, refined_grid);
	}

	////////////////////////////////////////////////////////////////////////

	auto solution = stepper.solve(
		refined_grid,
		[] (Real) { return 0.; },
		discretization,
		solver
	);

	return ResultsTuple1(
		{(Real) refined_grid.size(), (Real) N2k },
		solution, S_0
	);

}

int main(int argc, char **argv) {
	// Parse configuration file
	Configuration configuration = getConfiguration(argc, argv);

	// Get options
	Real S_max, S_min, dS;
	kn = getInt(configuration, "maximum_refinement", 8);
	k0 = getInt(configuration, "minimum_refinement", 0);
	r = getReal(configuration, "interest_rate", .04);
	vol = getReal(configuration, "volatility", .2);
	S_0 = getReal(configuration, "asset_price", 100.);
	management_rate = getReal(configuration, "management_rate", .015);
	withdrawal_rate = getReal(configuration, "withdrawal_rate", .05);
	bonus_rate = getReal(configuration, "bonus_rate", .06);
	S_min = getReal(configuration, "print_asset_price_minimum", 0.);
	S_max = getReal(configuration, "print_asset_price_maximum", S_0 * 2.);
	dS = getReal(configuration, "print_asset_price_step_size", S_0 / 10.);
	N = getInt(configuration, "initial_number_of_timesteps", 12);

	// PDE grid
	RectilinearGrid1 default_grid(S_0 * Axis::special);
	RectilinearGrid1 tmp1 = getGrid(configuration, "initial_grid", default_grid);
	grid = &tmp1;

	// Mortality data (default is DAV2004R)
	RectilinearGrid1 default_mortality( Axis {
		0.008886, 0.009938, 0.011253, 0.012687, 0.014231,
		0.015887, 0.017663, 0.019598, 0.021698, 0.023990,
		0.026610, 0.029533, 0.032873, 0.036696, 0.041106,
		0.046239, 0.052094, 0.058742, 0.066209, 0.074583,
		0.083899, 0.094103, 0.105171, 0.116929, 0.129206,
		0.141850, 0.154860, 0.168157, 0.181737, 0.195567,
		0.209614, 0.223854, 0.238280, 0.252858, 0.267526,
		0.278816, 0.293701, 0.308850, 0.324261, 0.339936,
		0.355873, 0.372069, 0.388523, 0.405229, 0.422180,
		0.439368, 0.456782, 0.474411, 0.492237, 0.510241,
		0.528401, 0.546689, 0.565074, 0.583517, 0.601976,
		0.620400, 1.000000
	} );
	RectilinearGrid1 tmp2 = getGrid(configuration, "mortality_data", default_mortality);
	mortality = &tmp2;
	T = mortality->size();

	// Determine survival rate
	R = new Real[T+1]; // R[0], ..., R[T]
	R[0] = 1.;
	for(int n = 1; n <= T; ++n) {
		R[n] = R[n-1] * (1. - (*mortality)[0][n-1]);
	}

	// Default penalty rates
	RectilinearGrid1 default_penalty( Axis {
		0.03, 0.02, 0.01, 0., 0.,
		0., 0., 0., 0., 0.,
		0., 0., 0., 0., 0.,
		0., 0., 0., 0., 0.,
		0., 0., 0., 0., 0.,
		0., 0., 0., 0., 0.,
		0., 0., 0., 0., 0.,
		0., 0., 0., 0., 0.,
		0., 0., 0., 0., 0.,
		0., 0., 0., 0., 0.,
		0., 0., 0., 0., 0.,
		0., 0.
	} );
	RectilinearGrid1 tmp3 = getGrid(configuration, "penalty_data", default_penalty);
	penalty_rates = &tmp3;

	// Print configuration file
	cerr << configuration << endl << endl;

	// Run and print results
	ResultsBuffer1 buffer(
		run,
		{ "Nodes", "Steps" },
		kn, k0
	);
	buffer.setPrintGrid( RectilinearGrid1(Axis::range(S_min, dS, S_max)) );
	buffer.stream();

	// Deallocate
	delete [] R;

	return 0;
}

