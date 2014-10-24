////////////////////////////////////////////////////////////////////////////////
// guaranteed_minimum_withdrawal_benefit.cpp
// -----------------------------------------
//
// Computes the price of a Guaranteed Minimum Withdrawal Benefit (GMWB) using an
// implicit, impulse control formulation.
//
// Author: Parsiad Azimzadeh
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>
#include <QuantPDE/Modules/Lambdas>
#include <QuantPDE/Modules/Operators>

////////////////////////////////////////////////////////////////////////////////

#include <algorithm> // max, min, max_element
#include <cmath>     // sqrt
#include <iomanip>  // setw
#include <iostream>  // cout
#include <numeric>   // accumulate
#include <tuple>     // get

////////////////////////////////////////////////////////////////////////////////

using namespace QuantPDE;
using namespace QuantPDE::Modules;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

class Withdrawal final : public ControlledLinearSystem2, public IterationNode {

	RectilinearGrid2 &grid;
	Noncontrollable2 contractAmount, kappa;

	Controllable2 control;

public:

	template <typename G, typename F1, typename F2>
	Withdrawal(G &grid, F1 &&contractAmount, F2 &&kappa) noexcept :
		grid(grid),
		contractAmount(contractAmount),
		kappa(kappa),
		control( Control2(grid) ) {
		registerControl( control );
	}

	inline Real dt() const {
		return time(0) - nextTime();
	}

	inline Real amountWithdrawnPrepenalty(Real t, Real S, Real W) const {
		// Contract withdrawal amount for the period
		const Real G = contractAmount(t, S, W);
		const Real Gdt = G * dt();

		Real gamma;

		#if   defined(GMWB_CONSTANT_WITHDRAWAL)
			gamma = min(W, Gdt);
		#elif defined(GMWB_SURRENDER)
			gamma = W;
		#else
			// Control in [0,2]
			const Real q = control(t, S, W);

			assert(q >= 0.);
			assert(q <= 2.);

			if(W <= Gdt) {
				gamma = (q / 2.) * W;
			} else {
				if(q <= 1.) {
					gamma = q * Gdt;
				} else {
					gamma = Gdt + (q - 1.) * (W - Gdt);
				}
			}
		#endif

		assert(gamma >= 0);
		assert(gamma <= W);

		return gamma;
	}

	virtual Matrix A(Real t) {
		Matrix M(grid.size(), grid.size());
		M.reserve(IntegerVector::Constant(grid.size(), 4));

		Index i = 0;
		for(auto node : grid) {

			const Real S = node[0]; // Investment
			const Real W = node[1]; // Withdrawal

			// Amount withdrawn pre-penalty
			const Real gamma = amountWithdrawnPrepenalty(t, S, W);

			// Interpolation data
			std::array<Real, 2> coordinates {{
				max(S - gamma, 0.),
				W - gamma
			}};
			auto data = interpolationData<2>( grid, coordinates );

			const Index i0 = get<0>( data[0] );
			const Index i1 = get<0>( data[1] );
			const Real  w0 = get<1>( data[0] );
			const Real  w1 = get<1>( data[1] );

			const Index j = grid.index(i0, i1);

			M.insert(i, j                     ) =    w0  *    w1 ;
			M.insert(i, j     + grid[0].size()) =    w0  * (1-w1);
			M.insert(i, j + 1                 ) = (1-w0) *    w1 ;
			M.insert(i, j + 1 + grid[0].size()) = (1-w0) * (1-w1);

			++i;

		}

		M.makeCompressed();
		return grid.identity() - M;
	}

	virtual Vector b(Real t) {
		Vector b = grid.vector();

		for(auto node : accessor(grid, b)) {

			const Real S = (&node)[0]; // Investment
			const Real W = (&node)[1]; // Withdrawal

			// Amount withdrawn, pre-penalty
			const Real gamma = amountWithdrawnPrepenalty(t, S, W);

			// Contract withdrawal amount for the period
			const Real G = contractAmount(t, S, W);
			const Real Gdt = G * dt();

			// Cashflow (including penalty if gamma > Gdt)
			*node = gamma - kappa(t, S, W) * max(gamma - Gdt, 0.);

		}

		return b;
	}

};

////////////////////////////////////////////////////////////////////////////////

int main() {

	// 2014-10-15: Tested without withdrawals; closed-form is
	// dm :=   (log(S0 / (1-kappa) * (r - alpha - 1/2 * sigma * sigma) * T))
	//       / (sigma * sqrt(T))
	// V   =   S0 * exp(-alpha T) * normcdf(dm + sigma * sqrt(T))
	//       + W0   exp(-r   * T) * (1 - kappa) * (1 - normcdf(dm))

	int n = 10; // Initial optimal control partition size
	int N = 32; // Initial number of timesteps

	Real T = 14.28; //10.;
	Real r = .05;
	Real v = .2;

	Real w_0 = 100.;

	Real alpha = 0.036; //0.01389; // Hedging fee

	Real G = 7.; //10.; // Contract rate
	Real kappa = 0.1; // Penalty rate

	int refinement = 5;

	////////////////////////////////////////////////////////////////////////
	// Solution grid
	////////////////////////////////////////////////////////////////////////

	/*
	RectilinearGrid2 grid(
		Axis {
			0., 10., 20.,
			30., 40.,
			50., 60., 70., 75., 80., 84.,
			86., 90., 92., 94.,
			96., 98., 100.,
			102., 104., 106.,
			108., 110., 114.,
			118., 123.,
			130., 140., 150., 175., 225.,
			300., 750., 1000.
		},
		Axis::range(0., 4., 100.)
	);
	*/

	RectilinearGrid2 grid(
		Axis {
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
		Axis::range(0., 2., 100.)
	);

	////////////////////////////////////////////////////////////////////////
	// Table headers
	////////////////////////////////////////////////////////////////////////

	cout.precision(6);
	Real previousValue = nan(""), previousChange = nan("");

	const int td = 20;
	cout
		<< setw(td) << "Nodes"                           << "\t"
		<< setw(td) << "Control Nodes"                   << "\t"
		<< setw(td) << "Time Steps"                      << "\t"
		<< setw(td) << "Value"                           << "\t"
		<< setw(td) << "Mean Inner Iterations"           << "\t"
		<< setw(td) << "Std Inner Iterations"            << "\t"
		<< setw(td) << "Max Inner Iterations"            << "\t"
		<< setw(td) << "Change"                          << "\t"
		<< setw(td) << "Ratio"
		<< endl
	;

	for(
		int l = 0, outer = N, partitionSize = n;
		l < refinement;
		++l, outer *= 2, partitionSize *= 2
	) {

		////////////////////////////////////////////////////////////////
		// Control grid
		////////////////////////////////////////////////////////////////

		// Control partition 0 : 1/n : 1 (MATLAB notation)
		RectilinearGrid1 controls(Axis::range(
			0.,
			2. / (partitionSize - 1),
			2.
		));

		////////////////////////////////////////////////////////////////
		// Iteration tree
		////////////////////////////////////////////////////////////////

		ReverseConstantStepper stepper(
			0.,       // Initial time
			T,        // Expiry time
			T / outer // Timestep size
		);
		ToleranceIteration tolerance;
		stepper.setInnerIteration(tolerance);

		////////////////////////////////////////////////////////////////
		// Linear system tree
		////////////////////////////////////////////////////////////////

		BlackScholes<2, 0> bs(grid, r, v, alpha);
		ReverseRannacher discretization(grid, bs);
		discretization.setIteration(stepper);

		Withdrawal impulse(grid, G, kappa);
		impulse.setIteration(stepper);

		MinPolicyIteration2_1 policy(grid, controls, impulse);
		PenaltyMethod penalty(grid, discretization, policy);

		// TODO: It currently matters what order each linear system is
		//       associated with an iteration; fix this.

		policy.setIteration(tolerance);
		penalty.setIteration(tolerance);

		////////////////////////////////////////////////////////////////
		// Payoff
		////////////////////////////////////////////////////////////////

		Function2 payoff = [=] (Real S, Real W) {
			return max(S, (1 - kappa) * W);
		};

		////////////////////////////////////////////////////////////////
		// Running
		////////////////////////////////////////////////////////////////

		BiCGSTABSolver solver;

		auto V = stepper.solve(
			grid,    // Domain
			payoff,  // Initial condition
			penalty, // Root of linear system tree
			solver   // Linear system solver
		);

		////////////////////////////////////////////////////////////////
		// Print table rows
		////////////////////////////////////////////////////////////////

		auto its = tolerance.iterations();

		Real
			value = V(w_0, w_0),
			var = 0.,
			mean = accumulate(its.begin(),its.end(),0.)/its.size(),
			change = value - previousValue,
			ratio = previousChange / change
		;
		for(auto x : its) { var += (x - mean) * (x - mean); }
		int max = ( *max_element(its.begin(), its.end()) );

		cout
			<< setw(td) << grid.size()   << "\t"
			<< setw(td) << partitionSize << "\t"
			<< setw(td) << outer         << "\t"
			<< setw(td) << value         << "\t"
			<< setw(td) << mean          << "\t"
			<< setw(td) << sqrt(var)     << "\t"
			<< setw(td) << max           << "\t"
			<< setw(td) << change        << "\t"
			<< setw(td) << ratio
			<< endl
		;

		previousChange = change;
		previousValue = value;

		////////////////////////////////////////////////////////////////
		// Refine Solution grid
		////////////////////////////////////////////////////////////////

		grid.refine( RectilinearGrid2::NewTickBetweenEachPair() );
	}

	return 0;
}
