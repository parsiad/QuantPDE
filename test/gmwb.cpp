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
#include <iomanip>   // setw
#include <iostream>  // cout
#include <memory>    // unique_ptr
#include <numeric>   // accumulate
#include <tuple>     // get

////////////////////////////////////////////////////////////////////////////////

using namespace QuantPDE;
using namespace QuantPDE::Modules;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

class ContinuousWithdrawal final : public ControlledLinearSystem2,
		public IterationNode {

	RectilinearGrid2 &grid;
	Noncontrollable2 contractAmount;

	Controllable2 control;

public:

	template <typename G, typename F1>
	ContinuousWithdrawal(G &grid, F1 &&contractAmount) noexcept
			: grid(grid), contractAmount(contractAmount),
			control( Control2(grid) ) {
		registerControl(control);
	}

	inline Real dt() const {
		return time(0) - nextTime();
	}

	virtual Matrix A(Real t) {
		Matrix M = grid.matrix();
		M.reserve(IntegerVector::Constant(grid.size(), 4));

		Index k = 0;
		for(auto node : grid) {

			const Real S = node[0]; // Investment
			const Real W = node[1]; // Withdrawal

			const Real G = contractAmount(t, S, W);
			const Real Gdt = G * dt();
			const Real gamma = control(t, S, W) * Gdt;

			// TODO: Remove branching
			if(W > epsilon) {
				if(S > epsilon) {
					M.insert(k, k) =  2. * gamma;
					M.insert(k, k - 1) = -1. * gamma;
					M.insert(k, k - grid[0].size()) = -1.
							* gamma;
				} else {
					// S ~= 0
					M.insert(k, k) =  1. * gamma;
					M.insert(k, k - grid[0].size()) = -1.
							* gamma;
				}
			}

			++k;

		}

		M.makeCompressed();
		return M;
	}

	virtual Vector b(Real t) {
		Vector b = grid.vector();

		for(auto node : accessor(grid, b)) {
			const Real S = (&node)[0]; // Investment
			const Real W = (&node)[1]; // Withdrawal

			const Real G = contractAmount(t, S, W);
			const Real Gdt = G * dt();

			*node = control(t, S, W) * Gdt;
		}

		return b;
	}

};

class ImpulseWithdrawal final : public ControlledLinearSystem2,
		public IterationNode {

	RectilinearGrid2 &grid;
	Noncontrollable2 contractAmount, kappa;

	Controllable2 control;

public:

	template <typename G, typename F1, typename F2>
	ImpulseWithdrawal(G &grid, F1 &&contractAmount, F2 &&kappa) noexcept
			: grid(grid), contractAmount(contractAmount),
			kappa(kappa), control( Control2(grid) ) {
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

		// Control in [0,2]
		const Real q = control(t, S, W);

		assert(q >= 0.);
		assert(q <= 2.);

		if( q <= 1. ) {
			// Nonpenalty
			gamma = q * min(W, Gdt);
		} else {
			if( W > Gdt ) {
				// Penalty
				gamma = Gdt + (q - 1.) * (W - Gdt);
			} else {
				gamma = min(W, Gdt);
			}
		}

		assert(gamma >= 0);
		assert(gamma <= W);

		return gamma;
	}

	virtual Matrix A(Real t) {
		Matrix M = grid.matrix();
		M.reserve(IntegerVector::Constant(grid.size(), 4));

		Index k = 0;
		for(auto node : grid) {

			const Real S = node[0]; // Investment
			const Real W = node[1]; // Withdrawal

			// Amount withdrawn pre-penalty
			const Real gamma = amountWithdrawnPrepenalty(t, S, W);

			const Real Splus = max(S - gamma, 0.);
			const Real Wplus = W - gamma;

			// Interpolation data
			auto data = interpolationData<2>(grid, {{Splus,Wplus}});

			const Index i0 = get<0>( data[0] );
			const Index i1 = get<0>( data[1] );
			const Real  w0 = get<1>( data[0] );
			const Real  w1 = get<1>( data[1] );

			assert( (grid[0][i0+1] - Splus)
					/ (grid[0][i0+1] - grid[0][i0]) == w0 );
			assert( (grid[1][i1+1] - Wplus)
					/ (grid[1][i1+1] - grid[1][i1]) == w1 );

			const Index j = grid.index(i0, i1);

			M.insert(k, j                     ) =    w0  *    w1 ;
			M.insert(k, j     + grid[0].size()) =    w0  * (1-w1);
			M.insert(k, j + 1                 ) = (1-w0) *    w1 ;
			M.insert(k, j + 1 + grid[0].size()) = (1-w0) * (1-w1);

			++k;

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
			*node = gamma - kappa(t, S, W) * max(gamma - Gdt, 0.)
					- epsilon;

			assert( *node >= -epsilon );

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
	Real kappa = 0.; // Penalty rate

	int refinement = 5;

	////////////////////////////////////////////////////////////////////////
	// Solution grid
	////////////////////////////////////////////////////////////////////////

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
		#if   defined(GMWB_SURRENDER)
			RectilinearGrid1 impulseControls(Axis { 2. });
		#elif defined(GMWB_CONTRACT_WITHDRAWAL)
			RectilinearGrid1 impulseControls(Axis { 1. });
		#else
			RectilinearGrid1 impulseControls(Axis::range(
				0.,
				2. / (partitionSize - 1),
				2.
			));
		#endif

		RectilinearGrid1 continuousControls( Axis { 0., 1. } );

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

		// Impulse withdrawal policy iteration
		ImpulseWithdrawal impulseWithdrawal(grid, G, kappa);
		impulseWithdrawal.setIteration(stepper);
		MinPolicyIteration2_1 impulsePolicy(
			grid,
			impulseControls,
			impulseWithdrawal
		);

		// Continuous withdrawal policy iteration
		ContinuousWithdrawal *cw = new ContinuousWithdrawal(grid, G);
		cw->setIteration(stepper);
		MinPolicyIteration2_1 continuousPolicy(
			grid,
			continuousControls,
			*cw
		);

		// Linear system sum
		unique_ptr<LinearSystem> bs(
			new BlackScholes<2, 0>(
				grid,
				r, v, alpha
			)
		);
		unique_ptr<LinearSystem> continuousWithdrawal(cw); cw = nullptr;
		auto sum = std::move(bs) + std::move(continuousWithdrawal);

		// Discretization
		ReverseRannacher discretization(grid, sum);
		discretization.setIteration(stepper);

		// Penalty method
		PenaltyMethod penalty(grid, discretization, impulsePolicy);

		// TODO: It currently matters what order each linear system is
		//       associated with an iteration; fix this.

		impulsePolicy.setIteration(tolerance);
		continuousPolicy.setIteration(tolerance);
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

		RectilinearGrid2 printGrid(
			Axis::range(0., 25., 200.),
			Axis { 100. }
		);
		cout << accessor( printGrid, V ) << endl;

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
			<< setw(td) << grid.size()            << "\t"
			<< setw(td) << impulseControls.size() << "\t"
			<< setw(td) << outer                  << "\t"
			<< setw(td) << value                  << "\t"
			<< setw(td) << mean                   << "\t"
			<< setw(td) << sqrt(var)              << "\t"
			<< setw(td) << max                    << "\t"
			<< setw(td) << change                 << "\t"
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
