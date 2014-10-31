////////////////////////////////////////////////////////////////////////////////
// gmwb.cpp
// --------
//
// Computes the price of a GMWB using several different formulations.
//
// Author: Parsiad Azimzadeh
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>
#include <QuantPDE/Modules/Lambdas>
#include <QuantPDE/Modules/Operators>

////////////////////////////////////////////////////////////////////////////////

#include <algorithm> // max, min
#include <climits>   // INT_MAX
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

class ImpulseWithdrawal final : public ControlledLinearSystem2,
		public IterationNode {

	RectilinearGrid2 &grid;
	Noncontrollable2 kappa;

	Controllable2 control;

public:

	template <typename G, typename F1>
	ImpulseWithdrawal(G &grid, F1 &&kappa) noexcept
			: grid(grid), kappa(kappa), control( Control2(grid) ) {
		registerControl( control );
	}

	virtual Matrix A(Real t) {
		Matrix M = grid.matrix();
		M.reserve(IntegerVector::Constant(grid.size(), 4));

		Index k = 0;
		for(auto node : grid) {
			const Real S = node[0]; // Investment
			const Real W = node[1]; // Withdrawal

			// Amount withdrawn pre-penalty
			const Real gamma = control(t, S, W) * W;

			const Real Splus = std::max(S - gamma, 0.);
			const Real Wplus = W - gamma;

			// Interpolation data
			auto data = interpolationData<2>(grid, {{Splus,Wplus}});

			const Index i0 = std::get<0>( data[0] );
			const Index i1 = std::get<0>( data[1] );
			const Real  w0 = std::get<1>( data[0] );
			const Real  w1 = std::get<1>( data[1] );

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
			const Real gamma = control(t, S, W) * W;

			// Cashflow minus adjustment
			*node = (1 - kappa(t, S, W)) * gamma - epsilon;
		}

		return b;
	}

};

class ContinuousWithdrawal final : public ControlledLinearSystem2,
		public IterationNode {

	RectilinearGrid2 &grid;
	Noncontrollable2 contractRate;

	Controllable2 control;

public:

	template <typename G, typename F1>
	ContinuousWithdrawal(G &grid, F1 &&contractRate) noexcept
			: grid(grid), contractRate(contractRate),
			control( Control2(grid) ) {
		registerControl(control);
	}

	inline Real dt() const {
		return time(0) - nextTime();
	}

	virtual Matrix A(Real t) {
		Matrix M = grid.matrix();

		auto M_G = grid.indexer(M);

		const Axis &S = grid[0];
		const Axis &W = grid[1];

		M.reserve(IntegerVector::Constant(grid.size(), 3));

		// W > 0
		for(Index j = 1; j < W.size(); ++j) {
			// S = 0
			{
				const Real G  = contractRate(t, S[0], W[j]);
				const Real g  =      control(t, S[0], W[j]) * G;

				const Real tW = g / (W[j] - W[j-1]);

				M_G(0, j, 0, j    ) =   tW;
				M_G(0, j, 0, j - 1) = - tW;
			}

			// S > 0
			for(Index i = 1; i < S.size(); ++i) {
				const Real G  = contractRate(t, S[i], W[j]);
				const Real g  =      control(t, S[i], W[j]) * G;

				const Real tW = g / (W[j] - W[j-1]);
				const Real tS = g / (S[i] - S[i-1]);

				M_G(i, j, i    , j    ) = + tW + tS;
				M_G(i, j, i    , j - 1) = - tW;
				M_G(i, j, i - 1, j    ) =      - tS;
			}
		}

		M.makeCompressed();
		return M;
	}

	virtual Vector b(Real t) {
		Vector b = grid.vector();

		auto b_G = grid.indexer(b);

		const Axis &S = grid[0];
		const Axis &W = grid[1];

		// W = 0 (no withdrawal)
		for(Index i = 0; i < S.size(); ++i) {
			b_G(i, 0) = 0.;
		}

		// W > 0
		for(Index j = 1; j < W.size(); ++j) {
			// S >= 0
			for(Index i = 0; i < S.size(); ++i) {
				const Real G = contractRate(t, S[i], W[j]);
				const Real g =      control(t, S[i], W[j]) * G;

				b_G(i, j) = g;
			}
		}

		return b;
	}

};

////////////////////////////////////////////////////////////////////////////////

int main() {

	////////////////////////////////////////////////////////////////////////
	// Methods
	////////////////////////////////////////////////////////////////////////

	constexpr int SEMI_LAGRANGIAN_WITHDRAWAL_NO_PENALTY = 1 << 0;
	constexpr int SEMI_LAGRANGIAN_WITHDRAWAL_AT_PENALTY = 1 << 1;

	constexpr int EXPLICIT =
			  SEMI_LAGRANGIAN_WITHDRAWAL_AT_PENALTY
			| SEMI_LAGRANGIAN_WITHDRAWAL_NO_PENALTY;

	constexpr int IMPLICIT = 0;

	////////////////////////////////////////////////////////////////////////

	int method;

	//method = SEMI_LAGRANGIAN_WITHDRAWAL_NO_PENALTY;
	method = SEMI_LAGRANGIAN_WITHDRAWAL_AT_PENALTY;
	//method = EXPLICIT;
	//method = IMPLICIT;

	////////////////////////////////////////////////////////////////////////

	Real T = 10.; // 14.28;
	Real r = .05;
	Real v = .2;

	Real alpha = 0.01389; // 0.036; // Hedging fee

	Real G = 10.; // 7.; // Contract rate
	Real kappa = 0.1; // Penalty rate

	Real w_0 = 100.; // Initial value of the account

	int N = 32; // Initial number of timesteps
	int M = 2; // Initial control set partition
	int Mmax = INT_MAX; //16; // Maximum control set partition size

	int R = 1; // Minimum level of refinement
	int Rmax = 10; // Maximum level of refinement

	assert(r >= 0.);

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
		Axis::range(0., 5., 100.)
	);

	////////////////////////////////////////////////////////////////////////
	// Table headers
	////////////////////////////////////////////////////////////////////////

	cout.precision(6);
	Real previousValue = nan(""), previousChange = nan("");

	const int td = 20;
	cout
		<< setw(td) << "Nodes"           << "\t"
		<< setw(td) << "Control Nodes"   << "\t"
		<< setw(td) << "Time Steps"      << "\t"
		<< setw(td) << "Value"           << "\t"
		<< setw(td) << "Mean Iterations" << "\t"
		<< setw(td) << "Std Iterations"  << "\t"
		<< setw(td) << "Max Iterations"  << "\t"
		<< setw(td) << "Change"          << "\t"
		<< setw(td) << "Ratio"
		<< endl
	;

	////////////////////////////////////////////////////////////////////////
	// Refinement loop
	////////////////////////////////////////////////////////////////////////

	for(
		int l = R;
		l <= Rmax;
		++l, N *= 2, M *= 2
	) {

		M = min(M, Mmax);

		////////////////////////////////////////////////////////////////
		// Iteration tree
		////////////////////////////////////////////////////////////////

		ReverseConstantStepper stepper(
			0.,    // Initial time
			T,     // Expiry time
			T / N  // Timestep size
		);

		// Tolerance iteration
		unique_ptr<ToleranceIteration> tolerance;
		if(method != EXPLICIT) {
			tolerance = unique_ptr<ToleranceIteration>(
					new ToleranceIteration());
			stepper.setInnerIteration(*tolerance);
		}

		////////////////////////////////////////////////////////////////
		// Linear system tree
		////////////////////////////////////////////////////////////////

		typedef ReverseLinearBDFOne Discretization;

		unique_ptr<BlackScholes<2, 0>> blackScholes(
			new BlackScholes<2,0>(
				grid,
				r, v, alpha
			)
		);

		unique_ptr<RectilinearGrid1> continuousControls;
		unique_ptr<ContinuousWithdrawal> continuousWithdrawal;
		unique_ptr<MinPolicyIteration2_1> continuousPolicy;
		unique_ptr<LinearSystemSum> sum;
		LinearSystem *discretizee;
		if(method & SEMI_LAGRANGIAN_WITHDRAWAL_NO_PENALTY) {
			// Using semi-lagrangian for continuous withdrawal
			discretizee = blackScholes.get();
		} else {
			// Continuous control grid
			continuousControls = unique_ptr<RectilinearGrid1>(
					new RectilinearGrid1(Axis { 0., 1. }));

			// Continuous withdrawal
			continuousWithdrawal = unique_ptr<ContinuousWithdrawal>(
					new ContinuousWithdrawal(grid, G));
			continuousWithdrawal->setIteration(stepper);

			// Continuous withdrawal policy iteration
			continuousPolicy = unique_ptr<MinPolicyIteration2_1>(
				new MinPolicyIteration2_1(
					grid,
					*continuousControls,
					*continuousWithdrawal
				)
			);
			continuousPolicy->setIteration(*tolerance);

			// Sum of linear systems
			sum = unique_ptr<LinearSystemSum>(
				new LinearSystemSum(
					  std::move(blackScholes    )
					+ std::move(continuousPolicy)
				)
			);

			discretizee = sum.get();
		}

		Discretization discretization(grid, *discretizee);
		discretization.setIteration(stepper);

		unique_ptr<RectilinearGrid1> impulseControls;
		unique_ptr<ImpulseWithdrawal> impulseWithdrawal;
		unique_ptr<MinPolicyIteration2_1> impulsePolicy;
		unique_ptr<PenaltyMethod> penalty;
		IterationNode *root;
		if(method & SEMI_LAGRANGIAN_WITHDRAWAL_AT_PENALTY) {
			// No impulse root
			root = &discretization;
		} else {
			// Using semi-lagrangian for withdrawal at a penalty

			// Impulse control grid
			impulseControls = unique_ptr<RectilinearGrid1>(
				new RectilinearGrid1(
					Axis::range(
						1. / M,
						1. / M,
						1.
					)
				)
			);

			// Impulse withdrawal
			impulseWithdrawal = unique_ptr<ImpulseWithdrawal>(
					new ImpulseWithdrawal(grid, kappa));

			// Impulse withdrawal policy iteration
			impulsePolicy = unique_ptr<MinPolicyIteration2_1>(
				new MinPolicyIteration2_1(
					grid,
					*impulseControls,
					*impulseWithdrawal
				)
			);
			impulsePolicy->setIteration(*tolerance);

			// Penalty method
			penalty = unique_ptr<PenaltyMethod>(
				new PenaltyMethod(
					grid,
					discretization,
					*impulsePolicy
				)
			);
			penalty->setIteration(*tolerance);

			// Impulse root
			root = penalty.get();
		}

		////////////////////////////////////////////////////////////////
		// Exercise events
		////////////////////////////////////////////////////////////////

		auto withdrawal = [=] (const Interpolant2 &V, Real S, Real W) {
			Real best = V(S, W);

			// Contract withdrawal amount
			const Real Gdt = G * T / N;

			/*
			#if   defined(GMWB_CONTRACT_WITHDRAWAL)
				// Constant withdrawal
				const Real gamma = min(W, Gdt);
				const Real interp = V(
					max(S - gamma, 0.),
					W - gamma
				);
				const Real cashflow = gamma;

				best = interp + cashflow;
			#elif defined(GMWB_SURRENDER)
				const Real gamma = W;

				const Real interp = V(
					max(S - gamma, 0.),
					W - gamma
				);
				const Real cashflow = gamma  - kappa
						* max(gamma - Gdt, 0.);

				best = interp + cashflow;
			#else
			*/

			if(method & SEMI_LAGRANGIAN_WITHDRAWAL_NO_PENALTY) {
				// Nonpenalty
				const Real beta = min(W, Gdt);
				//for(int i = 1; i <= 1; ++i) {
					//const Real gamma = beta * i / M;
					const Real gamma = beta;
					const Real interp = V(
						max(S - gamma, 0.),
						W - gamma
					);
					const Real cashflow = gamma;
					const Real newValue = interp + cashflow;
					if(newValue > best) {
						best = newValue;
					}
				//}
			}

			// Penalty
			if(
				method & SEMI_LAGRANGIAN_WITHDRAWAL_AT_PENALTY
				&& W > Gdt
			) {
				for(int i = 1; i <= M; ++i) {
					const Real gamma = Gdt + (W-Gdt) * i/M;
					const Real interp = V(
						max(S - gamma, 0.),
						W - gamma
					);
					const Real cashflow = gamma - kappa
							* (gamma - Gdt);
					const Real newValue = interp + cashflow;
					if(newValue > best) {
						best = newValue;
					}
				}
			}

			//#endif

			return best;
		};

		if(method != IMPLICIT) {
			for(int m = 0; m < N; ++m) {
				stepper.add(
					// Time at which the event takes place
					T / N * m,

					withdrawal,

					// Spatial grid to interpolate on
					grid
				);
			}
		}

		////////////////////////////////////////////////////////////////
		// Payoff
		////////////////////////////////////////////////////////////////

		Function2 payoff = [=] (Real S, Real W) {
			return max(S, (1 - kappa) * W - epsilon);
		};

		////////////////////////////////////////////////////////////////
		// Running
		////////////////////////////////////////////////////////////////

		BiCGSTABSolver solver;

		auto V = stepper.solve(
			grid,   // Domain
			payoff, // Initial condition
			*root,  // Root of linear system tree
			solver  // Linear system solver
		);

		////////////////////////////////////////////////////////////////
		// Print table rows
		////////////////////////////////////////////////////////////////

		/*
		RectilinearGrid2 printGrid(
			Axis::range(0., 25., 200.),
			Axis { 100. }
		);
		cout << accessor( printGrid, V ) << endl;
		*/

		Real
			value = V(w_0, w_0),
			change = value - previousValue,
			ratio = previousChange / change,
			var = 0., mean = 1.
		;

		int max = 1;

		if(method != EXPLICIT) {
			auto its = tolerance->iterations();

			var = 0.;
			mean = accumulate(its.begin(),its.end(),0.)/its.size();

			for(auto x : its) { var += (x - mean) * (x - mean); }
			max = ( *max_element(its.begin(), its.end()) );
		}

		cout
			<< setw(td) << grid.size() << "\t"
			<< setw(td) << M           << "\t"
			<< setw(td) << N           << "\t"
			<< setw(td) << value       << "\t"
			<< setw(td) << mean        << "\t"
			<< setw(td) << sqrt(var)   << "\t"
			<< setw(td) << max         << "\t"
			<< setw(td) << change      << "\t"
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
