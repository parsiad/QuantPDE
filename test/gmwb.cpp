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
#include <cstdlib>   // abs
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

class ImpulseWithdrawal final : public ControlledLinearSystem2 {

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

			const Real Splus = max(S - gamma, 0.);
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
			*node = (1 - kappa(t, S, W)) * gamma;
		}

		return b;
	}

};

class ContinuousWithdrawal final : public ControlledLinearSystem2 {

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

	virtual Matrix A(Real t) {
		Matrix M = grid.matrix();
		M.reserve(IntegerVector::Constant(grid.size(), 4));

		const Axis &S = grid[0];
		const Axis &W = grid[1];

		// Control as a vector
		const Vector &raw = ((const Control2 *) control.get())->raw();
		Index k = S.size();

		// W > 0
		for(Index j = 1; j < W.size(); ++j) {
			// S = 0
			{
				const Real G  = contractRate(t, S[0], W[j]);
				const Real g = raw(k) * G;

				const Real tW = g / (W[j] - W[j-1]);

				M.insert(k, k           ) = +tW;
				M.insert(k, k - S.size()) = -tW;

				++k;
			}

			#if 0
			// S > 0
			for(Index i = 1; i < S.size(); ++i) {
				const Real G  = contractRate(t, S[i], W[j]);
				const Real g = raw(k) * G;

				const Real tW = g / (W[j] - W[j-1]);
				const Real tS = g / (S[i] - S[i-1]);

				M.insert(k, k           ) = + tW + tS;
				M.insert(k, k - S.size()) = - tW     ;
				M.insert(k, k - 1       ) =      - tS;

				++k;
			}
			#endif

			//#if 0
			// 0 < S < S_max
			for(Index i = 1; i < S.size() - 1; ++i) {
				const Real G  = contractRate(t, S[i], W[j]);
				const Real g = raw(k) * G;

				const Real tW = g / (W[j  ] - W[j-1]);
				const Real tS = g / (S[i+1] - S[i-1]);

				M.insert(k, k + 1       ) =      + tS;
				M.insert(k, k           ) = + tW     ;
				M.insert(k, k - S.size()) = - tW     ;
				M.insert(k, k - 1       ) =      - tS;

				++k;
			}

			// S = S_max
			{
				const Index i = S.size() - 1;

				const Real G  = contractRate(t, S[i], W[j]);
				const Real g = raw(k) * G;

				const Real tW = g / (W[j] - W[j-1]);
				const Real tS = g / (S[i] - S[i-1]);

				M.insert(k, k           ) = + tW + tS;
				M.insert(k, k - S.size()) = - tW     ;
				M.insert(k, k - 1       ) =      - tS;

				++k;
			}
			//#endif
		}

		M.makeCompressed();
		return M;
	}

	virtual Vector b(Real t) {
		Vector b = grid.vector();

		const Axis &S = grid[0];
		const Axis &W = grid[1];

		// Control as a vector
		const Vector &raw = ((const Control2 *) control.get())->raw();
		Index k = 0;

		// W = 0 (no withdrawal)
		for(Index i = 0; i < S.size(); ++i) {
			b(k) = 0.;
			++k;
		}

		// W > 0
		for(Index j = 1; j < W.size(); ++j) {
			// S >= 0
			for(Index i = 0; i < S.size(); ++i) {
				const Real G = contractRate(t, S[i], W[j]);
				const Real g = raw(k) * G;
				b(k) = g;
				++k;
			}
		}

		return b;
	}

};

inline Real Vminus(const Interpolant2 &V, Real S, Real W, Real gamma) {
	return V(
		max(S - gamma, 0.),
		W - gamma
	);
}

////////////////////////////////////////////////////////////////////////////////
// Methods
////////////////////////////////////////////////////////////////////////////////

constexpr int SEMI_LAGRANGIAN_WITHDRAWAL_NO_PENALTY = 1 << 0;
constexpr int SEMI_LAGRANGIAN_WITHDRAWAL_AT_PENALTY = 1 << 1;

constexpr int EXPLICIT =
		  SEMI_LAGRANGIAN_WITHDRAWAL_AT_PENALTY
		| SEMI_LAGRANGIAN_WITHDRAWAL_NO_PENALTY;

constexpr int IMPLICIT = 0;

////////////////////////////////////////////////////////////////////////////////
// Options
////////////////////////////////////////////////////////////////////////////////

//int method = SEMI_LAGRANGIAN_WITHDRAWAL_NO_PENALTY;
//int method = SEMI_LAGRANGIAN_WITHDRAWAL_AT_PENALTY;
//int method = EXPLICIT;
int method = IMPLICIT;

Real T = 10.; // 14.28;
Real r = .05;
Real v = .2;

Real alpha = 0.01389; // 0.036; // Hedging fee

Real G = 10.; // 7.; // Contract rate
Real kappa = 0.1; // Penalty rate

Real w0 = 100.; // Initial value of the account

int N = 32; // Initial number of timesteps
int M = 2; // Initial control set partition
int Mmax = INT_MAX; //16; // Maximum control set partition size

int Rmin = 0;
int Rmax = 10; // Maximum level of refinement

bool newton = true;

////////////////////////////////////////////////////////////////////////////////
// Solution grid
////////////////////////////////////////////////////////////////////////////////

// Peter's grid
RectilinearGrid2 grid(
	Axis {
		0., 5., 10., 15., 20., 25.,
		30., 35., 40., 45.,
		50., 55., 60., 65., 70., 72.5, 75., 77.5, 80., 82., 84.,
		86., 88., 90.,91., 92., 93., 94., 95.,
		96., 97., 98., 99., 100.,
		101., 102., 103., 104., 105., 106.,
		107., 108., 109., 110., 112., 114.,
		116., 118., 120., 123., 126.,
		130., 135., 140., 145., 150., 160., 175., 200., 225.,
		250., 300., 500.,750., 1000.
	},
	Axis::range(0., 2., 100.)
);

/*
constexpr int points1 = 64;
constexpr int points2 = 50;

RectilinearGrid2 grid(
	Axis::cluster(0., 1000., points1, w0, w0 / 5.),
	Axis::cluster(0.,  100., points2, w0, w0 / 5.)
);
*/

////////////////////////////////////////////////////////////////////////////////

std::tuple<Real, Real, Real, int> solve(Real alpha) {

	////////////////////////////////////////////////////////////////////////
	// Iteration tree
	////////////////////////////////////////////////////////////////////////

	ReverseConstantStepper stepper(
		0.,    // Initial time
		T,     // Expiry time
		T / N  // Timestep size
	);

	// Tolerance iteration
	unique_ptr<ToleranceIteration> toleranceIteration;
	if(method != EXPLICIT) {
		toleranceIteration = unique_ptr<ToleranceIteration>(
				new ToleranceIteration());
		stepper.setInnerIteration(*toleranceIteration);
	}

	////////////////////////////////////////////////////////////////////////
	// Linear system tree
	////////////////////////////////////////////////////////////////////////

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

		// Continuous withdrawal policy iteration
		continuousPolicy = unique_ptr<MinPolicyIteration2_1>(
			new MinPolicyIteration2_1(
				grid,
				*continuousControls,
				*continuousWithdrawal
			)
		);
		continuousPolicy->setIteration(*toleranceIteration);

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
		impulsePolicy->setIteration(*toleranceIteration);

		// Penalty method
		penalty = unique_ptr<PenaltyMethod>(
			new PenaltyMethod(
				grid,
				discretization,
				*impulsePolicy
			)
		);
		penalty->setIteration(*toleranceIteration);

		// Impulse root
		root = penalty.get();
	}

	////////////////////////////////////////////////////////////////////////
	// Exercise events
	////////////////////////////////////////////////////////////////////////

	auto withdrawal = [=] (const Interpolant2 &V, Real S, Real W) {

		// No withdrawal
		Real best = V(S, W);

		// Contract withdrawal amount
		const Real Gdt = G * T / N;

		if(method & SEMI_LAGRANGIAN_WITHDRAWAL_NO_PENALTY) {
			// Nonpenalty

			const Real beta = min(W, Gdt);

			//for(int i = 1; i <= M; ++i) {
				//const Real gamma = beta * i/M;
				const Real gamma = beta;

				const Real newValue =
					Vminus(V, S, W, gamma)
					+ gamma;

				if(newValue > best) {
					best = newValue;
				}
			//}

		}

		// Penalty
		if(
			(method & SEMI_LAGRANGIAN_WITHDRAWAL_AT_PENALTY)
			&& W > Gdt
		) {

			for(int i = 1; i <= M; ++i) {

				const Real gamma = Gdt + (W-Gdt) * i/M;

				const Real newValue =
					Vminus(V, S, W, gamma)
					+ gamma - kappa*(gamma-Gdt);

				if(newValue > best) {
					best = newValue;
				}

			}

		}

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

	////////////////////////////////////////////////////////////////////////
	// Payoff
	////////////////////////////////////////////////////////////////////////

	Function2 payoff = [=] (Real S, Real W) {
		return max(S, (1 - kappa) * W);
	};

	////////////////////////////////////////////////////////////////////////
	// Running
	////////////////////////////////////////////////////////////////////////

	BiCGSTABSolver solver;

	auto V = stepper.solve(
		grid,   // Domain
		payoff, // Initial condition
		*root,  // Root of linear system tree
		solver  // Linear system solver
	);

	Real mean = 1., var = 0.;
	int max = 1;

	if(method != EXPLICIT) {
		auto its = toleranceIteration->iterations();

		mean = accumulate(its.begin(),its.end(),0.)/its.size();

		var = 0.;
		for(auto x : its) { var += (x - mean) * (x - mean); }

		max = ( *max_element(its.begin(), its.end()) );
	}

	return make_tuple( V(w0, w0), mean, var, max );

}

constexpr int td = 20;

void printHeaders() {
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
}

int main() {

	////////////////////////////////////////////////////////////////////////
	// Table headers
	////////////////////////////////////////////////////////////////////////

	cout.precision(6);
	Real previousValue = nan(""), previousChange = nan("");

	if(!newton) {
		printHeaders();
	}

	////////////////////////////////////////////////////////////////////////
	// Refinement loop
	////////////////////////////////////////////////////////////////////////

	for(
		int l = 0;
		l <= Rmax;
		++l, N *= 2, M *= 2//, points1 *= 2, points2 *= 2, far *= 2.
	) {

		if( l < Rmin ) {
			// Refine grid
			grid.refine(RectilinearGrid2::NewTickBetweenEachPair());
			continue;
		}

		M = min(M, Mmax);

		////////////////////////////////////////////////////////////////
		// Outermost Newton iteration to find fair fee
		////////////////////////////////////////////////////////////////

		Real value, mean, var;
		int max;

		if(newton) {
			// Spacing
			cout << endl;

			// Newton iteration
			cout
				<< setw(td) << "\t"
				<< setw(td) << "Fair Fee"   << "\t"
				<< setw(td) << "Value"      << "\t"
				<< setw(td) << "Residual"   << "\t"
				<< endl
			;

			while(true) {
				// f(alpha)
				std::tie(value, mean, var, max) = solve(alpha);
				const Real f0 = value - w0;

				cout
					<< setw(td) << "\t"
					<< setw(td) << alpha << "\t"
					<< setw(td) << value << "\t"
					<< setw(td) << f0    << "\t"
					<< endl
				;

				if(abs(f0) < tolerance) {
					break;
				}

				// f(alpha + epsilon)
				auto tmp = solve(alpha + epsilon);
				const Real f1 = std::get<0>(tmp) - w0;

				// f'(alpha)
				const Real fp = (f1 - f0) / epsilon;

				// Next iterand
				alpha -= f1 / fp;
			}

			// Spacing
			cout << endl;

			// Print headers
			printHeaders();
		} else {
			// No Newton iteration
			std::tie(value, mean, var, max) = solve(alpha);
		}

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
			change = value - previousValue,
			ratio = previousChange / change
		;

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

		// Refine grid
		grid.refine(RectilinearGrid2::NewTickBetweenEachPair());

	}

	return 0;
}
