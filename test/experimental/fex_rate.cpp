////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>
#include "hjbqvi.hpp"

using namespace QuantPDE;
using namespace QuantPDE::HJBQVI;

////////////////////////////////////////////////////////////////////////////////

#include <array> // array
#include <cmath> // fabs

using namespace std;

////////////////////////////////////////////////////////////////////////////////

int main() {

	constexpr int Dimension = 1; // Spatial dimension

	const Axis X = Axis::cluster(
		-6., // Left-hand boundary
		0.,  // Feature to cluster around
		+6., // Right-hand boundary
		32,  // Number of points
		10.  // Clustering intensity
	);

	// Problem description
	Description<Dimension> hjbqvi(
		// Initial spatial grid
		{ X },

		// Initial stochastic control grid
		{ Axis::uniform( 0., 0.05, 32 ) },

		// Initial impulse control grid
		{ X },

		// Expiry
		10. /* numeric_limits<double>::infinity() */ ,

		// Discount
		0.02,

		// Volatility
		{
			[] (Real t, Real x) {
				return 0.3;
			}
		},

		// Controlled drift
		{
			[] (Real t, Real x, Real q) {
				const Real a = 0.25;
				return -a * q;
			}
		},

		// Uncontrolled drift
		{
			[] (Real t, Real x) {
				return 0.;
			}
		},

		// Controlled continuous flow
		[] (Real t, Real x, Real q) {
			const Real b = 3.;
			return -q * q * b;
		},

		// Uncontrolled continuous flow
		[] (Real t, Real x) {
			const Real m = 0.;
			const Real tmp = max(x - m, 0.);
			return -tmp * tmp;
		},

		// Transition
		{
			[] (Real t, Real x, Real x_new) {
				return x_new;
			}
		},

		// Impulse flow
		[] (Real t, Real x, Real x_new) {
			const Real lambda = 1., c = 0.;
			const Real zeta = x_new - x;
			return - lambda * fabs(zeta) - c;
		},

		// Exit function
		[] (Real t, Real x) {
			return 0.;
		},

		// How to handle the control
		ControlMethod::SEMI_LAGRANGIAN,

		// Initial number of timesteps
		32
	);

	Options<Dimension> opts(
		// Convergence test point
		{ 0. },

		// Max refinement
		6
	);

	// Use common routine
	run(hjbqvi, opts);

	return 0;

}

