#include <QuantPDE/Core>

#include "hjbqvi.hpp"

#include <cmath> // fabs

int main() {

	const RectilinearGrid1 spatial_grid = RectilinearGrid1(
		Axis::cluster(
			-6., // Left-hand boundary
			0.,  // Feature to cluster around
			+6., // Right-hand boundary
			32,  // Number of points
			10.  // Clustering intensity
		)
	);

	constexpr int dimension = 1; // Spatial dimension
	HJBQVI< dimension > hjbqvi(
		// Spatial grid
		spatial_grid,

		// Stochastic control grid
		RectilinearGrid1(
			Axis::uniform(
				0.,   // Left-hand boundary
				0.05, // Right-hand boundary
				32    // Number of points
			)
		),

		// Impulse control grid
		spatial_grid,

		// Expiry
		10., //numeric_limits<double>::infinity(),

		// Discount
		0.02,

		// Volatility
		[] (Real t, Real x) {
			const array<Real, 1> res { 0.3 };
			return res;
		},

		// Controlled drift
		[] (Real t, Real x, Real q) {
			const Real a = 0.25;
			const array<Real, 1> res { -a * q };
			return res;
		},

		// Uncontrolled drift
		[] (Real t, Real x) {
			const array<Real, 1> res { 0. };
			return res;
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
		[] (Real t, Real x, Real x_new) {
			const array<Real, 1> res { x_new };
			return res;
		},

		// Impulse flow
		[] (Real t, Real x, Real x_new) {
			const Real lambda = 1., c = 0.;
			const Real zeta = x_new - x;
			return - lambda * fabs(zeta) - c;
		},

		// Exit function
		[] (Real t, Real x) { return 0.; },

		// How to handle the control
		ControlMethod::FULLY_IMPLICIT,

		// Number of timesteps
		32
	);

	// Use common routine
	return hjbqvi_main(hjbqvi);

}

