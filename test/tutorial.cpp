////////////////////////////////////////////////////////////////////////////////
// tutorial.cpp
// ------------
//
// Code from the tutorial featured at
// https://github.com/parsiad/QuantPDE/wiki/Getting-Started
//
// Author: Parsiad Azimzadeh
////////////////////////////////////////////////////////////////////////////////

#include <QuantPDE/Core>
#include <QuantPDE/Modules/Operators>

using namespace QuantPDE;
using namespace QuantPDE::Modules;

#include <algorithm> // std::max
#include <iostream>  // std::cout, std::endl

using namespace std;

int main() {

// Strike price
const Real K = 100.;

// Automatic grid
/*
RectilinearGrid1 grid(
    Axis::cluster(
        0.,   // Left-hand boundary
        K,    // Feature to cluster nodes around
        200., // Right-hand boundary
        64,   // Number of nodes
        5.    // Intensity with which to cluster around K
    )
    + Axis { 225., 300., 750., 2000., 10000. }
);
*/

// Hand-picked grid (this is generally not a robust approach)
RectilinearGrid1 grid(
    Axis {
        0., 10., 20., 30., 40., 50., 60., 70.,
        75., 80.,
        84., 88., 92.,
        94., 96., 98., 100., 102., 104., 106., 108., 110.,
        114., 118.,
        123.,
        130., 140., 150.,
        175.,
        225.,
        300.,
        750.,
        2000.,
        10000.
    }
);

// Payoff for a put option
auto payoff = [K] (Real S) {
    return max( K - S, 0. );
};

// Initial time
const Real t0 = 0.;

// Expiry time
const Real T = 1.;

// Timestep size
const Real dt = 0.01;

// Number of exercise times
const int E = 10;

ReverseConstantStepper stepper(t0, T, dt);

// Discrete dividends
/*
// A discrete dividend of one dollar
const Real D = 1.;

// Relates the value of V- to V (see the above derivation)
auto dividendEvent = [D] (const Interpolant1 &V, Real S) {
    return V( max( S - D, 0. ) );
};

// Create uniformly spaced dividend payouts
for(int e = 0; e < E; ++e) {
    // Time at which the event takes place
    const Real te = T / E * e;

    // Add the event
    stepper.add(te, dividendEvent, grid);
}
*/

// Relates the value of V- to V
auto exerciseEvent = [K] (const Interpolant1 &V, Real S) {
    return max( V(S) , K - S );
};

// Create uniformly spaced exercise times
for(int e = 0; e < E; ++e) {
    // Time at which the event takes place
    const Real te = T / E * e;

    // Add the event
    stepper.add(te, exerciseEvent, grid);
}

// Interest rate
const Real r = 0.04;

// Local volatility function
const Real alpha = 20.;
auto v = [alpha] (Real S) {
    return alpha / S;
};

// A function of time and space could have been specified by using
// [] (Real t, Real S) { ... }

// Dividends
const Real q = 0.;

// Black-Scholes operator
BlackScholes1 bs(grid, r, v, q);

// Backward-differentiation formula of order two
ReverseBDFTwo1 bdf2(grid, bs);
bdf2.setIteration(stepper);

// Linear system solver
BiCGSTABSolver solver;

// Get the solution
auto V = stepper.solve(
    grid,   // Domain
    payoff, // Initial condition
    bdf2,   // Backward-differentiation formula
    solver  // Linear system solver
);

// Axis::range(a, b, c) creates ticks a to c using a step-size of b
// This is equivalent to the a:b:c notation used in MATLAB/Octave
RectilinearGrid1 printGrid( Axis::range(0., 10., 200.) );
cout << accessor( printGrid, V );

} // main()
