#ifndef QUANT_PDE_MODULES_RESULTS_HPP
#define QUANT_PDE_MODULES_RESULTS_HPP

#include <iomanip>  // std::setw
#include <iostream> // std::cout, std::endl, std::scientific
#include <cmath>    // std::abs, std::modf
#include <chrono>   // std::chrono
#include <tuple>    // std::tuple
#include <vector>   // std::vector

namespace QuantPDE {

template <Index Dimension>
using ResultsTuple = std::tuple<
	std::vector<Real>,
	InterpolantWrapper<Dimension>,
	Real
>;

typedef ResultsTuple<1> ResultsTuple1;
typedef ResultsTuple<2> ResultsTuple2;
typedef ResultsTuple<3> ResultsTuple3;

/**
 * Prints a table of results to a specified stream.
 * @tparam Type of data.
 * @param run Function to run.
 * @param headers Table headers.
 * @param kn Maximum level of refinement (inclusive).
 * @param k0 Minimum level of refinement (inclusive).
 * @param os The output stream to print to.
 * @param precision Level of numerical precision to display.
 * @param spacing Amount of whitespace in table.
 * @param ratio Whether or not to print the convergence ratio using the last
 *              column of the data.
 * @return A string.
 */
template <Index Dimension>
void streamResults(
	std::function<ResultsTuple<Dimension> (int k)> run,
	const std::vector<std::string> &headers,
	int kn, int k0,
	std::ostream &os,
	int precision, int spacing,
	bool ratio
) {
	// Store format flags to reset them later
	std::ios::fmtflags f(os.flags());

	// Headers
	for(auto it = headers.begin(); it != headers.end(); ++it) {
		os << std::setw(spacing) << *it;
	}
	os << std::setw(spacing) << "Value";
	if(ratio) {
		os
			<< std::setw(spacing) << "Change"
			<< std::setw(spacing) << "Ratio"
		;
	}
	os << std::setw(spacing) << "Timing (sec.)" << std::endl;

	os.precision(precision);

	Real previousValue = nan(""), previousChange = nan("");

	// Results
	for(int k = k0; k <= kn; ++k) {
		// Run
		Real seconds;
		auto start = std::chrono::steady_clock::now();
		auto ret = run(k); // Unpack
		auto end = std::chrono::steady_clock::now();
		auto diff = end - start;
		seconds = std::chrono::duration<Real>(diff).count();

		auto results = std::get<0>(ret);
		auto solution = std::get<1>(ret);
		const Real point = std::get<2>(ret);

		// Display results
		auto it = results.begin();
		for(; it != results.end(); ++it) {
			Real intpart;
			const Real fractpart = modf(*it, &intpart);
			if(std::abs(fractpart) < epsilon) {
				os << std::fixed;
			} else {
				os << std::scientific;
			}
			os << std::setw(spacing) << *it;
		}

		const Real value = solution(point);
		os << std::scientific;
		os << std::setw(spacing) << value;

		// Calculate and display ratio
		if(ratio) {
			const Real change = value - previousValue;
			const Real ratio = previousChange / change;

			os
				<< std::setw(spacing) << change
				<< std::setw(spacing) << ratio
			;

			previousChange = change;
			previousValue = value;
		}

		// Timing information
		os << std::setw(spacing) << seconds << std::endl;
	}

	// Reset format flags
	os.flags(f);
}

#define QUANT_PDE_MODULES_UTILITIES_RESULTS_FUNCTION(name, dimension) \
	void name( \
	std::function<ResultsTuple<dimension> (int k)> run, \
	const std::vector<std::string> &headers = {}, \
	int kn = 5, int k0 = 0, \
	std::ostream &os = std::cout, \
	int precision = 6, int spacing = 23, \
	bool ratio = true \
) { \
	streamResults<dimension>(run, headers, kn, k0, os, precision, spacing, \
		ratio); \
}

QUANT_PDE_MODULES_UTILITIES_RESULTS_FUNCTION(streamResults1, 1)
QUANT_PDE_MODULES_UTILITIES_RESULTS_FUNCTION(streamResults2, 2)
QUANT_PDE_MODULES_UTILITIES_RESULTS_FUNCTION(streamResults3, 3)

#undef QUANT_PDE_MODULES_UTILITIES_RESULTS_FUNCTION

}

#endif

