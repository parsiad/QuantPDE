#ifndef QUANT_PDE_CORE_RESULTS
#define QUANT_PDE_CORE_RESULTS

#include <iomanip>  // std::setw
#include <iostream> // std::cout, std::endl, std::scientific
#include <cmath>    // std::abs, std::modf
#include <vector>   // std::vector

namespace QuantPDE {

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
void results(
	std::function<std::vector<Real> (int k)> run,
	const std::vector<std::string> &headers = {},
	int kn = 5, int k0 = 0,
	std::ostream &os = std::cout,
	int precision = 6, int spacing = 23,
	bool ratio = true
) {
	// Store format flags to reset them later
	std::ios::fmtflags f(os.flags());

	// Headers
	for(auto it = headers.begin(); it != headers.end(); ++it) {
		os << std::setw(spacing) << *it;
	}
	os
		<< std::setw(spacing) << "Change"
		<< std::setw(spacing) << "Ratio"
		<< std::endl
	;

	os.precision(precision);

	Real previousValue = nan(""), previousChange = nan("");

	// Results
	for(int k = k0; k <= kn; ++k) {
		auto results = run(k);
		auto it = results.begin();
		Real value, fractpart, intpart;
		for(; it != results.end(); ++it) {
			value = *it;
			fractpart = modf(value, &intpart);
			if(std::abs(fractpart) < epsilon) {
				os << std::fixed;
			} else {
				os << std::scientific;
			}
			os << std::setw(spacing) << *it;
		}

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

		os << std::endl;
	}

	// Reset format flags
	os.flags(f);
}

}

#endif

