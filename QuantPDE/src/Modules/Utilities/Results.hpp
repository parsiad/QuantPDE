#ifndef QUANT_PDE_MODULES_RESULTS_HPP
#define QUANT_PDE_MODULES_RESULTS_HPP

#include <iomanip>  // std::setw
#include <iostream> // std::cout, std::endl, std::scientific
#include <cmath>    // std::abs, std::modf
#include <chrono>   // std::chrono
#include <memory>   // std::forward, std::unique_ptr
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

template <Index Dimension>
class ResultsBuffer {

	typedef std::vector<std::string> Headers;

	std::function<ResultsTuple<Dimension> (int k)> run;
	const Headers headers;
	const int kn; const int k0;
	const int precision; const int spacing;
	const bool ratio;
	std::unique_ptr<RectilinearGrid<Dimension>> grid;

public:

	/**
	 * @tparam Type of data.
	 * @param run Function to run.
	 * @param headers Table headers.
	 * @param kn Maximum level of refinement (inclusive).
	 * @param k0 Minimum level of refinement (inclusive).
	 * @param precision Level of numerical precision to display.
	 * @param spacing Amount of whitespace in table.
	 * @param ratio Whether or not to print the convergence ratio.
	 * @return A string.
	 */
	template <typename F>
	ResultsBuffer(
		F &&run,
		const Headers &headers = {},
		int kn = 5, int k0 = 0,
		int precision = 6, int spacing = 23,
		bool ratio = true
	) noexcept :
		run( std::forward<F>(run) ),
		headers(headers),
		kn(kn), k0(k0),
		precision(precision), spacing(spacing),
		ratio(ratio),
		grid(nullptr)
	{}

	// TODO: Might need to implement these in the future
	ResultsBuffer(const ResultsBuffer &that) = delete;
	ResultsBuffer &operator=(const ResultsBuffer &) = delete;

	ResultsBuffer &addPrintGrid(const RectilinearGrid<Dimension> &grid) {
		this->grid = std::unique_ptr<RectilinearGrid<Dimension>>(
				new RectilinearGrid<Dimension>(grid));
		return *this;
	}

	/**
	 * Runs and pushes output to the specified stream.
	 * @param os The output stream to print to.
	 */
	void stream(std::ostream &os = std::cout) {
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
		os << std::setw(spacing) << "Timing (Seconds)" << std::endl;

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

			// Reset format flags
			os.flags(f);

			// Print solution to grid
			if(k == kn && grid) {
				os << std::endl;
				for(int i = 0; i < Dimension; ++i) {
					os
						<< std::setw(spacing)
						<< ("x_" + std::to_string(i+1))
					;
				}
				os
					<< std::setw(spacing)
					<< "Solution at x"
					<< std::endl
					<< accessor( *grid, solution, spacing )
				;
			}
		}
	}

};

typedef ResultsBuffer<1> ResultsBuffer1;
typedef ResultsBuffer<2> ResultsBuffer2;
typedef ResultsBuffer<3> ResultsBuffer3;

}

#endif

