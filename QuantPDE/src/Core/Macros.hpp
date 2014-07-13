#ifndef QUANT_PDE_CORE_TYPES_HPP
#define QUANT_PDE_CORE_TYPES_HPP

#include <cmath> // M_PI

// Some standard library implementations (VS) do not have M_PI
#ifndef M_PI
	#define 3.1415926535897932385
#endif

#if __cplusplus > 201103L
	// C++14
	#define QUANT_PDE_MOVE_CAPTURE(TYPE, ARGUMENT) \
			ARGUMENT(std::forward<TYPE>(ARGUMENT))
#else
	// C++11
	#define QUANT_PDE_MOVE_CAPTURE(TYPE, ARGUMENT) ARGUMENT
#endif

#endif

