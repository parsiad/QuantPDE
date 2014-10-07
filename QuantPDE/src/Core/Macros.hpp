#ifndef QUANT_PDE_CORE_MACROS_HPP
#define QUANT_PDE_CORE_MACROS_HPP

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

#ifdef __NVCC__
	// workaround issue between gcc >= 4.7 and cuda 5.5
	#if (defined __GNUC__) && (__GNUC__>4 || __GNUC_MINOR__>=7)
		#undef _GLIBCXX_ATOMIC_BUILTINS
		#undef _GLIBCXX_USE_INT128
	#endif

	#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int
#endif

#endif

