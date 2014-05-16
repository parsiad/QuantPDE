#ifndef QUANT_PDE_METAFUNCTIONS
#define QUANT_PDE_METAFUNCTIONS

#include <cstdint> // std::intmax_t

namespace QuantPDE {

// TODO: To intmax_t or not to intmax_t?

/**
 * Used to compute an integer powers as constant expressions.
 */
template <std::intmax_t Base, std::intmax_t Exponent>
struct IntegerPower {
	typedef IntegerPower<Base, Exponent/2> HalfPower;

	static constexpr std::intmax_t temporary = HalfPower::value;

	static constexpr std::intmax_t value = temporary * temporary
			* (Exponent % 2 == 1 ? Base : 1);

	static constexpr bool overflow = HalfPower::overflow ? true :
			(temporary > std::numeric_limits<intmax_t>::max() /
			(temporary * (Exponent % 2 == 1 ? Base : 1)) ? true
			: false);
};

template <std::intmax_t Base>
struct IntegerPower<Base, 0>
{
	static constexpr std::intmax_t value = 1;
	static constexpr bool overflow = false;
};

////////////////////////////////////////////////////////////////////////////////

namespace GenerateSequenceHelpers {

template <int... Is>
struct Sequence {
};

} // GenerateSequenceHelpers

/**
 * Used to generate a sequence of integers.
 */
template <int N, int... Is>
struct GenerateSequence : GenerateSequence<N - 1, N - 1, Is...> {
};

template <int... Is>
struct GenerateSequence<0, Is...> : GenerateSequenceHelpers::Sequence<Is...> {
};

////////////////////////////////////////////////////////////////////////////////

// Did not know where else to put this; not really a metafunction

namespace UnpackAndCall {

template <typename F, typename T, Index... indices>
Real call(F &&function, const T *array,
		GenerateSequenceHelpers::Sequence<indices...>) {
	return (std::forward<F>(function))( array[indices]... );
}

}

#define QUANT_PDE_UNPACK_AND_CALL(function, array, N)          \
		QuantPDE::UnpackAndCall::call(function, array, \
		GenerateSequence<N>())

////////////////////////////////////////////////////////////////////////////////

namespace NaryFunctionSignatureHelpers {

template <template <class...> class Target, unsigned N, class R, class T,
		class... Ts>
struct Recurse : Recurse <Target, N - 1, R, T, T, Ts...> {
};

template <template <class...> class Target, class R, class T, class... Ts>
struct Recurse<Target, 0, R, T, Ts...> {
	typedef Target<R, Ts...> type;
};

template <template <class...> class Target, unsigned N, class R, class T>
using Type = typename Recurse<Target, N, R, T>::type;

template <typename R, typename... Ts>
using Target = R (Ts...);

} // NaryFunctionSignatureHelpers

/**
 * NaryFunctionSignature is synonymous to the type R (T...), where T... is
 * repeated N times.
 */
template<unsigned N, class R, class T>
using NaryFunctionSignature = NaryFunctionSignatureHelpers::Type<
		NaryFunctionSignatureHelpers::Target, N, R, T>;

}

#endif

