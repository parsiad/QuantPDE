#ifndef QUANT_PDE_CORE_FUNCTION
#define QUANT_PDE_CORE_FUNCTION

#include <functional> // std::function

namespace QuantPDE {

/**
 * @see QuantPDE::NaryFunctionSignature
 */
template<Index N>
using NRealToReal = NaryFunctionSignature<Real, N, Real>;

typedef std::function< NRealToReal<1> > Function1;
typedef std::function< NRealToReal<2> > Function2;
typedef std::function< NRealToReal<3> > Function3;
typedef std::function< NRealToReal<4> > Function4;

template <Index N> // TODO: Check if nonnegative
using Function = std::function< NRealToReal<N> >;

////////////////////////////////////////////////////////////////////////////////

#define QUANT_PDE_TMP ( std::forward<F>(function) )( array[Indices]... )
template <typename F, typename T, int ...Indices>
inline auto packAndCall(F &&function, const T *array,
		Metafunctions::Sequence<Indices...>)
		-> decltype( QUANT_PDE_TMP ) {
	return QUANT_PDE_TMP;
}
#undef QUANT_PDE_TMP

/**
 * Call a function with each element of an array as an argument.
 * @tparam N Number of input arguments.
 * @param function The function.
 * @param array The array.
 */
#define QUANT_PDE_TMP packAndCall( std::forward<F>(function), array, \
		GenerateSequence<N>() )
template <int N, typename F, typename T>
inline auto packAndCall(F &&function, const T *array)
		-> decltype( QUANT_PDE_TMP ) {
	return QUANT_PDE_TMP;
}
#undef QUANT_PDE_TMP

////////////////////////////////////////////////////////////////////////////////

#define QUANT_PDE_TMP \
		( std::forward<F>(function) )( std::move(array[Indices])... )
template <typename F, typename T, int ...Indices>
inline auto packMoveAndCall(F &&function, T *array,
		Metafunctions::Sequence<Indices...>)
		-> decltype( QUANT_PDE_TMP ) {
	return QUANT_PDE_TMP;
}
#undef QUANT_PDE_TMP

/**
 * Call a function with each element of an array as an argument, invoking move
 * semantics on each element in the function call.
 * @tparam N Number of input arguments.
 * @param function The function.
 * @param array The array.
 */
#define QUANT_PDE_TMP packMoveAndCall( std::forward<F>(function), array, \
		GenerateSequence<N>() )
template <int N, typename F, typename T>
inline auto packMoveAndCall(F &&function, T *array)
		-> decltype( QUANT_PDE_TMP ) {
	return QUANT_PDE_TMP;
}
#undef QUANT_PDE_TMP

////////////////////////////////////////////////////////////////////////////////

#define QUANT_PDE_TMP ( std::forward<C>(caller).*method )( array[Indices]... )
template <typename C, typename M, typename T, int ...Indices>
inline auto packAndCall(C &&caller, M method, const T *array,
		Metafunctions::Sequence<Indices...>)
		-> decltype( QUANT_PDE_TMP ) {
	return QUANT_PDE_TMP;
}
#undef QUANT_PDE_TMP

/**
 * Call a method with each element of an array as an argument.
 * @tparam N Number of input arguments.
 * @param caller The object.
 * @param method The function.
 * @param array The array.
 */
#define QUANT_PDE_TMP packAndCall( std::forward<C>(caller), method, array, \
		GenerateSequence<N>() )
template <int N, typename C, typename M, typename T>
inline auto packAndCall(C &&caller, M method, const T *array)
		-> decltype( QUANT_PDE_TMP ) {
	return QUANT_PDE_TMP;
}
#undef QUANT_PDE_TMP

////////////////////////////////////////////////////////////////////////////////

#define QUANT_PDE_TMP \
		(std::forward<C>(caller).*method)(std::move(array[Indices])...)
template <typename C, typename M, typename T, int ...Indices>
inline auto packMoveAndCall(C &&caller, M method, T *array,
		Metafunctions::Sequence<Indices...>)
		-> decltype( QUANT_PDE_TMP ) {
	return QUANT_PDE_TMP;
}
#undef QUANT_PDE_TMP

/**
 * Call a method with each element of an array as an argument.
 * @tparam N Number of input arguments.
 * @param caller The object.
 * @param method The function.
 * @param array The array.
 */
#define QUANT_PDE_TMP packMoveAndCall( std::forward<C>(caller), method, array, \
		GenerateSequence<N>() )
template <int N, typename C, typename M, typename T>
inline auto packMoveAndCall(C &&caller, M method, T *array)
		-> decltype( QUANT_PDE_TMP ) {
	return QUANT_PDE_TMP;
}
#undef QUANT_PDE_TMP

}

#endif

