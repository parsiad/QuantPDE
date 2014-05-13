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

/**
 * A sequence of integers.
 */
template <int... Is>
struct Sequence {
};

/**
 * Used to generate a sequence of integers.
 */
template <int N, int... Is>
struct GenerateSequence : GenerateSequence<N - 1, N - 1, Is...> {
};

template <int... Is>
struct GenerateSequence<0, Is...> : Sequence<Is...> {
};

/**
 * Calls a function by unpacking an array.
 */
template <typename F, typename T, Index... indices>
Real unpackAndCall(F &&function, const T *array, Sequence<indices...>) {
	return (std::forward<F>(function))( array[indices]... );
}

/*
template <std::size_t N, typename ReturnType, typename Type, typename ...Types>
struct NaryFunction
{
	using type = typename NaryFunction<N - 1, ReturnType, Type, Type,
			Types...>::type;
};

template <typename ReturnType, typename Type, typename ...Types>
struct NaryFunction<0, ReturnType, Type, Types...>
{
	using type = std::function <ReturnType (Types...)>;
};
*/

template<template<class...>class target, unsigned N, class R, class T,
		class... Ts>
struct RepeatTypeN: RepeatTypeN<target, N-1, R, T, T, Ts...> {};

template<template<class...>class target, class R, class T, class... Ts>
struct RepeatTypeN<target, 0, R, T, Ts...> {
	typedef target<R, Ts...> type;
};

template<template<class...>class target, unsigned N, class R, class T>
using RepeatTypeNTimes = typename RepeatTypeN<target, N, R, T>::type;

template <typename R, typename... Ts>
using Operation = R (Ts...);

template<unsigned N, class R, class T>
using NaryOperation = RepeatTypeNTimes<Operation, N, R, T>;

template<Index N>
using NRealToReal = NaryOperation<N, Real, Real>;

}

#endif

