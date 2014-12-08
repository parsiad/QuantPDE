#ifndef QUANT_PDE_CORE_METAFUNCTIONS_HPP
#define QUANT_PDE_CORE_METAFUNCTIONS_HPP

#include <cstdint> // std::intmax_t
#include <cstdlib> // std::size_t
#include <utility> // std::forward

namespace QuantPDE {

// TODO: To intmax_t or not to intmax_t?

/**
 * Used to compute integer powers at compile-time. For example,
 * \code{.cpp}
 * // Compute 2^20 at compile-time
 * std::intmax_t ipow = IntegerPower<2, 20>::value;
 * \endcode
 * @tparam Base The base.
 * @tparam Exponent The exponent.
 */
template <std::intmax_t Base, std::intmax_t Exponent>
struct IntegerPower {
	typedef IntegerPower<Base, Exponent/2> HalfPower;

	static constexpr std::intmax_t temporary = HalfPower::value;

	/**
	 * The value of base to the exponent, assuming no overflow occurred.
	 * @see QuantPDE::IntegerPower::overflow
	 */
	static constexpr std::intmax_t value = temporary * temporary
			* (Exponent % 2 == 1 ? Base : 1);

	/**
	 * Logical true if and only if overflow occurred.
	 */
	static constexpr bool overflow = HalfPower::overflow ? true :
			(temporary > std::numeric_limits<intmax_t>::max() /
			(temporary * (Exponent % 2 == 1 ? Base : 1)) ? true
			: false);

	static_assert(!overflow, "Overflow detected");
};

/** @cond QUANT_PDE_HIDDEN */
template <std::intmax_t Base>
struct IntegerPower<Base, 0> {
	static constexpr std::intmax_t value = 1;
	static constexpr bool overflow = false;
};
/** @endcond */

////////////////////////////////////////////////////////////////////////////////

/**
 * Computes the integer product of factors. For example,
 * \code{.cpp}
 * // Compute the product 2 x 3 x 5  at compile-time
 * std::intmax_t iprod = IntegerProduct<2,3,5>::value;
 * \endcode
 * @tparam Factor The first factor.
 * @tparam Factors The remaining factors.
 */
template <std::intmax_t Factor, std::intmax_t ...Factors>
struct IntegerProduct {
	/**
	 * The product.
	 */
	static constexpr std::intmax_t value = IntegerProduct<Factors...>::value
			* Factor;
};

/** @cond QUANT_PDE_HIDDEN */
template <std::intmax_t Factor>
struct IntegerProduct<Factor> {
	static constexpr std::intmax_t value = Factor;
	//static constexpr bool overflow = false;
};
/** @endcond */

////////////////////////////////////////////////////////////////////////////////

/**
 * Selects the N-th type and value out of a parameter pack.
 *
 * \code{.cpp}
 * template <typename ...Ts>
 * void f(Ts ...args) {
 * 	Select<0, Ts...>::type firstType;
 * 	Select<1, Ts...>::type secondType;
 *
 * 	firstType arg1  = Select<0, Ts...>::get(args...);
 * 	secondType arg2 = Select<1, Ts...>::get(args...);
 * }
 * \endcode
 *
 * @tparam N A nonnegative integer.
 * @tparam Args The parameter pack.
 */
template <std::size_t N, typename ...Args>
struct Select {
	static_assert(N < sizeof...(Args), "Index out of bounds");
	typedef void type;
};

/** @cond QUANT_PDE_HIDDEN */
template <std::size_t N, typename T, typename ...Args>
struct Select<N, T, Args...> {
	typedef typename Select<N - 1, Args...>::type type;

	static inline const type &get(const T &, const Args &...args) {
		return Select<N - 1, Args...>::get(args...);
	}
};

template <typename T, typename ...Args>
struct Select<0, T, Args...> {
	typedef T type;

	static inline const type &get(const T &a, const Args &...) {
		return a;
	}
};
/** @endcond */

////////////////////////////////////////////////////////////////////////////////

/**
 * A sequence of integers.
 * @tparam Integers The sequence.
 */
template <int ...Integers>
struct Sequence {
};

/**
 * Generates the integers 0 to N - 1.
 * @tparam N The last integer (exclusive).
 * @see QuantPDE::Sequence
 */
template <int N, int ...Integers>
struct GenerateSequence : GenerateSequence<N - 1, N - 1, Integers...> {
};

/** @cond QUANT_PDE_HIDDEN */
template <int ...Integers>
struct GenerateSequence<0, Integers...> : Sequence<Integers...> {
};
/** @endcond */

////////////////////////////////////////////////////////////////////////////////

/** @cond QUANT_PDE_HIDDEN */
namespace NaryFunctionSignatureHelpers {

template <template <class...> class Target, class R, unsigned N, class T,
		class ...Ts>
struct Recurse : Recurse <Target, R, N - 1, T, T, Ts...> {
};

template <template <class...> class Target, class R, class T, class ...Ts>
struct Recurse<Target, R, 0, T, Ts...> {
	typedef Target<R, Ts...> type;
};

template <template <class...> class Target, class R, unsigned N, class T>
using Type = typename Recurse<Target, R, N, T>::type;

template <typename R, typename ...Ts>
using Target = R (Ts...);

} // NaryFunctionSignatureHelpers
/** @endcond */

/** @cond QUANT_PDE_HIDDEN */
namespace NaryMethodHelpers {

template <template <class...> class Target, class R, class C, unsigned N,
		class T, class ...Ts>
struct Recurse : Recurse <Target, R, C, N - 1, T, T, Ts...> {
};

template <template <class...> class Target, class R, class C, class T,
		class ...Ts>
struct Recurse<Target, R, C, 0, T, Ts...> {
	typedef Target<R, C, Ts...> type;
};

template <template <class...> class Target, class R, class C, unsigned N,
		class T>
using Type = typename Recurse<Target, R, C, N, T>::type;

template <typename R, class C, typename ...Ts>
using TargetNonConst = R (C::*)(Ts...);

template <typename R, class C, typename ...Ts>
using TargetConst = R (C::*)(Ts...) const;

} // NaryFunctionSignatureHelpers
/** @endcond */

/**
 * Type definition for `R(T, ..., T)`, where `T` appears `N` times.
 */
template<class R, unsigned N, class T>
using NaryFunctionSignature = NaryFunctionSignatureHelpers::Type<
	NaryFunctionSignatureHelpers::Target, R, N, T
>;

/**
 * Type definition for `R (C::*)(T, ..., T)` where `T` appears `N` times.
 */
template<class R, class C, unsigned N, class T>
using NaryMethodNonConst = NaryMethodHelpers::Type<
	NaryMethodHelpers::TargetNonConst, R, C, N, T
>;

/**
 * Type definition for `R (C::*)(T, ..., T) const` where `T` appears `N` times.
 */
template<class R, class C, unsigned N, class T>
using NaryMethodConst = NaryMethodHelpers::Type<
		NaryMethodHelpers::TargetConst, R, C, N, T
>;

////////////////////////////////////////////////////////////////////////////////

/** @cond QUANT_PDE_HIDDEN */
template <int N>
struct NaryFunctionPlaceholder {
	static NaryFunctionPlaceholder ph;
};

template<int N>
NaryFunctionPlaceholder<N> NaryFunctionPlaceholder<N>::ph;
/** @endcond */

/*
template <class R, class T, class ...Types, class U, int ...Indices>
std::function<R (Types...)> curry(std::function<R (T, Types...)> f, U val,
		Sequence<Indices...>) {
	return std::bind(f, val, NaryFunctionPlaceholder<Indices+1>::ph...);
}

template <class R, class T, class ...Types, class U>
std::function<R (Types...)> curry(std::function<R (T, Types...)> f, U val)
		{
	return curry(f, val, GenerateSequence<sizeof...(Types)>());
}
*/

#define QUANT_PDE_TMP ( std::bind( std::forward<F>(f), \
		value, NaryFunctionPlaceholder<Indices+1>::ph... ) )
template <int ...Indices, typename F, typename U>
auto curry(F &&f, U value, Sequence<Indices...>) -> decltype(QUANT_PDE_TMP) {
	return QUANT_PDE_TMP;
}
#undef QUANT_PDE_TMP

#define QUANT_PDE_TMP ( curry( std::forward<F>(f), value, \
		GenerateSequence<N-1>() ) )
/**
 * Binds the first value of a function-like object to a particular value.
 * @param f The function-like object.
 * @param value The value.
 * @tparam N The total number of arguments to the function-like object.
 */
template <int N, typename F, typename U>
auto curry(F &&f, U value) -> decltype(QUANT_PDE_TMP) {
	return QUANT_PDE_TMP;
}
#undef QUANT_PDE_TMP

////////////////////////////////////////////////////////////////////////////////

/** @cond QUANT_PDE_HIDDEN */
namespace IsLValueHelpers {

template <typename T>
struct Nondeducible {
	typedef T type;
};

char (& Helper(...))[1];

template <typename T>
char (& Helper(T&, typename Nondeducible<const volatile T&>::type))[2];

} // IsLValueHelpers
/** @endcond */

#define QUANT_PDE_IS_LVALUE(X) \
	(sizeof(QuantPDE::IsLValueHelpers::Helper((X), (X))) == 2)

////////////////////////////////////////////////////////////////////////////////

/**
 * This class (and the corresponding function makeRRef) are used as a workaround
 * for the lack of generalized capture (providing move captures) in lambda
 * functions as of C++11.
 * \code{.cpp}
 * std::unique_ptr<int> p{new int(0)};
 * auto rref = makeRRef( std::move(p) );
 * auto lambda = [rref]() mutable { return rref.get(); };
 * \endcode
 */
template <typename T>
class RRef {

	T x;

	#ifndef NDEBUG
	bool copied = false;
	#endif

public:

	RRef() = delete;

	RRef(T &&x) noexcept : x{std::move(x)} {
	}

	RRef(RRef &that) noexcept : x{std::move(that.x)} {
		#ifndef NDEBUG
		copied = true;
		#endif
		assert(that.copied == false);
	}

	RRef(RRef &&that) noexcept : x{std::move(that.x)} {
		#ifndef NDEBUG
		copied = that.copied;
		#endif
	}

	RRef &operator=(RRef that) = delete;

	T &&move() {
		return std::move(x);
	}

	T &operator*() const {
		return x;
	}

};

/**
 * @see QuantPDE::RRef
 */
template <typename T>
RRef<T> makeRRef(T &&x) {
	return RRef<T>{ std::move(x) };
}

} // QuantPDE

////////////////////////////////////////////////////////////////////////////////

/** @cond QUANT_PDE_HIDDEN */
// Note: adding members to the std namespace is not standards-compliant
namespace std {

template <int N>
struct is_placeholder< ::QuantPDE::NaryFunctionPlaceholder<N>> :
		std::integral_constant<int, N>
{};

} // std
/** @endcond */

#endif

