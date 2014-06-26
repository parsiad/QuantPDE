#ifndef QUANT_PDE_CORE_METAFUNCTIONS
#define QUANT_PDE_CORE_METAFUNCTIONS

#include <cstdint> // std::intmax_t
#include <utility> // std::forward

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

namespace Metafunctions {

template <int ...Is>
struct Sequence {
};

} // Metafunctions

/**
 * Used to generate a sequence of integers.
 */
template <int N, int ...Is>
struct GenerateSequence : GenerateSequence<N - 1, N - 1, Is...> {
};

template <int ...Is>
struct GenerateSequence<0, Is...> : Metafunctions::Sequence<Is...> {
};

////////////////////////////////////////////////////////////////////////////////

namespace Metafunctions {

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

} // Metafunctions

/**
 * Type definition for `R(T, ..., T)`, where `T` appears `N` times.
 */
template<class R, unsigned N, class T>
using NaryFunctionSignature = Metafunctions::NaryFunctionSignatureHelpers::Type<
		Metafunctions::NaryFunctionSignatureHelpers::Target, R, N, T>;

/**
 * Type definition for `R (C::*)(T, ..., T)` where `T` appears `N` times.
 */
template<class R, class C, unsigned N, class T>
using NaryMethodNonConst = Metafunctions::NaryMethodHelpers::Type<
		Metafunctions::NaryMethodHelpers::TargetNonConst, R, C, N, T>;

/**
 * Type definition for `R (C::*)(T, ..., T) const` where `T` appears `N` times.
 */
template<class R, class C, unsigned N, class T>
using NaryMethodConst = Metafunctions::NaryMethodHelpers::Type<
		Metafunctions::NaryMethodHelpers::TargetConst, R, C, N, T>;

////////////////////////////////////////////////////////////////////////////////

namespace Metafunctions {

namespace IsLValueHelpers {

template <typename T>
struct Nondeducible {
	typedef T type;
};

char (& Helper(...))[1];

template <typename T>
char (& Helper(T&, typename Nondeducible<const volatile T&>::type))[2];

} // IsLValueHelpers

} // Metafunctions

#define QUANT_PDE_IS_LVALUE(X) \
		(sizeof(QuantPDE::Metafunctions::IsLValueHelpers::Helper((X), \
		(X))) == 2)

#define QUANT_PDE_ASSERT_LVALUE(X) static_assert(QUANT_PDE_IS_LVALUE(X), \
		"Passing by rvalue is not supported")

////////////////////////////////////////////////////////////////////////////////

namespace Metafunctions {

/**
 * This class (and the corresponding function makeRRef) are used as a workaround
 * for the lack of generalized capture (providing move captures) in lambda
 * functions as of C++11.
 * \code{.cpp}
 * std::unique_ptr<int> p{new int(0)};
 * auto rref = makeRRef( std::move(p) );
 * auto lambda = [rref]() mutable {
 * 	return rref.get();
 * };
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

	RRef(T && x) noexcept : x{std::move(x)} {
	}

	RRef(RRef &that) noexcept : x{std::move(that.x)} {
		#ifndef NDEBUG
		copied = true;
		#endif
		assert(that.copied == false);
	}

	RRef(RRef && that) noexcept : x{std::move(that.x)} {
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

} // Metafunctions

template <typename T>
Metafunctions::RRef<T> makeRRef(T &&x) {
	return Metafunctions::RRef<T>{ std::move(x) };
}

} // QuantPDE

#endif

