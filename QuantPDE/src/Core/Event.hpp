#ifndef QUANT_PDE_CORE_EVENT_HPP
#define QUANT_PDE_CORE_EVENT_HPP

#include <functional> // std::function
#include <memory>     // std::unique_ptr
#include <utility>    // std::forward, std::move

namespace QuantPDE {

namespace NaryFunctionSignatureHelpers {

template <typename R, typename ...Ts>
using TransformTarget = R (const Interpolant<sizeof...(Ts)> &, Ts...);

template <unsigned N>
using TransformSignature = Type<TransformTarget, Real, N, Real>;

} // Metafunctions

template <unsigned N>
using Transform = std::function<
		NaryFunctionSignatureHelpers::TransformSignature<N>>;

typedef Transform<1> Transform1;
typedef Transform<2> Transform2;
typedef Transform<3> Transform3;

/**
 * Transforms a vector.
 */
class EventBase {

	virtual Vector doEvent(const Vector &vector) const = 0;
	virtual Vector doEvent(Vector &&vector) const = 0;

public:

	template <typename V>
	Vector operator()(V &&vector) const {
		return doEvent( std::forward<V>(vector) );
	}

};

/**
 * An event describes a manner in which a solution on a domain is transformed.
 * In particular, an event creates an interpolant from the solution data,
 * transforms the interpolant using a lambda function, and maps the result back
 * onto a vector (whose components correspond to the domain nodes).
 *
 * For example, a Bermudan put event has the associated lambda function:
 * \code{.cpp}
 * [K] (const Function1 &V, Real S) {
 * 	// V is the solution prior to the event
 * 	// S is the coordinate (asset price)
 * 	max( V(S), K - S );
 * }
 * \endcode
 */
template <Index Dimension>
class Event : public EventBase {

private:

	typedef InterpolantFactoryWrapper<Dimension> In;
	typedef MapWrapper<Dimension> Out;

	Transform<Dimension> transform;

	In in;
	Out out;

	////////////////////////////////////////////////////////////////////////

	// TODO: Use C++1y move-capture in the future; return a bound function

	template <Index N, typename T, typename ...Ts>
	struct Driver : public Driver<N - 1, T, T, Ts...> {
	};

	template <typename T, typename ...Ts>
	struct Driver<0, T, Ts...> {
		static_assert( sizeof...(Ts) == Dimension,
				"The number of binding arguments in an event "
				"should be equal to the dimension");

		static inline Vector function(
			const Transform<Dimension> &transform,
			const Interpolant<Dimension> &interpolant,
			const Map<Dimension> &out
		) {
			return out([&] (Ts ...coordinates) {
				return transform(interpolant, coordinates...);
			});
		}
	};

	////////////////////////////////////////////////////////////////////////

	template <typename V>
	Vector _doEvent(V &&vector) const {
		// 1. Make interpolated function from vector
		// 2. Create a new the function by binding the first argument of
		//    the transform function to the interpolated function
		// 3. Map the function to the domain (producing a vector)

		return Driver<Dimension, Real>::function(
			transform,
			in.make( std::forward<V>(vector) ),
			out
		);
	}

	virtual Vector doEvent(const Vector &vector) const {
		return _doEvent(vector);
	}

	virtual Vector doEvent(Vector &&vector) const {
		return _doEvent(std::move(vector));
	}

public:

	/**
	 * Constructor.
	 * @param transform A function that transforms the solution.
	 * @param in An interpolant factory that handles how the solution is
	 *           interpolated off the domain nodes.
	 * @param out A map that handles how to transfer the transformed
	 *            solution back onto the domain.
	 */
	template <typename T>
	Event(
		T &&transform,
		In in,
		Out out
	) noexcept :
		transform(std::forward<T>(transform)),
		in(std::move(in)),
		out(std::move(out)) {
	}

	/**
	 * Default constructor for rectilinear grids.
	 * @param transform A function that transforms the solution.
	 * @param grid A rectilinear grid.
	 */
	template <typename T>
	Event(
		T &&transform,
		const RectilinearGrid<Dimension> &grid
	) noexcept :
		transform(std::forward<T>(transform)),
		in(grid.defaultInterpolantFactory()),
		out( std::unique_ptr<Map<Dimension>>(
				new PointwiseMap<Dimension>(grid)) )
	{
	}

	// Disable copy constructor and copy assignment operator.
	Event(const Event &) = delete;
	Event &operator=(const Event &) = delete;

	/**
	 * Move constructor.
	 */
	Event(Event &&that) noexcept :
		transform(std::move(that.transform)),
		in(std::move(that.in)),
		out(std::move(that.out))
	{
	}

	/**
	 * Move assignment operator.
	 */
	Event &operator=(Event &&that) & noexcept {
		transform = std::move(that.transform);
		in = std::move(that.in);
		out = std::move(that.out);
		return *this;
	}

};

typedef Event<1> Event1;
typedef Event<2> Event2;
typedef Event<3> Event3;

/**
 * The null event returns the original vector (no transformation occurs).
 */
class NullEvent : public EventBase {

	virtual Vector doEvent(const Vector &vector) const {
		return vector;
	}

	virtual Vector doEvent(Vector &&vector) const {
		return vector;
	}

};

/**
 * The do-event lets the user perform a specified action without modifying
 * the vector (no transformation occurs).
 */
class DoEvent : public EventBase {

	virtual void onCall(const Vector &vector) const = 0;

	virtual Vector doEvent(const Vector &vector) const {
		onCall(vector);
		return vector;
	}

	virtual Vector doEvent(Vector &&vector) const {
		onCall(vector);
		return vector;
	}

};

}

#endif
