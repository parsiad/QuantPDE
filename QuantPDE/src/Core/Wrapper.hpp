#ifndef QUANT_PDE_CORE_WRAPPER
#define QUANT_PDE_CORE_WRAPPER

#define QUANT_PDE_CORE_WRAPPER_HEAD(CLASS_NAME)                                \
	class Cloneable : public CLASS_NAME {                                  \
	public:                                                                \
		typedef std::unique_ptr<Cloneable> P;                          \
		virtual P clone() const = 0;                                   \
	};                                                                     \
	class Wrapper : public CLASS_NAME {                                    \
		typename Cloneable::P p;                                       \
	public:                                                                \
		Wrapper(typename Cloneable::P p) noexcept : p( std::move(p) ) {\
		}                                                              \
		Wrapper(const Wrapper &that) noexcept : p(that.p->clone()) {   \
		}                                                              \
		Wrapper &operator=(const Wrapper &that) noexcept {             \
			p = that.p->clone();                                   \
			return *this;                                          \
		}                                                              \
		Wrapper(Wrapper &&that) noexcept : p( std::move(that.p) ) {    \
		}                                                              \
		Wrapper &operator=(Wrapper &&that) noexcept {                  \
			p = std::move(that.p);                                 \
			return *this;                                          \
		}

#define QUANT_PDE_CORE_WRAPPER_TAIL };

#endif

