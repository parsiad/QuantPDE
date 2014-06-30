#ifndef QUANT_PDE_CORE_WRAPPER
#define QUANT_PDE_CORE_WRAPPER

#define QUANT_PDE_CORE_WRAPPER_BODY(CLASS_NAME)                                \
		CLASS_NAME##Wrapper(P p) noexcept : p( std::move(p) ) {        \
		}                                                              \
		CLASS_NAME##Wrapper(const CLASS_NAME##Wrapper &that) noexcept  \
				: p(that.p->clone()) {                         \
		}                                                              \
		CLASS_NAME##Wrapper &operator=(const CLASS_NAME##Wrapper &that)\
				noexcept {                                     \
			p = that.p->clone();                                   \
			return *this;                                          \
		}                                                              \
		CLASS_NAME##Wrapper(CLASS_NAME##Wrapper &&that) noexcept       \
				: p( std::move(that.p) ) {                     \
		}                                                              \
		CLASS_NAME##Wrapper &operator=(CLASS_NAME##Wrapper &&that)     \
				noexcept {                                     \
			p = std::move(that.p);                                 \
			return *this;                                          \
		}                                                              \
		inline operator bool() const {                                 \
			return (bool)p;                                        \
		}                                                              \
		virtual P clone() const {                                      \
			return p->clone();                                     \
		}

#endif

