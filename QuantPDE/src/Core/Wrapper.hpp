#ifndef QUANT_PDE_CORE_WRAPPER
#define QUANT_PDE_CORE_WRAPPER

#define QUANT_PDE_CORE_WRAPPER_BODY(P, p)                                      \
		Wrapper(P p) noexcept : p( std::move(p) ) {                    \
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
		}                                                              \
		virtual P clone() const {                                      \
			return p->clone();                                     \
		}

#endif

