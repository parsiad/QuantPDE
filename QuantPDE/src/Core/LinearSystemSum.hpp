#ifndef QUANT_PDE_CORE_LINEAR_SYSTEM_SUM_HPP
#define QUANT_PDE_CORE_LINEAR_SYSTEM_SUM_HPP

#include <forward_list> // std::forward_list

namespace QuantPDE {

class LinearSystemSum : public LinearSystem {

	std::forward_list<LinearSystem *> systems;

public:

	template <typename ...Ts>
	LinearSystemSum(Ts &...args) noexcept : systems( {(&args)...} ) {
	}

	virtual bool isATheSame() const {
		for(auto system : systems) {
			if(!system->isATheSame()) {
				return false;
			}
		}
		return true;
	}

	virtual Matrix A(Real time) {
		auto it = systems.begin();
		Matrix A = (*it++)->A(time);
		for(; it != systems.end(); ++it) {
			A += (*it)->A(time);
		}
		return A;
	}

	virtual Vector b(Real time) {
		auto it = systems.begin();
		Vector b = (*it++)->b(time);
		for(; it != systems.end(); ++it) {
			b += (*it)->b(time);
		}
		return b;
	}

};

}

#endif

