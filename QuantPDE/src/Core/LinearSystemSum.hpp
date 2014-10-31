#ifndef QUANT_PDE_CORE_LINEAR_SYSTEM_SUM_HPP
#define QUANT_PDE_CORE_LINEAR_SYSTEM_SUM_HPP

#include <forward_list> // std::forward_list
#include <memory>       // std::unique_ptr

namespace QuantPDE {

class LinearSystemSum : public LinearSystem {

	typedef std::unique_ptr<LinearSystem> Ptr;
	std::forward_list<Ptr> systems;

	LinearSystemSum() noexcept {
	}

public:

	LinearSystemSum(LinearSystemSum &&sum) noexcept
			: systems(std::move(sum.systems)) {
	}

	virtual bool isATheSame() const {
		for(auto &system : systems) {
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

	friend LinearSystemSum operator+(Ptr a, Ptr b);
	friend LinearSystemSum &&operator+(LinearSystemSum &&O, Ptr a);
	friend LinearSystemSum &&operator+(Ptr a, LinearSystemSum &&O);

};

LinearSystemSum operator+(LinearSystemSum::Ptr a, LinearSystemSum::Ptr b) {
	LinearSystemSum O;
	O.systems.push_front(std::move(a));
	O.systems.push_front(std::move(b));
	return O;
}

LinearSystemSum &&operator+(LinearSystemSum &&O, LinearSystemSum::Ptr a) {
	O.systems.push_front(std::move(a));
	return std::move(O);
}

LinearSystemSum &&operator+(LinearSystemSum::Ptr a, LinearSystemSum &&O) {
	O.systems.push_front(std::move(a));
	return std::move(O);
}

}

#endif

