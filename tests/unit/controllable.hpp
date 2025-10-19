#ifndef CONTROLLABLE_HPP
#define CONTROLLABLE_HPP

#include <iostream>
#include <stdexcept>
#include <utility>

struct Controllable {
    static int control_;
    int* resource_;

    Controllable() : resource_(new int(42)) {}

    Controllable(Controllable &&rhs) noexcept : resource_(rhs.resource_) {
        rhs.resource_ = nullptr;
    }

    Controllable &operator=(Controllable &&rhs) noexcept {
        std::swap(resource_, rhs.resource_);
        return *this;
    }

    Controllable(const Controllable &rhs) : resource_(new int(*rhs.resource_)) {
        if (control_ == 0) {
            control_ = 1;
            std::cout << "Control reached" << std::endl;
            throw std::bad_alloc{};
        }
        control_ -= 1;
    }

    Controllable& operator=(const Controllable &rhs) {
        Controllable tmp(rhs);
        std::swap(*this, tmp);
        return *this;
    }

    ~Controllable() { 
        delete resource_; 
    }
};

#endif // CONTROLLABLE_HPP