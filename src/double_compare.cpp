#include <limits>
#include <cmath>
#include <array>

#include "double_compare.hpp"

bool Compare::is_equal(const double x, const double y) {
    return std::fabs(x - y) < epsilon;
}

bool Compare::is_greater_or_equal(const double x, const double y) {
    return Compare::is_equal(x, y) || (x > y);
}

bool Compare::is_less_or_equal(const double x, const double y) {
    return Compare::is_equal(x, y) || (x < y);
}