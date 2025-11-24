#pragma once

namespace Compare 
{
    
inline constexpr double epsilon = 1e-6;

template <typename T>
bool is_equal(const T& x, const T& y)
{
    if constexpr (std::is_floating_point_v<T>)
    {
        using std::fabs;
        return fabs(x - y) < epsilon;
    }
    else
    {
        return x == y;
    }
}

template <typename T>
bool is_greater_or_equal(const T& x, const T& y)
{
    return is_equal(x, y) || (x > y);
}

template <typename T>
bool is_less_or_equal(const T& x, const T& y)
{
    return is_equal(x, y) || (x < y);
}


} // namespace Compare
