#pragma once

namespace Compare 
{    
const double epsilon = 1.0e-7;

bool is_equal(const double x, const double y);
bool is_greater_or_equal(const double x, const double y);
bool is_less_or_equal(const double x, const double y);

} // namespace Compare
