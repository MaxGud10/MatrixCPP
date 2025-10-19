#pragma once

namespace Compare 
{
const double _epsilon = 1e-6;

bool is_equal           (const double x, const double y);
bool is_greater_or_equal(const double x, const double y);
bool is_less_or_equal   (const double x, const double y);

} // namespace Compare
