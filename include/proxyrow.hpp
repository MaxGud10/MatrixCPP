#pragma once
#include <cstdlib> 

namespace matrix
{

namespace detail
{

template<class ElemT>
class ProxyRow
{
private:
    ElemT* row_;

public:
    ProxyRow(ElemT* row) : row_(row) {}

    const ElemT& operator[](size_t index) const { return row_[index]; }
          ElemT& operator[](size_t index)       { return row_[index]; }
};

} // namespace detail

} // namespace matrix
