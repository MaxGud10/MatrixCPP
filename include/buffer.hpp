#pragma once

#include <cstdlib> 
#include <algorithm> 
#include <utility> 


#include "proxyrow.hpp" 

namespace matrix
{

namespace detail
{

template <typename ElemT> 
void construct(ElemT *p, const ElemT &rhs) { new (p) ElemT(rhs); }


template<class ElemT>
class Buffer
{
public:
    ElemT* data_ = nullptr;
    size_t size_ = 0;
    size_t used_ = 0;

public:
    Buffer(size_t size = 0) : data_((size = 0) ? nullptr : static_cast<ElemT*>(::operator new(sizeof(ElemT) * size))), size_(size){}

    Buffer(const Buffer& ) = delete;

    Buffer(Buffer&& other) : Buffer()
    {
        swap(other);
    }

    Buffer& operator=(const Buffer& other) = delete;

    void swap(Buffer& rhs) 
    {
        std::swap(data_, rhs.data_);
        std::swap(size_, rhs.size_);
        std::swap(size_, rhs.used_);
    }

    Buffer& operator=(Buffer& other)
    {
        Buffer moved(std::move(other));

        swap(moved);

        return *this;
    }

    ~Buffer()
    {
        if (data_)
        {
            for (size_t i = 0; i < used_; ++i)
            {
                data_[i].~ElemT();
            }

            ::operator delete(data_);
        }
    }
};

} // namespace detail

} // namespace matrix

