#pragma once

#include <cstdlib> // size_t
#include <algorithm> // fill, copy
#include <utility> // move


#include "proxyrow.hpp" // ProxyRow

namespace matrix
{

namespace detail
{

template<class ElemT>
class Buffer
{
private:
    ElemT* data_;
    size_t size_;
    size_t cols_;

public:
    Buffer(size_t rows, size_t cols, const ElemT& value) : data_(new ElemT[rows * cols]), size_(rows * cols), cols_(cols)
    {
        std::fill(data_, data_ + size_, value);
    }
   
    Buffer(size_t dim, const ElemT& value) : data_(new ElemT[dim * dim]), size_(dim * dim), cols_(dim)
    {
        std::fill(data_, data_ + size_, value);
    }

    Buffer(const Buffer& other) : data_(new ElemT[other.size_]), size_(other.size_), cols_(other.cols_)
    {
        std::copy(other.data_, other.data_ + size_, data_);
    }

    Buffer& operator=(const Buffer& other)
    {
        if (this == &other)
            return *this;

        delete[] data_;
        
        size_ = other.size_;
        cols_ = other.cols_;

        std::copy(other.data_, other.data_ + size_, data_);
    }

    Buffer(Buffer&& other) noexcept : data_(std::move(other.data_)), size_(other.size_), cols_(other.cols_)
    {
        other.data_ = nullptr;
    }

    Buffer& operator=(Buffer&& other) noexcept
    {
        if (this == &other)
            return *this;

        delete[] data_;

        data_ = std::move(other.data_); 
        size_ = other.size_;
        cols_ = other.cols_;

        other.data_ = nullptr;
    }

    ~Buffer()
    {
        delete[] data_;
    }
    
    const ProxyRow<ElemT> operator[] (size_t index) const
    {
        return ProxyRow<ElemT>(data_ + index * cols_);
    }

    ProxyRow<ElemT> operator[] (size_t index)
    {
        return ProxyRow<ElemT>(data_ + index * cols_);
    }
};

} // namespace detail

} // namespace matrix

