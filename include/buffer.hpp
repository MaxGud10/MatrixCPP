#pragma once
#include <algorithm>
#include <utility>

namespace matrix 
{


    
template <typename ElemT> 
void Constructor(ElemT* ptr, const ElemT& value) 
{
    std::construct_at(ptr, value);
    // new (ptr) ElemT{value};
}

template <typename ElemT> 
void Destructor(ElemT* ptr) 
{
    ptr->~ElemT();
}

template <typename Iterator> 
void Destructor(Iterator begin, Iterator end) 
{
    while (begin != end) 
    {
        Destructor(&(*begin));
        ++begin;
    }
}


template <typename ElemT> class Buffer 
{
protected:
    ElemT* buf_  = nullptr;
    size_t rows_ = 0;
    size_t cols_ = 0;
    size_t used_ = 0;

protected:
    Buffer(size_t rows, size_t cols): 
    buf_{static_cast<ElemT*>(::operator new(sizeof(ElemT) * rows * cols))}, rows_{rows}, cols_(cols), used_{0} {}

    
    Buffer& operator=(const Buffer& other_buf) = delete;
    Buffer(const Buffer& other_buf)            = delete;

    Buffer(Buffer&& other_buf) noexcept: buf_ {other_buf.buf_},
                                                           rows_{other_buf.rows_},
                                                           cols_{other_buf.cols_},
                                                           used_{other_buf.used_} 
    {
        other_buf.buf_  = nullptr;                                                           
        other_buf.rows_ = 0;
        other_buf.cols_ = 0;
        other_buf.used_ = 0;
    }

    Buffer& operator=(Buffer&& other_buf) 
    {
        if (this == &other_buf)
            return *this; 

        std::swap(buf_,  other_buf.buf_ );
        std::swap(rows_, other_buf.rows_);
        std::swap(cols_, other_buf.cols_);
        std::swap(used_, other_buf.used_);

        return *this;
    }

    ~Buffer() 
    {
        Destructor(buf_, buf_ + used_);
        ::operator delete(buf_);
    }
};

} // namespace matrix