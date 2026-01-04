#pragma once

#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <type_traits>
#include <vector>

#include "buffer.hpp"
#include "double_compare.hpp" 


namespace matrix 
{

template <typename ElemT>
class Matrix final : private Buffer<ElemT> 
{
    using Buffer<ElemT>::buf_ ;
    using Buffer<ElemT>::rows_;
    using Buffer<ElemT>::cols_;
    using Buffer<ElemT>::used_;

    struct ProxyRow 
    {
        ElemT *row;

        const ElemT& operator[](size_t n) const { return row[n]; }
              ElemT& operator[](size_t n)       { return row[n]; }
    };


    struct ConstProxyRow
    {
        const ElemT *row;

        const ElemT& operator[](size_t n) const { return row[n]; }
    };

    static constexpr long double zero_pivot_threshold_     = 1e-12L;
    static constexpr long double negligible_factor_thresh_ = 1e-18L;

    std::vector<long double> make_gauss_data_copy_() const
    {
        const size_t size = rows_;

        std::vector<long double> data(size * size);
        for (size_t idx = 0; idx < size * size; ++idx)
            data[idx] = static_cast<long double>(buf_[idx]);

        return data;
    }

    long double& gauss_element_(std::vector<long double>& data, size_t row, size_t col) const
    {
        const size_t size = rows_;
        return data[row * size + col];
    }

    size_t find_pivot_row_(std::vector<long double>& data, size_t pivot_col, long double& max_abs_in_col) const
    {
        const size_t size     = rows_;

        size_t pivot_row      = pivot_col;
        max_abs_in_col        = std::fabs(gauss_element_(data, pivot_col, pivot_col));

        for (size_t candidate_row = pivot_col + 1; candidate_row < size; ++candidate_row)
        {
            const long double candidate_abs = std::fabs(gauss_element_(data, candidate_row, pivot_col));
            if (candidate_abs > max_abs_in_col)
            {
                max_abs_in_col = candidate_abs;
                pivot_row      = candidate_row;
            }
        }

        return pivot_row;
    }

    void swap_gauss_rows_(std::vector<long double>& data, size_t first_row, size_t second_row) const
    {
        const size_t size = rows_;

        for (size_t col = 0; col < size; ++col)
            std::swap(gauss_element_(data, first_row, col), gauss_element_(data, second_row, col));
    }

    void eliminate_below_pivot_(std::vector<long double>& data, size_t pivot_col) const
    {
        const size_t size = rows_;

        const long double pivot_value = gauss_element_(data, pivot_col, pivot_col);

        for (size_t row = pivot_col + 1; row < size; ++row)
        {
            const long double elimination_factor = gauss_element_(data, row, pivot_col) / pivot_value;

            if (std::fabs(elimination_factor) <= negligible_factor_thresh_)
            {
                gauss_element_(data, row, pivot_col) = 0.0L;
                continue;
            }

            gauss_element_(data, row, pivot_col) = 0.0L;
            for (size_t col = pivot_col + 1; col < size; ++col)
                gauss_element_(data, row, col) -= elimination_factor * gauss_element_(data, pivot_col, col);
        }
    }

    long double multiply_gauss_diag_(std::vector<long double>& data) const
    {
        const size_t size = rows_;

        long double determinant = 1.0L;
        for (size_t d = 0; d < size; ++d)
            determinant *= gauss_element_(data, d, d);

        return determinant;
    }

    ElemT convert_gauss_det_(long double determinant, int row_swap_parity) const
    {
        if (row_swap_parity)
            determinant = -determinant;

        if constexpr (std::is_floating_point_v<ElemT>)
            return static_cast<ElemT>(determinant);

        return static_cast<ElemT>(std::llround(determinant));
    }

public:
    size_t get_used() const  { return used_; }
    size_t get_rows() const  { return rows_; }
    size_t get_cols() const  { return cols_; }

    ElemT multiply_diag() const 
    {
        assert(rows_ == cols_        ); 
        assert(buf_  != nullptr      );
        assert(used_ == rows_ * cols_);

        ElemT result = 1;

        for (size_t i = 0; i < rows_; ++i) 
            result *= buf_[i * cols_ + i];
        
        return result;
    }

    ElemT trace() const 
    {
        assert(rows_ == cols_        );
        assert(buf_  != nullptr      );
        assert(used_ == rows_ * cols_);

        ElemT result = 0;

        for (size_t i = 0; i < rows_; ++i) 
            result += buf_[i * cols_ + i];

        return result;
    }

    ProxyRow operator[](size_t n) 
    {
        assert(n     <  rows_        );
        assert(buf_  != nullptr      );
        assert(used_ == rows_ * cols_);

        return ProxyRow{buf_ + n * cols_};
    }

    ConstProxyRow operator[](size_t n) const
    {
        assert(n     <  rows_        );
        assert(buf_  != nullptr      );
        assert(used_ == rows_ * cols_);

        return ConstProxyRow{buf_ + n * cols_};
    }

    ElemT get_det_by_gauss_algorithm() const
    {
        assert(rows_ == cols_        );
        assert(buf_  != nullptr      );
        assert(used_ == rows_ * cols_);

        const size_t size = rows_;
        if (size == 0)
            return static_cast<ElemT>(1);

        std::vector<long double> data = make_gauss_data_copy_();

        int row_swap_parity = 0;

        for (size_t pivot_col = 0; pivot_col < size; ++pivot_col)
        {
            long double max_abs_in_col = 0.0L;
            const size_t pivot_row     = find_pivot_row_(data, pivot_col, max_abs_in_col);

            if (max_abs_in_col <= zero_pivot_threshold_)
                return static_cast<ElemT>(0);

            if (pivot_row != pivot_col)
            {
                swap_gauss_rows_(data, pivot_row, pivot_col);
                row_swap_parity ^= 1;
            }

            eliminate_below_pivot_(data, pivot_col);
        }

        const long double determinant = multiply_gauss_diag_(data);
        return convert_gauss_det_(determinant, row_swap_parity);
    }


    // ElemT get_det_by_gauss_algorithm() const
    // {
    //     assert(rows_ == cols_);

    //     const size_t size = rows_;
    //     if (size == 0)
    //         return static_cast<ElemT>(1);

    //     std::vector<long double> data(size * size);
    //     for (size_t idx = 0; idx < size * size; ++idx)
    //         data[idx] = static_cast<long double>(buf_[idx]);

    //     auto element = [&](size_t row, size_t col) -> long double& { return data[row * size + col]; };

        
    //     const long double zero_pivot_threshold     = 1e-12L;
    //     const long double negligible_factor_thresh = 1e-18L;

    //     int row_swap_parity = 0;

    //     for (size_t pivot_col = 0; pivot_col < size; ++pivot_col)
    //     {
    //         size_t      pivot_row       = pivot_col;
    //         long double max_abs_in_col  = std::fabs(element(pivot_col, pivot_col));

    //         for (size_t candidate_row = pivot_col + 1; candidate_row < size; ++candidate_row)
    //         {
    //             long double candidate_abs = std::fabs(element(candidate_row, pivot_col));
    //             if (candidate_abs > max_abs_in_col)
    //             {
    //                 max_abs_in_col = candidate_abs;
    //                 pivot_row      = candidate_row;
    //             }
    //         }

    //         if (max_abs_in_col <= zero_pivot_threshold)
    //             return static_cast<ElemT>(0);

    //         if (pivot_row != pivot_col)
    //         {
    //             for (size_t col = 0; col < size; ++col)
    //                 std::swap(element(pivot_col, col), element(pivot_row, col));
                    
    //             row_swap_parity ^= 1;
    //         }

    //         const long double pivot_value = element(pivot_col, pivot_col);
    //         for (size_t row = pivot_col + 1; row < size; ++row)
    //         {
    //             long double elimination_factor = element(row, pivot_col) / pivot_value;
    //             if (std::fabs(elimination_factor) <= negligible_factor_thresh)
    //             {
    //                 element(row, pivot_col) = 0.0L;
    //                 continue;
    //             }

    //             element(row, pivot_col) = 0.0L;
    //             for (size_t col = pivot_col + 1; col < size; ++col)
    //                 element(row, col) -= elimination_factor * element(pivot_col, col);
    //         }
    //     }

    //     long double determinant = 1.0L;
    //     for (size_t d = 0; d < size; ++d)
    //         determinant *= element(d, d);

    //     if (row_swap_parity)
    //         determinant = -determinant;

    //     if  (std::is_floating_point_v<ElemT>)
    //         return static_cast<ElemT>(determinant);
            
    //     else
    //         return static_cast<ElemT>(std::llround(determinant));
    // }

//===================================================================================================
public:
    void swap_rows(size_t first_row_number, size_t second_row_number) 
    {
        assert(first_row_number  < rows_);
        assert(second_row_number < rows_);
        assert(buf_ != nullptr          );

        if (first_row_number == second_row_number) 
            return;

        ElemT* first_row  = buf_ + first_row_number  * cols_;
        ElemT* second_row = buf_ + second_row_number * cols_;

        for (size_t i = 0; i < cols_; ++i) 
            std::swap(first_row[i], second_row[i]);
    }

    Matrix& negate() & 
    {
        assert(buf_  != nullptr      );
        assert(used_ == rows_ * cols_);

        for (size_t i = 0; i < rows_; ++i) 
        {
            for (size_t j = 0; j < cols_; ++j) 
            {
                buf_[i * cols_ + j] *= -1;
            }
        }

        return *this;
    }

    Matrix& transpose() & 
    {
        assert(buf_ != nullptr);
        assert(used_ == rows_ * cols_);

        Matrix transposed{cols_, rows_};
        for (size_t i = 0; i < rows_; ++i) 
        {
            for (size_t j = 0; j < cols_; ++j) 
            {
                transposed.buf_[j * rows_ + i] = buf_[i * cols_ + j];
            }
        }

        std::swap(buf_, transposed.buf_);
        std::swap(rows_, cols_);
        
        return *this;
    }

//===================================================================================================
public:
    explicit Matrix(size_t rows, size_t cols) : Buffer<ElemT>{rows, cols} 
    {
        while (used_ < rows_ * cols_) 
        {
            Constructor(buf_ + used_, ElemT{}); 
            ++used_;
        }
    }

    template <typename Iterator>
    Matrix(size_t rows, size_t cols, Iterator start, Iterator end): Buffer<ElemT>{rows, cols} 
    {
        for (auto iter = start; iter != end; ++iter) 
        {
            assert(used_ < rows_ * cols_);
            Constructor(buf_ + used_, static_cast<ElemT>(*iter));
            ++used_;
        }

        assert(used_ == rows_ * cols_);
    }

    template <typename AnotherElemT> 
    explicit Matrix(const Matrix<AnotherElemT>& other): Buffer<ElemT>{other.get_rows(), other.get_cols()} 
    {
        assert(rows_ == other.get_rows());
        assert(cols_ == other.get_cols());

        for (size_t i = 0; i < rows_; ++i)
        {
            for (size_t j = 0; j < cols_; ++j)
            {
                Constructor(buf_ + (i * cols_ + j),
                            static_cast<ElemT>(other[i][j]));
                ++used_;
            }
        }

        assert(used_ == rows_ * cols_);
    }
    
    Matrix(const Matrix& other): Buffer<ElemT>{other.get_rows(), other.get_cols()} 
    {
        assert(rows_ == other.get_rows());
        assert(cols_ == other.get_cols());

        for (size_t i = 0; i < rows_; ++i)
        {
            for (size_t j = 0; j < cols_; ++j)
            {
                Constructor(buf_ + (i * cols_ + j), other[i][j]);
                ++used_;
            }
        }

        assert(used_ == rows_ * cols_);
    }

    Matrix& operator=(const Matrix& other) 
    { 
        Matrix<ElemT> tmp{other};
        std::swap(*this, tmp);
        
        return *this;
    }

    Matrix(Matrix&& other)            noexcept = default;
    Matrix& operator=(Matrix&& other) noexcept = default;

    ~Matrix() = default;
};



//=================================================================================================
template <typename ElemT>
bool operator==(const Matrix<ElemT>& lhs, const Matrix<ElemT>& rhs) 
{
    if ((lhs.get_cols() != rhs.get_cols()) || (lhs.get_rows() != rhs.get_rows()))
        return false;

    for (size_t i = 0; i < lhs.get_rows(); ++i) 
    {
        for (size_t j = 0; j < lhs.get_cols(); ++j) 
        {
            if (!Compare::is_equal(lhs[i][j], rhs[i][j]))
                return false;
        }
    }

    return true;
}

template <typename ElemT>
bool operator!=(const Matrix<ElemT>& lhs, const Matrix<ElemT>& rhs) 
{
    return !(lhs == rhs);
}

template <typename ElemT> 
std::istream& operator>>(std::istream& inp_stream, Matrix<ElemT>& matrix) 
{
    size_t rows = matrix.get_rows();
    size_t cols = matrix.get_cols();

    for (size_t i = 0; i < rows; ++i) 
    {
        for (size_t j = 0; j < cols; ++j) 
        {
            inp_stream >> matrix[i][j];
        }
    }

    return inp_stream;
}

template <typename ElemT> 
std::ostream& operator<<(std::ostream& out_stream, const Matrix<ElemT>& matrix) 
{
    size_t rows = matrix.get_rows();
    size_t cols = matrix.get_cols();
    
    for (size_t i = 0; i < rows; ++i) 
    {
        for (size_t j = 0; j < cols; ++j) 
        {
            out_stream << matrix[i][j] << " ";
        }

        out_stream << std::endl;
    }

    out_stream << std::endl;

    return out_stream;
}

} // namespace matrix 