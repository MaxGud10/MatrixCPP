#pragma once

#include <iostream>
#include <cmath>
#include <algorithm>
#include <type_traits>

#include "buffer.hpp"
#include "proxyrow.hpp"

namespace matrix 
{

template<typename ElemT>
class Matrix
{
private:
    detail::Buffer<ElemT> data_;
    size_t rows_;
    size_t cols_;

public:
    Matrix(size_t rows, size_t cols                    ) : rows_(rows), cols_(cols), data_(rows, cols, ElemT{}) {}
    Matrix(size_t rows, size_t cols, const ElemT& value) : rows_(rows), cols_(cols), data_(rows, cols, value)   {}
    Matrix(size_t dim                                  ) : rows_(dim), cols_(dim),   data_(dim, dim, ElemT{})   {}
    Matrix(size_t dim, const ElemT& value              ) : rows_(dim), cols_(dim),   data_(dim, dim, value)     {}

    template<typename FromT>
    Matrix(const Matrix<FromT>& other) : rows_(other.get_rows()), cols_(other.get_cols()), data_(other.get_rows(), other.get_cols())
    {
        for (size_t i = 0; i < rows_; i++)
        {
            for (size_t j = 0; j < cols_; j++)
            {
                data_[i][j] = static_cast<ElemT>(other[i][j]);
            }
        }
    }

    detail::ProxyRow<ElemT> operator[](size_t index)
    {
        return detail::ProxyRow<ElemT>(data_[index]);
    }

    const detail::ProxyRow<ElemT> operator[](size_t index) const
    {
        return detail::ProxyRow<ElemT>(data_[index]);
    }

    size_t get_rows() const { return rows_; }
    size_t get_cols() const { return cols_; }

    double det()     const
    {
        if (rows_ != cols_)
        {
            std::cerr << "can't count determinant!";

            return 0;
        }

        #ifdef DEBUG 
            dump();
        #endif

        size_t         swap_count  = 0;
        Matrix<double> matrix_copy = *this;
        
        if (!matrix_copy.gaussian_elimination(swap_count))
        {
            return 0;
        }

        #ifdef DEBUG
            matrix_copy.dump();
        #endif

        double det = 1;

        for (size_t i = 0; i < cols_; ++i)
        {
            det *= matrix_copy[i][i];
        }

        return (swap_count % 2 == 0) ? det : -det;
    }

    bool gaussian_elimination(size_t& swap_count)
    {
        for (size_t i_row = 0; i_row < rows_; ++i_row)
        {
            size_t pivot_row = find_pivot_row(i_row);

            if (data_[pivot_row][i_row] == 0)
            {
                return false;
            }

            if (pivot_row != i_row)
            {
                swap_rows(i_row, pivot_row);
                swap_count++;
            }

            eliminate_column(i_row);
            //#ifdef DEBUG
                dump();
            //#endif
        }

        return true; 
    }

    void dump() const
    {
        for (size_t row = 0; row < rows_; ++row)
        {
            for (size_t col = 0; col < cols_; ++col)
            {
                 std::cerr << data_[row][col] << " ";
            }
            
            std::cerr << "\n";
        }

        std::cerr << "\n";
    }

private:
    void swap_rows(size_t row1, size_t row2)
    {
        if (row1 == row2)
            return;

        for (size_t col = 0; col < cols_; ++col)
        {
            std::swap(data_[row1][col], data_[row2][col]);
        }
    }

    size_t find_pivot_row(size_t col)
    {
        size_t max_row  = col;
        ElemT  max_elem = std::abs(data_[col][col]);

        for (size_t i_col = col + 1; i_col < cols_; i_col++)
        {
            ElemT cur_elem = std::abs(data_[i_col][col]);

            if (cur_elem > max_elem)
            {
                max_elem = cur_elem; 
                max_row  = i_col;
            }
        }

        return max_row;
    }

    void eliminate_column(size_t pivot_row) 
    {
        for (size_t i_row = pivot_row + 1; i_row < rows_; ++i_row)
        {
            double factor = data_[i_row][pivot_row] / data_[pivot_row][pivot_row];

            for (int i_col = pivot_row; i_col < cols_; ++i_col) 
            {
                data_[i_row][i_col] -= factor * data_[pivot_row][i_col];
            }
        }
    }
};


} // namespace matrix 