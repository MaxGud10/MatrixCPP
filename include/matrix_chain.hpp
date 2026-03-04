#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "matrix.hpp"

namespace matrix
{

template <typename T>
class MatrixChain
{
public:
    using MatrixType = Matrix<T>;

    MatrixChain() = default;

    void add_matrix(MatrixType matrix_to_add)
    {
        if (dimension_chain_.empty())
        {
            dimension_chain_.push_back(matrix_to_add.get_rows());
            dimension_chain_.push_back(matrix_to_add.get_cols());
        }
        else
        {
            const std::size_t expected_rows = dimension_chain_.back();
            if (expected_rows != matrix_to_add.get_rows())
                throw std::invalid_argument("MatrixChain: incompatible matrix dimensions in chain");

            dimension_chain_.push_back(matrix_to_add.get_cols());
        }

        matrix_storage_.push_back(std::move(matrix_to_add));
        dp_cache_is_dirty_ = true;
    }

    void add_dimensions(std::size_t rows_count, std::size_t cols_count)
    {
        if (dimension_chain_.empty())
        {
            dimension_chain_.push_back(rows_count);
            dimension_chain_.push_back(cols_count);
        }
        else
        {
            const std::size_t expected_rows = dimension_chain_.back();
            if (expected_rows != rows_count)
                throw std::invalid_argument("MatrixChain: incompatible dimensions in chain");

            dimension_chain_.push_back(cols_count);
        }

        dp_cache_is_dirty_ = true;
    }

    std::size_t get_matrix_count() const noexcept
    {
        return dimension_chain_.empty() ? 0 : (dimension_chain_.size() - 1);
    }

    std::uint64_t get_optimal_multiplication_cost()
    {
        build_dp_cache_if_needed_();

        const std::size_t matrix_count = get_matrix_count();
        if (matrix_count <= 1)
            return 0ULL;

        return dp_table_[0][matrix_count - 1].minimal_cost;
    }

    std::uint64_t get_naive_multiplication_cost() const
    {
        const std::size_t matrix_count = get_matrix_count();
        if (matrix_count <= 1)
            return 0ULL;

        std::uint64_t total_cost = 0;

        const std::size_t current_result_rows = dimension_chain_[0];
        std::size_t current_result_cols = dimension_chain_[1];

        for (std::size_t next_matrix_index = 1; next_matrix_index < matrix_count; ++next_matrix_index)
        {
            const std::size_t next_result_cols = dimension_chain_[next_matrix_index + 1];

            total_cost = add_multiplication_cost_checked_(
                total_cost,
                current_result_rows,
                current_result_cols,
                next_result_cols
            );

            current_result_cols = next_result_cols;
        }

        return total_cost;
    }

    std::vector<int> get_optimal_operation_order()
    {
        build_dp_cache_if_needed_();

        const std::size_t matrix_count = get_matrix_count();
        if (matrix_count <= 1)
            return {};

        return dp_table_[0][matrix_count - 1].operation_order;
    }

    double get_improvement_factor()
    {
        const std::uint64_t optimal_cost = get_optimal_multiplication_cost();
        const std::uint64_t naive_cost = get_naive_multiplication_cost();

        if (optimal_cost == 0)
            return 1.0;

        return static_cast<double>(naive_cost) / static_cast<double>(optimal_cost);
    }

    MatrixType multiply_in_optimal_order()
    {
        if (matrix_storage_.empty())
            throw std::logic_error("MatrixChain: no matrices stored; use add_matrix() to multiply");

        build_dp_cache_if_needed_();

        const std::size_t matrix_count = matrix_storage_.size();
        if (matrix_count == 0)
            throw std::logic_error("MatrixChain: empty chain");
        if (matrix_count == 1)
            return matrix_storage_[0]; // copy

        return multiply_range_recursively_(0, static_cast<int>(matrix_count - 1));
    }

private:
    struct DpCell
    {
        std::uint64_t minimal_cost = 0;
        int split_index = -1;
        std::vector<int> operation_order;
        bool is_initialized = false;
    };

    std::vector<std::size_t> dimension_chain_;
    std::vector<MatrixType>  matrix_storage_;

    bool dp_cache_is_dirty_ = true;
    std::vector<std::vector<DpCell>> dp_table_;

private:
    static std::uint64_t add_multiplication_cost_checked_(std::uint64_t base_cost,
                                                          std::size_t  left_rows,
                                                          std::size_t  shared_dimension,
                                                          std::size_t  right_cols)
    {
        __uint128_t multiplication_part =
            static_cast<__uint128_t>(left_rows) *
            static_cast<__uint128_t>(shared_dimension) *
            static_cast<__uint128_t>(right_cols);

        __uint128_t total_sum =
            static_cast<__uint128_t>(base_cost) + multiplication_part;

        if (total_sum > std::numeric_limits<std::uint64_t>::max())
            throw std::overflow_error("MatrixChain: multiplication cost overflow");

        return static_cast<std::uint64_t>(total_sum);
    }

    void build_dp_cache_if_needed_()
    {
        if (!dp_cache_is_dirty_)
            return;

        const std::size_t matrix_count = get_matrix_count();
        dp_table_.assign(matrix_count, std::vector<DpCell>(matrix_count));

        for (std::size_t single_index = 0; single_index < matrix_count; ++single_index)
        {
            dp_table_[single_index][single_index].minimal_cost   = 0;
            dp_table_[single_index][single_index].split_index    = -1;
            dp_table_[single_index][single_index].operation_order.clear();
            dp_table_[single_index][single_index].is_initialized = true;
        }

        for (std::size_t chain_length = 2; chain_length <= matrix_count; ++chain_length)
        {
            for (std::size_t left_index = 0; left_index + chain_length - 1 < matrix_count; ++left_index)
            {
                const std::size_t right_index = left_index + chain_length - 1;

                DpCell best_cell;
                best_cell.minimal_cost   = std::numeric_limits<std::uint64_t>::max();
                best_cell.is_initialized = false;

                for (std::size_t split_index = left_index; split_index < right_index; ++split_index)
                {
                    const DpCell& left_cell  = dp_table_[left_index     ][split_index];
                    const DpCell& right_cell = dp_table_[split_index + 1][right_index];

                    if (left_cell.minimal_cost > std::numeric_limits<std::uint64_t>::max() - right_cell.minimal_cost)
                        throw std::overflow_error("MatrixChain: cost overflow");

                    std::uint64_t candidate_cost = left_cell.minimal_cost + right_cell.minimal_cost;

                    candidate_cost = add_multiplication_cost_checked_(candidate_cost,
                                                                      dimension_chain_[left_index],
                                                                      dimension_chain_[split_index + 1],
                                                                      dimension_chain_[right_index + 1]);

                    std::vector<int> candidate_order;
                    candidate_order.reserve(left_cell .operation_order.size()  +
                                            right_cell.operation_order.size() + 1);

                    candidate_order.insert(candidate_order.end(),
                                           left_cell.operation_order.begin(),
                                           left_cell.operation_order.end());

                    candidate_order.insert(candidate_order.end(),
                                           right_cell.operation_order.begin(),
                                           right_cell.operation_order.end());

                    candidate_order.push_back(static_cast<int>(split_index));

                    const bool is_better_cost =
                        (!best_cell.is_initialized) || (candidate_cost < best_cell.minimal_cost);

                    const bool is_tie_but_lexicographically_smaller =
                        best_cell.is_initialized &&
                        (candidate_cost == best_cell.minimal_cost) &&
                        std::lexicographical_compare(
                            candidate_order.begin(),
                            candidate_order.end(),
                            best_cell.operation_order.begin(),
                            best_cell.operation_order.end());

                    if (is_better_cost || is_tie_but_lexicographically_smaller)
                    {
                        best_cell.minimal_cost    = candidate_cost;
                        best_cell.split_index     = static_cast<int>(split_index);
                        best_cell.operation_order = std::move(candidate_order);
                        best_cell.is_initialized  = true;
                    }
                }

                dp_table_[left_index][right_index] = std::move(best_cell);
            }
        }

        dp_cache_is_dirty_ = false;
    }

    MatrixType multiply_range_recursively_(int left_index, int right_index)
    {
        if (left_index == right_index)
            return matrix_storage_[static_cast<std::size_t>(left_index)]; // copy

        const int split_index =
            dp_table_[static_cast<std::size_t>(left_index)][static_cast<std::size_t>(right_index)].split_index;

        if (split_index < left_index || split_index >= right_index)
            throw std::logic_error("MatrixChain: invalid split");

        MatrixType left_product  = multiply_range_recursively_(left_index, split_index);
        MatrixType right_product = multiply_range_recursively_(split_index + 1, right_index);

        return left_product * right_product;
    }
};

} // namespace matrix