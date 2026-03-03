#pragma once

#include <cstdint>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <limits>

#include "matrix.hpp"

namespace matrix
{

template <typename T>
class MatrixChain
{
public:
    using MatrixT = Matrix<T>;

    MatrixChain() = default;

    void add_matrix(MatrixT m)
    {
        if (matrices_.empty())
        {
            dims_.push_back(m.get_rows());
            dims_.push_back(m.get_cols());
        }

        else
        {
            if (dims_.back() != m.get_rows())
                throw std::invalid_argument("MatrixChain: incompatible matrix dimensions in chain");
            dims_.push_back(m.get_cols());
        }
        matrices_.push_back(std::move(m));
        dirty_ = true;
    }

    void add_dims(std::size_t rows, std::size_t cols)
    {
        if (dims_.empty())
        {
            dims_.push_back(rows);
            dims_.push_back(cols);
        }
        else
        {
            if (dims_.back() != rows)
                throw std::invalid_argument("MatrixChain: incompatible dimensions in chain");
            dims_.push_back(cols);
        }
        dirty_ = true;
    }

    std::size_t matrix_count() const noexcept
    {
        return dims_.size() > 0 ? dims_.size() - 1 : 0;
    }

    std::uint64_t optimal_cost()
    {
        ensure_dp_();
        const auto n = matrix_count();

        return (n <= 1) ? 0ULL : dp_[0][n - 1].cost;
    }

    std::uint64_t naive_cost() const
    {
        const auto n = matrix_count();
        if (n <= 1)
            return 0ULL;

        std::uint64_t cost  = 0;
        std::size_t   cur_r = dims_[0];
        std::size_t   cur_c = dims_[1];

        for (std::size_t i = 1; i < n; ++i)
        {
            const std::size_t next_c = dims_[i + 1];
            cost  = add_cost_(cost, cur_r, cur_c, next_c);
            cur_c = next_c;
        }

        return cost;
    }

    std::vector<int> optimal_order()
    {
        ensure_dp_();
        const auto n = matrix_count();
        if (n <= 1)
            return {};
        return dp_[0][n - 1].order;
    }

    double factor()
    {
        const auto opt = optimal_cost();
        const auto nai = naive_cost();
        if (opt == 0)
            return 1.0;
        return static_cast<double>(nai) / static_cast<double>(opt);
    }

    MatrixT multiply_optimal()
    {
        if (matrices_.empty())
            throw std::logic_error("MatrixChain: no matrices stored; use add_matrix() to multiply");

        ensure_dp_();
        const auto n = matrices_.size();
        if (n == 0)
            throw std::logic_error("MatrixChain: empty chain");
        if (n == 1)
            return matrices_[0];

        return multiply_range_(0, static_cast<int>(n - 1));
    }

private:
    struct Cell
    {
        std::uint64_t    cost = 0;
        int              split = -1;
        std::vector<int> order;
        bool             set = false;
    };

    std::vector<std::size_t> dims_;
    std::vector<MatrixT>     matrices_;
    bool                     dirty_ = true;

    std::vector<std::vector<Cell>> dp_;

    static std::uint64_t add_cost_(std::uint64_t base, std::size_t a, std::size_t b, std::size_t c)
    {
        __uint128_t mult = static_cast<__uint128_t>(a) * static_cast<__uint128_t>(b) * static_cast<__uint128_t>(c);
        __uint128_t sum  = static_cast<__uint128_t>(base) + mult;
        if (sum > std::numeric_limits<std::uint64_t>::max())
            throw std::overflow_error("MatrixChain: multiplication cost overflow");
        return static_cast<std::uint64_t>(sum);
    }

    void ensure_dp_()
    {
        if (!dirty_)
            return;

        const auto n = matrix_count();
        dp_.assign(n, std::vector<Cell>(n));

        for (std::size_t i = 0; i < n; ++i)
        {
            dp_[i][i].cost  = 0;
            dp_[i][i].split = -1;
            dp_[i][i].order.clear();
            dp_[i][i].set   = true;
        }

        for (std::size_t len = 2; len <= n; ++len)
        {
            for (std::size_t i = 0; i + len - 1 < n; ++i)
            {
                const std::size_t j = i + len - 1;

                Cell best;
                best.cost = std::numeric_limits<std::uint64_t>::max();
                best.set  = false;

                for (std::size_t k = i; k < j; ++k)
                {
                    const Cell& L = dp_[i][k];
                    const Cell& R = dp_[k + 1][j];

                    std::uint64_t cand = L.cost;
                    if (cand > std::numeric_limits<std::uint64_t>::max() - R.cost)
                        throw std::overflow_error("MatrixChain: cost overflow");
                    cand += R.cost;

                    cand = add_cost_(cand, dims_[i], dims_[k + 1], dims_[j + 1]);

                    std::vector<int> cand_order;
                    cand_order.reserve(L.order   .size() + R.order.size() + 1);
                    cand_order.insert (cand_order.end(),  L.order.begin(), L.order.end());
                    cand_order.insert (cand_order.end(),  R.order.begin(), R.order.end());

                    cand_order.push_back(static_cast<int>(k));

                    const bool better_cost = (!best.set) || (cand < best.cost);
                    const bool tie_lexi    =   best.set  && (cand == best.cost) &&
                        std::lexicographical_compare(cand_order.begin(), cand_order.end(),
                                                     best.order.begin(), best.order.end());

                    if (better_cost || tie_lexi)
                    {
                        best.cost  = cand;
                        best.split = static_cast<int>(k);
                        best.order = std::move(cand_order);
                        best.set   = true;
                    }
                }

                dp_[i][j] = std::move(best);
            }
        }

        dirty_ = false;
    }

    MatrixT multiply_range_(int i, int j)
    {
        if (i == j)
            return matrices_[static_cast<std::size_t>(i)];

        const int k = dp_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)].split;
        if (k < i || k >= j)
            throw std::logic_error("MatrixChain: invalid split");

        MatrixT left  = multiply_range_(i, k);
        MatrixT right = multiply_range_(k + 1, j);
        return left * right;
    }
};

} // namespace matrix
