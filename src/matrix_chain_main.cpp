#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "matrix_chain.hpp"

// TODO: переписать README
namespace // TODO: убрать
{

void print_order_line(const std::vector<int>& operation_order)
{
    for (std::size_t index = 0; index < operation_order.size(); ++index)
    {
        if (index != 0)
            std::cout << ' ';
        std::cout << operation_order[index];
    }
    std::cout << '\n';
}

void print_naive_order(std::size_t matrix_count)
{
    if (matrix_count <= 1)
    {
        std::cout << '\n';
        return;
    }

    for (std::size_t split_index = 0; split_index + 1 < matrix_count; ++split_index)
    {
        if (split_index != 0)
            std::cout << ' ';
        std::cout << split_index;
    }
    std::cout << '\n';
}
} // namespace

int main(int argc, char* argv[])
{
    const bool verbose_mode = (argc >= 2) && (std::string(argv[1]) == "--verbose");

    std::size_t dimension_count = 0;
    if (!(std::cin >> dimension_count))
    {
        std::cerr << "Error input\n";
        return -1;
    }

    if (dimension_count < 2)
    {
        if (verbose_mode)
        {
            std::cout << "optimal_order:\n";
            std::cout << "naive_order:\n";
            std::cout << "naive_cost: 0\n";
            std::cout << "optimal_cost: 0\n";
            std::cout << "factor: " << std::fixed << std::setprecision(6) << 1.0 << '\n';
        }
        else
            std::cout << std::fixed << std::setprecision(6) << 1.0 << '\n';

        return 0;
    }

    std::vector<std::size_t> dimension_values(dimension_count);
    for (std::size_t dimension_index = 0; dimension_index < dimension_count; ++dimension_index)
    {
        if (!(std::cin >> dimension_values[dimension_index]) || dimension_values[dimension_index] == 0)
        {
            std::cerr << "Error input\n";
            return -1;
        }
    }

    try
    {
        matrix::MatrixChain<double> matrix_chain;

        for (std::size_t matrix_index = 0; matrix_index + 1 < dimension_count; ++matrix_index)
        {
            const std::size_t rows_count = dimension_values[matrix_index];
            const std::size_t cols_count = dimension_values[matrix_index + 1];
            matrix_chain.add_dimensions(rows_count, cols_count);
        }

        const std::size_t matrix_count = matrix_chain.get_matrix_count();

        const std::vector<int> optimal_operation_order = matrix_chain.get_optimal_operation_order();

        const std::uint64_t    naive_cost   = matrix_chain.get_naive_multiplication_cost();
        const std::uint64_t    optimal_cost = matrix_chain.get_optimal_multiplication_cost();
        const double           factor       = (optimal_cost == 0) ? 1.0 : static_cast<double>(naive_cost) / static_cast<double>(optimal_cost);
        if (verbose_mode)
        {
            std::cout << "optimal_order: ";
            for (std::size_t index = 0; index < optimal_operation_order.size(); ++index)
            {
                if (index != 0)
                    std::cout << ' ';
                std::cout << optimal_operation_order[index];
            }
            std::cout << '\n';

            std::cout << "naive_order: ";
            if (matrix_count > 1)
            {
                for (std::size_t split_index = 0; split_index + 1 < matrix_count; ++split_index)
                {
                    if (split_index != 0)
                        std::cout << ' ';
                    std::cout << split_index;
                }
            }
            std::cout << '\n';

            std::cout << "naive_cost: "   << naive_cost << '\n';
            std::cout << "optimal_cost: " << optimal_cost << '\n';
            std::cout << "factor: "       << std::fixed << std::setprecision(6) << factor << '\n';
        }
        else
        {
            print_order_line(optimal_operation_order);
            std::cout << std::fixed << std::setprecision(6) << factor << '\n';
        }
    }
    catch (const std::exception&)
    {
        std::cerr << "Error input\n";
        return -1;
    }

    return 0;
}