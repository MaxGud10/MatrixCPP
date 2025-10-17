#include <iostream>

#include "matrix.hpp"

int main() {
    size_t matrix_degree = 0;
    
    std::cin >> matrix_degree;
    if ((!std::cin.good()) || (matrix_degree <= 0)) {
        std::cerr << "Error input" << std::endl;
        return -1;
    }

    Matrix::matrix_t<double> matrix{matrix_degree, matrix_degree};

    std::cin >> matrix;

    std::cout << matrix.get_det_by_gauss_algorithm() << std::endl;

    return 0;
}