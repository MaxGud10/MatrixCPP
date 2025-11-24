#include <iostream>

#include "matrix.hpp"

int main() 
{
    size_t matrix_degree = 0;
    
    std::cin >> matrix_degree;

    const size_t MAX_MATRIX_SIZE = 2000; 

    if ((!std::cin.good()) || (matrix_degree == 0) || (matrix_degree > MAX_MATRIX_SIZE)) 
    {
        std::cerr << "Error input" << std::endl;
        return -1;
    }

    matrix::Matrix<double> matrix{matrix_degree, matrix_degree};

    std::cin >> matrix;

    if (!std::cin.good())  
    {
        std::cerr << "Error input" << std::endl;
        return -1;
    }

    std::cout << matrix.get_det_by_gauss_algorithm() << std::endl;

    return 0;
}