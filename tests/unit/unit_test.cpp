#include <gtest/gtest.h>
#include <vector>
#include <array>
#include <iostream>

#include "matrix.hpp"
#include "double_compare.hpp"
#include "controllable.hpp"

int Controllable::control_ = 0;

TEST(MATRIX_FUNCTIONS, negate) {
    std::vector<double> vector_for_test{};
    vector_for_test.push_back(2.5);
    vector_for_test.push_back(7.9);
    vector_for_test.push_back(0.0);
    vector_for_test.push_back(-15.7);

    matrix::Matrix<double> matrix_for_test{2, 2, vector_for_test.begin(), vector_for_test.end()};

    std::vector<double> expected_vector{};
    expected_vector.push_back(-2.5);
    expected_vector.push_back(-7.9);
    expected_vector.push_back(-0.0);
    expected_vector.push_back(15.7);

    matrix::Matrix<double> expected_matrix{2, 2, expected_vector.begin(), expected_vector.end()};

    ASSERT_TRUE(matrix_for_test.negate() == expected_matrix);
} 

TEST(MATRIX_FUNCTIONS, transpose) {
    std::vector<double> vector_for_test{};
    vector_for_test.push_back(2.5);
    vector_for_test.push_back(7.9);
    vector_for_test.push_back(0.0);
    vector_for_test.push_back(-15.7);
    vector_for_test.push_back(100.0);
    vector_for_test.push_back(-90.0);    

    matrix::Matrix<double> matrix_for_test{2, 3, vector_for_test.begin(), vector_for_test.end()};

    std::vector<double> expected_vector{};
    expected_vector.push_back(2.5);
    expected_vector.push_back(-15.7);
    expected_vector.push_back(7.9);
    expected_vector.push_back(100.0);
    expected_vector.push_back(0.0);
    expected_vector.push_back(-90.0);

    matrix::Matrix<double> expected_matrix{3, 2, expected_vector.begin(), expected_vector.end()};

    ASSERT_TRUE(matrix_for_test.transpose() == expected_matrix);
} 

TEST(MATRIX_FUNCTIONS, swap_rows_1) {
    std::vector<double> vector_for_test{};
    vector_for_test.push_back(2.5);
    vector_for_test.push_back(7.9);
    vector_for_test.push_back(0.0);
    vector_for_test.push_back(-15.7);
    vector_for_test.push_back(100.0);
    vector_for_test.push_back(-90.0);    

    matrix::Matrix<double> matrix_for_test{2, 3, vector_for_test.begin(), vector_for_test.end()};

    std::vector<double> expected_vector{};
    expected_vector.push_back(-15.7);
    expected_vector.push_back(100.0);
    expected_vector.push_back(-90.0);       
    expected_vector.push_back(2.5);
    expected_vector.push_back(7.9);
    expected_vector.push_back(0.0);    
 
    matrix::Matrix<double> expected_matrix{2, 3, expected_vector.begin(), expected_vector.end()};

    matrix_for_test.swap_rows(0, 1);

    ASSERT_TRUE(matrix_for_test == expected_matrix);
} 

TEST(MATRIX_FUNCTIONS, swap_rows_2) {
    std::vector<double> vector_for_test{};
    vector_for_test.push_back(2.5);
    vector_for_test.push_back(7.9);
    vector_for_test.push_back(0.0);
    vector_for_test.push_back(-15.7);
    vector_for_test.push_back(100.0);
    vector_for_test.push_back(-90.0);

    matrix::Matrix<double> matrix_for_test{2, 3, vector_for_test.begin(), vector_for_test.end()};

    std::vector<double> expected_vector{};
    expected_vector.push_back(-15.7);
    expected_vector.push_back(100.0);
    expected_vector.push_back(-90.0);       
    expected_vector.push_back(2.5);
    expected_vector.push_back(7.9);
    expected_vector.push_back(0.0);
 
    matrix::Matrix<double> expected_matrix{2, 3, expected_vector.begin(), expected_vector.end()};

    matrix_for_test.swap_rows(1, 0);

    ASSERT_TRUE(matrix_for_test == expected_matrix);
}

TEST(MATRIX_FUNCTIONS, trace) {
    std::vector<double> vector_for_test{};
    vector_for_test.push_back(2.5);
    vector_for_test.push_back(7.9);
    vector_for_test.push_back(0.0);
    vector_for_test.push_back(-15.7);
    vector_for_test.push_back(100.0);
    vector_for_test.push_back(-90.0);
    vector_for_test.push_back(-100.7);
    vector_for_test.push_back(60.0);
    vector_for_test.push_back(-55.0);

    matrix::Matrix<double> matrix_for_test{3, 3, vector_for_test.begin(), vector_for_test.end()};

    ASSERT_DOUBLE_EQ(matrix_for_test.trace(), 47.5);
} 

TEST(MATRIX_FUNCTIONS, multiply_diag) {
    std::vector<double> vector_for_test{};
    vector_for_test.push_back(2.5);
    vector_for_test.push_back(7.9);
    vector_for_test.push_back(0.0);
    vector_for_test.push_back(-15.7);
    vector_for_test.push_back(100.0);
    vector_for_test.push_back(-90.0);
    vector_for_test.push_back(-100.7);
    vector_for_test.push_back(60.0);
    vector_for_test.push_back(-55.0);

    matrix::Matrix<double> matrix_for_test{3, 3, vector_for_test.begin(), vector_for_test.end()};

    ASSERT_DOUBLE_EQ(matrix_for_test.multiply_diag(), -13'750);
} 

TEST(MATRIX_FUNCTIONS, get_determinant_1) {
    std::vector<int> vector_for_test{};
    vector_for_test.push_back(1);
    vector_for_test.push_back(2);
    vector_for_test.push_back(3);
    vector_for_test.push_back(4);
    vector_for_test.push_back(5);
    vector_for_test.push_back(6);
    vector_for_test.push_back(7);
    vector_for_test.push_back(8);
    vector_for_test.push_back(9);

    matrix::Matrix<int> matrix_for_test{3, 3, vector_for_test.begin(), vector_for_test.end()};

    ASSERT_EQ(matrix_for_test.get_det_by_gauss_algorithm(), 0);
} 

TEST(MATRIX_FUNCTIONS, get_determinant_2) {
    std::vector<int> vector_for_test{};
    vector_for_test.push_back(2);
    vector_for_test.push_back(5);
    vector_for_test.push_back(-3);
    vector_for_test.push_back(1);
    vector_for_test.push_back(4);
    vector_for_test.push_back(-2);
    vector_for_test.push_back(-7);
    vector_for_test.push_back(3);
    vector_for_test.push_back(0);

    matrix::Matrix<int> matrix_for_test{3, 3, vector_for_test.begin(), vector_for_test.end()};

    ASSERT_EQ(matrix_for_test.get_det_by_gauss_algorithm(), -11);
} 

TEST(MATRIX_FUNCTIONS, get_determinant_3) {
    std::vector<int> vector_for_test{};
    vector_for_test.push_back(2);
    vector_for_test.push_back(0);
    vector_for_test.push_back(3);
    vector_for_test.push_back(0);
    vector_for_test.push_back(4);
    vector_for_test.push_back(0);
    vector_for_test.push_back(5);
    vector_for_test.push_back(0);
    vector_for_test.push_back(7);

    matrix::Matrix<int> matrix_for_test{3, 3, vector_for_test.begin(), vector_for_test.end()};

    ASSERT_EQ(matrix_for_test.get_det_by_gauss_algorithm(), -4);
} 

TEST(MATRIX_FUNCTIONS, get_determinant_4) {
    std::vector<int> vector_for_test{};
    vector_for_test.push_back(7);
    vector_for_test.push_back(-2);
    vector_for_test.push_back(9);
    vector_for_test.push_back(2);
    vector_for_test.push_back(1);
    vector_for_test.push_back(6);
    vector_for_test.push_back(-3);
    vector_for_test.push_back(-2);
    vector_for_test.push_back(7);

    matrix::Matrix<int> matrix_for_test{3, 3, vector_for_test.begin(), vector_for_test.end()};

    ASSERT_EQ(matrix_for_test.get_det_by_gauss_algorithm(), 188);
} 

TEST(MATRIX_FUNCTIONS, get_determinant_5) {
    std::vector<int> vector_for_test{};
    vector_for_test.push_back(1);
    vector_for_test.push_back(2);
    vector_for_test.push_back(1);
    vector_for_test.push_back(0);
    vector_for_test.push_back(-2);
    vector_for_test.push_back(3);
    vector_for_test.push_back(3);
    vector_for_test.push_back(1);
    vector_for_test.push_back(1);

    matrix::Matrix<int> matrix_for_test{3, 3, vector_for_test.begin(), vector_for_test.end()};
    ASSERT_EQ(matrix_for_test.get_det_by_gauss_algorithm(), 19);
} 

TEST(MATRIX_FUNCTIONS, get_determinant_6) {
    std::vector<int> vector_for_test{};
    vector_for_test.push_back(3);
    vector_for_test.push_back(-2);
    vector_for_test.push_back(1);
    vector_for_test.push_back(0);
    vector_for_test.push_back(2);
    vector_for_test.push_back(-1);
    vector_for_test.push_back(-3);
    vector_for_test.push_back(1);
    vector_for_test.push_back(0);

    matrix::Matrix<int> matrix_for_test{3, 3, vector_for_test.begin(), vector_for_test.end()};
    ASSERT_EQ(matrix_for_test.get_det_by_gauss_algorithm(), 3);
} 

TEST(MATRIX_FUNCTIONS, get_determinant_7) {
    std::vector<int> vector_for_test{};
    vector_for_test.push_back(-1);
    vector_for_test.push_back(2);
    vector_for_test.push_back(5);
    vector_for_test.push_back(7);
    vector_for_test.push_back(-4);
    vector_for_test.push_back(3);
    vector_for_test.push_back(-5);
    vector_for_test.push_back(0);
    vector_for_test.push_back(10);

    matrix::Matrix<int> matrix_for_test{3, 3, vector_for_test.begin(), vector_for_test.end()};

    ASSERT_EQ(matrix_for_test.get_det_by_gauss_algorithm(), -230);
} 

TEST(MATRIX_FUNCTIONS, get_determinant_8) {
    std::vector<int> vector_for_test{};
    vector_for_test.push_back(55);
    vector_for_test.push_back(10);
    vector_for_test.push_back(1);
    vector_for_test.push_back(1);
    vector_for_test.push_back(1);
    vector_for_test.push_back(2);
    vector_for_test.push_back(3);
    vector_for_test.push_back(2);
    vector_for_test.push_back(4);

    matrix::Matrix<int> matrix_for_test{3, 3, vector_for_test.begin(), vector_for_test.end()};

    ASSERT_EQ(matrix_for_test.get_det_by_gauss_algorithm(), 19);
} 

TEST(MATRIX_FUNCTIONS, get_determinant_9) {
    std::vector<int> vector_for_test{};
    vector_for_test.push_back(54);
    vector_for_test.push_back(2);
    vector_for_test.push_back(77);
    vector_for_test.push_back(89);
    vector_for_test.push_back(77);
    vector_for_test.push_back(86);
    vector_for_test.push_back(54);
    vector_for_test.push_back(54);
    vector_for_test.push_back(65);

    matrix::Matrix<int> matrix_for_test{3, 3, vector_for_test.begin(), vector_for_test.end()};

    ASSERT_EQ(matrix_for_test.get_det_by_gauss_algorithm(), 67'108);
} 

TEST(MATRIX_FUNCTIONS, get_determinant_10) {
    std::vector<int> vector_for_test{};
    vector_for_test.push_back(57);
    vector_for_test.push_back(4);
    vector_for_test.push_back(9);
    vector_for_test.push_back(1);
    vector_for_test.push_back(87);
    vector_for_test.push_back(54);
    vector_for_test.push_back(5);
    vector_for_test.push_back(0);
    vector_for_test.push_back(776);

    matrix::Matrix<int> matrix_for_test{3, 3, vector_for_test.begin(), vector_for_test.end()};

    ASSERT_EQ(matrix_for_test.get_det_by_gauss_algorithm(), 3'842'245);
} 

TEST(MATRIX_FUNCTIONS, get_determinant_11) {
    std::vector<int> vector_for_test{};
    vector_for_test.push_back(57);
    vector_for_test.push_back(4);
    vector_for_test.push_back(9);
    vector_for_test.push_back(1);
    vector_for_test.push_back(87);
    vector_for_test.push_back(54);
    vector_for_test.push_back(5);
    vector_for_test.push_back(0);
    vector_for_test.push_back(776);

    matrix::Matrix<int> matrix_for_test{3, 3, vector_for_test.begin(), vector_for_test.end()};

    ASSERT_EQ(matrix_for_test.get_det_by_gauss_algorithm(), 3'842'245);
} 

TEST(MATRIX_FUNCTIONS, get_determinant_12) {
    std::vector<int> vector_for_test{};
    vector_for_test.push_back(7);
    vector_for_test.push_back(1);
    vector_for_test.push_back(3);
    vector_for_test.push_back(333);

    vector_for_test.push_back(9);
    vector_for_test.push_back(99);
    vector_for_test.push_back(2);
    vector_for_test.push_back(3);

    vector_for_test.push_back(6);
    vector_for_test.push_back(1);
    vector_for_test.push_back(7);
    vector_for_test.push_back(132);

    vector_for_test.push_back(12);
    vector_for_test.push_back(4);
    vector_for_test.push_back(7);
    vector_for_test.push_back(4333);

    matrix::Matrix<int> matrix_for_test{4, 4, vector_for_test.begin(), vector_for_test.end()};

    ASSERT_EQ(matrix_for_test.get_det_by_gauss_algorithm(), 11'631'847);
}

TEST(MATRIX_FUNCTIONS, get_determinant_13) {
    std::vector<int> vector_for_test{};
    vector_for_test.push_back(777);
    vector_for_test.push_back(3);
    vector_for_test.push_back(6);
    vector_for_test.push_back(1);
    vector_for_test.push_back(4);
    vector_for_test.push_back(7);
    vector_for_test.push_back(2);
    vector_for_test.push_back(5);
    vector_for_test.push_back(8);

    matrix::Matrix<int> matrix_for_test{3, 3, vector_for_test.begin(), vector_for_test.end()};

    ASSERT_EQ(matrix_for_test.get_det_by_gauss_algorithm(), -2'331);
} 

TEST(MATRIX_FUNCTIONS, get_determinant_14) {
    std::vector<int> vector_for_test{};
    vector_for_test.push_back(1111);
    vector_for_test.push_back(123);
    vector_for_test.push_back(123);
    vector_for_test.push_back(7);
    vector_for_test.push_back(1);
    vector_for_test.push_back(77);
    vector_for_test.push_back(7);
    vector_for_test.push_back(9);
    vector_for_test.push_back(89);

    matrix::Matrix<int> matrix_for_test{3, 3, vector_for_test.begin(), vector_for_test.end()};

    ASSERT_EQ(matrix_for_test.get_det_by_gauss_algorithm(), -674'488);
} 

TEST(MATRIX_FUNCTIONS, get_determinant_15) {
    std::vector<double> vector_for_test{};
    vector_for_test.push_back(1.5);
    vector_for_test.push_back(67.0);
    vector_for_test.push_back(32.2);

    vector_for_test.push_back(-34.8);
    vector_for_test.push_back(43.0);
    vector_for_test.push_back(54.0);

    vector_for_test.push_back(0.5);
    vector_for_test.push_back(0.8);
    vector_for_test.push_back(0.7);

    matrix::Matrix<double> matrix_for_test{3, 3, vector_for_test.begin(), vector_for_test.end()};

    ASSERT_DOUBLE_EQ(matrix_for_test.get_det_by_gauss_algorithm(), 1'832.722);
} 

TEST(MATRIX_FUNCTIONS, get_determinant_16) {
    std::array<int, 100> array_for_test{
        12, 2, -170, 1742, -8, -2, -9, -10, -3, 5,
        8, 2, -17, 536, -2, 0, -2, -3, -1, 1,
        22, -8, -306, 1608, -14, 3, -5, -7, 1, 3,
        8, -4, -34, 0, -2, 2, 1, 0, 1, -1,
        20, -2, -153, 1139, -8, 3, -3, -6, 0, 2,
        58, -6, -459, 3283, -23, 10, -8, -17, 0, 6,
        -24, -12, -51, -1675, 3, 2, 7, 10, 5, -3,
        -38, 2, 493, -4154, 23, 1, 19, 23, 5, -11,
        -2, 0, 17, -335, 1, 1, 2, 2, 1, -1,
        6, -2, -68, 402, -3, 1, -1, -2, 0, 1,
    };

    matrix::Matrix<int> matrix_for_test{10, 10, array_for_test.begin(), array_for_test.end()};
    ASSERT_EQ(matrix_for_test.get_det_by_gauss_algorithm(), 4556);
}

TEST(MATRIX_FUNCTIONS, get_determinant_17) {
    std::array<double, 100> array_for_test{
        12, 2, -170, 1742, -8, -2, -9, -10, -3, 5,
        8, 2, -17, 536, -2, 0, -2, -3, -1, 1,
        22, -8, -306, 1608, -14, 3, -5, -7, 1, 3,
        8, -4, -34, 0, -2, 2, 1, 0, 1, -1,
        20, -2, -153, 1139, -8, 3, -3, -6, 0, 2,
        58, -6, -459, 3283, -23, 10, -8, -17, 0, 6,
        -24, -12, -51, -1675, 3, 2, 7, 10, 5, -3,
        -38, 2, 493, -4154, 23, 1, 19, 23, 5, -11,
        -2, 0, 17, -335, 1, 1, 2, 2, 1, -1,
        6, -2, -68, 402, -3, 1, -1, -2, 0, 1,
    };

    matrix::Matrix<double> matrix_for_test{10, 10, array_for_test.begin(), array_for_test.end()};
    ASSERT_TRUE(Compare::is_equal(matrix_for_test.get_det_by_gauss_algorithm(), 4556.0));
}

TEST(MATRIX_FUNCTIONS, copy_ctor) {
    matrix::Matrix<Controllable> matrix1{3, 3};

    Controllable::control_ = 0;
    
    bool exception_thrown = false;

    try {
        matrix::Matrix<Controllable> matrix2{matrix1};
    } catch (std::bad_alloc&) {
        exception_thrown = true;
    }

    ASSERT_FALSE(exception_thrown);
}

TEST(MATRIX_FUNCTIONS, copy_assignment) {
    matrix::Matrix<Controllable> matrix1{3, 3};

    Controllable::control_ = 0;

    bool exception_thrown = false;
    
    try {
        matrix::Matrix<Controllable> matrix2 = matrix1;
    } catch (std::bad_alloc&) {
        exception_thrown = true;
    }

    ASSERT_FALSE(exception_thrown); 
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}