#include <gtest/gtest.h>
#include <vector>
#include <array>
#include <iostream>

#include "matrix.hpp"
#include "double_compare.hpp"
#include "controllable.hpp"
#include "matrix_chain.hpp"

int Controllable::control_ = 0;

namespace {

template <typename T>
matrix::Matrix<T> make_matrix(std::size_t rows, std::size_t cols, std::initializer_list<T> values)
{
    std::vector<T> data(values);
    return matrix::Matrix<T>(rows, cols, data.begin(), data.end());
}

} // namespace

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

TEST(MATRIX_FUNCTIONS, copy_ctor_throws) {
    int old = Controllable::control_;
    Controllable::control_ = 1000;
    matrix::Matrix<Controllable> m1{3, 3};

    Controllable::control_ = 0;
    EXPECT_THROW( ([&]{ matrix::Matrix<Controllable> tmp{m1}; }()), std::bad_alloc );

    EXPECT_EQ(m1.get_rows(), 3u);
    EXPECT_EQ(m1.get_cols(), 3u);
    Controllable::control_ = old;
}

TEST(MATRIX_FUNCTIONS, copy_ctor_nothrow) {
    int old = Controllable::control_;
    Controllable::control_ = 1000;
    matrix::Matrix<Controllable> m1{3, 3};

    EXPECT_NO_THROW( ([&]{ matrix::Matrix<Controllable> tmp{m1}; }()) );
    Controllable::control_ = old;
}

TEST(MATRIX_FUNCTIONS, copy_assignment_throws) {
    int old = Controllable::control_;
    Controllable::control_ = 1000;
    matrix::Matrix<Controllable> m1{3, 3};
    matrix::Matrix<Controllable> m2{3, 3};

    Controllable::control_ = 0;
    EXPECT_THROW( (m2 = m1), std::bad_alloc );

    EXPECT_EQ(m2.get_rows(), 3u);
    EXPECT_EQ(m2.get_cols(), 3u);
    Controllable::control_ = old;
}

TEST(MATRIX_FUNCTIONS, copy_assignment_nothrow) {
    int old = Controllable::control_;
    Controllable::control_ = 1000;
    matrix::Matrix<Controllable> m1{3, 3};
    matrix::Matrix<Controllable> m2{3, 3};

    EXPECT_NO_THROW( (m2 = m1) );
    Controllable::control_ = old;
}

TEST(MATRIX_CONSTRUCTOR, iterator_range_too_short) {
    std::vector<int> data{};
    data.push_back(1);
    data.push_back(2);
    data.push_back(3);

    ASSERT_THROW(
        (matrix::Matrix<int>{2, 2, data.begin(), data.end()}),
        std::length_error
    );
}

TEST(MATRIX_CONSTRUCTOR, iterator_range_too_long) {
    std::vector<int> data{};
    data.push_back(1);
    data.push_back(2);
    data.push_back(3);
    data.push_back(4);
    data.push_back(5);

    ASSERT_THROW(
        (matrix::Matrix<int>{2, 2, data.begin(), data.end()}),
        std::length_error
    );
}


TEST(MATRIX_CONSTRUCTOR, iterator_range_exact_size) {
    std::vector<int> data{};
    data.push_back(1);
    data.push_back(2);
    data.push_back(3);
    data.push_back(4);

    ASSERT_NO_THROW(
        (matrix::Matrix<int>{2, 2, data.begin(), data.end()})
    );
}

TEST(MATRIX_CHAIN, add_dimensions_valid_chain)
{
    matrix::MatrixChain<double> chain;
    chain.add_dimensions(10, 30);
    chain.add_dimensions(30, 5);
    chain.add_dimensions(5, 60);

    EXPECT_EQ(chain.get_matrix_count(), 3u);
}

TEST(MATRIX_CHAIN, add_dimensions_invalid_chain_throws)
{
    matrix::MatrixChain<double> chain;
    chain.add_dimensions(10, 30);

    EXPECT_THROW(chain.add_dimensions(31, 5), std::invalid_argument);
}

TEST(MATRIX_CHAIN, naive_cost_example_from_text)
{
    matrix::MatrixChain<double> chain;
    chain.add_dimensions(10, 30);
    chain.add_dimensions(30, 5);
    chain.add_dimensions(5, 60);

    EXPECT_EQ(chain.get_naive_multiplication_cost(), 4500u);
    EXPECT_EQ(chain.get_optimal_multiplication_cost(), 4500u);
    EXPECT_DOUBLE_EQ(chain.get_improvement_factor(), 1.0);
}

TEST(MATRIX_CHAIN, prefers_better_parenthesization_cost)
{
    matrix::MatrixChain<double> chain;
    chain.add_dimensions(30, 35);
    chain.add_dimensions(35, 15);
    chain.add_dimensions(15, 5);
    chain.add_dimensions(5, 10);

    EXPECT_EQ(chain.get_naive_multiplication_cost(), 19500u);
    EXPECT_EQ(chain.get_optimal_multiplication_cost(), 9375u);

    const auto optimal_order = chain.get_optimal_operation_order();
    ASSERT_EQ(optimal_order.size(), 3u);
    EXPECT_EQ(optimal_order[0], 1);
    EXPECT_EQ(optimal_order[1], 0);
    EXPECT_EQ(optimal_order[2], 2);

    EXPECT_NEAR(chain.get_improvement_factor(), 2.08, 1e-12);
}

TEST(MATRIX_CHAIN, lexicographically_smallest_order_on_tie)
{
    matrix::MatrixChain<double> chain;
    chain.add_dimensions(10, 10);
    chain.add_dimensions(10, 10);
    chain.add_dimensions(10, 10);

    EXPECT_EQ(chain.get_optimal_multiplication_cost(), 2000u);

    const auto order = chain.get_optimal_operation_order();
    ASSERT_EQ(order.size(), 2u);
    EXPECT_EQ(order[0], 0);
    EXPECT_EQ(order[1], 1);
}

TEST(MATRIX_CHAIN, multiply_in_optimal_order_throws_if_only_dimensions_added)
{
    matrix::MatrixChain<int> chain;
    chain.add_dimensions(2, 2);
    chain.add_dimensions(2, 2);

    EXPECT_THROW(chain.multiply_in_optimal_order(), std::logic_error);
}

TEST(MATRIX_MULTIPLY, operator_multiply_basic)
{
    auto A = make_matrix<int>(2, 2, {1, 2, 3, 4});
    auto B = make_matrix<int>(2, 2, {5, 6, 7, 8});

    auto expected = make_matrix<int>(2, 2, {19, 22, 43, 50});
    auto result = A * B;

    EXPECT_TRUE(result == expected);
}

TEST(MATRIX_MULTIPLY, operator_multiply_incompatible_dimensions_throws)
{
    auto A = make_matrix<int>(2, 3, {1, 2, 3, 4, 5, 6});
    auto B = make_matrix<int>(2, 2, {1, 2, 3, 4});

    EXPECT_THROW((void)(A * B), std::invalid_argument);
}

TEST(MATRIX_CHAIN, multiply_in_optimal_order_correct_result_small)
{
    auto A = make_matrix<int>(2, 3, {
        1, 2, 3,
        4, 5, 6
    });

    auto B = make_matrix<int>(3, 2, {
        7,  8,
        9,  10,
        11, 12
    });

    auto C = make_matrix<int>(2, 2, {
        1, 0,
        0, 1
    });

    auto expected = make_matrix<int>(2, 2, {
        58, 64,
        139, 154
    });

    matrix::MatrixChain<int> chain;
    chain.add_matrix(A);
    chain.add_matrix(B);
    chain.add_matrix(C);

    auto result = chain.multiply_in_optimal_order();
    EXPECT_TRUE(result == expected);
}

TEST(MATRIX_CHAIN, dp_cache_updates_after_adding_more_matrices)
{
    matrix::MatrixChain<double> chain;
    chain.add_dimensions(10, 30);
    chain.add_dimensions(30, 5);

    EXPECT_EQ(chain.get_optimal_multiplication_cost(), 1500u);

    chain.add_dimensions(5, 60);
    EXPECT_EQ(chain.get_optimal_multiplication_cost(), 4500u);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}