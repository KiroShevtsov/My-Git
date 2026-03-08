#include <gtest/gtest.h>
#include "../Source/Slae.hxx"
#include <cmath>
TEST(VectorClass, VectorOperations){
    Vector vec1({1, 2});
    Vector vec2({3, 4});

    const double a = 2;
    const Vector result_mult_1({2, 4});
    const Vector result_mult_2({6, 8});

    /*сложение*/
    ASSERT_EQ(vec1 + vec2, Vector({4, 6}));
    ASSERT_EQ(vec2 + vec1, Vector({4, 6})) ;

    /*умножение на число*/
    ASSERT_EQ(vec1 * a, result_mult_1) ;
    ASSERT_EQ(vec2 * a, result_mult_2) ;

    EXPECT_EQ(vec1 * vec2, 11) ;
    EXPECT_EQ(vec2 * vec1, 11);

    EXPECT_THROW(Vector({}) * Vector({}), std::invalid_argument) ;
    EXPECT_THROW(vec1 + Vector({}), std::invalid_argument);
    EXPECT_THROW(vec1 + Vector({}) * (-1), std::invalid_argument);
    EXPECT_THROW(vec1 * Vector({}), std::invalid_argument);
}
TEST(DenseMatrixClass, MatrixOperation){
    ASSERT_THROW(Matrix(1, 2, {}), std::invalid_argument);
    ASSERT_THROW(Matrix(0, 0, {1, 3, 4}), std::invalid_argument);

    Matrix I(2, 2, {
        1, 0, 
        0, 1
    });
    Vector v({std::sqrt(2), std::sqrt(34)});
    ASSERT_EQ(I * v, v);
    ASSERT_EQ(I * (I * v), v);

    Matrix F(2, 2, {
        0, 1,
        1, 0
    });
    EXPECT_EQ(F * v, Vector({std::sqrt(34), std::sqrt(2)}));

    EXPECT_THROW(Matrix(0, 0, {}) * v, std::invalid_argument);
    Matrix M(3, 2, {
        1, 2, 3, 
        4, 5, 6
    });
    EXPECT_EQ(M(0, 0), 1);
    EXPECT_EQ(M(0, 1), 2);
    EXPECT_EQ(M(0, 2), 3);
    EXPECT_EQ(M(1, 0), 4);
    EXPECT_EQ(M(1, 1), 5);
    EXPECT_EQ(M(1, 2), 6);

    Matrix shift_down(3, 3, {
        0, 0, 1,
        1, 0, 0,
        0, 1, 0
    });
    
    Vector u({1, 2, 3});
    Vector expected({3, 1, 2});
    
    EXPECT_EQ(shift_down * u, expected);
    EXPECT_EQ(shift_down * (shift_down * (shift_down * u)), u);
}