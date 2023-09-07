#include <gtest/gtest.h>

#include "../s21_matrix_oop.h"

void CompareMatrix(S21Matrix matrix_1, S21Matrix matrix_2) {
  EXPECT_EQ(matrix_2.GetRows(), matrix_1.GetRows());
  EXPECT_EQ(matrix_2.GetCols(), matrix_1.GetCols());
  for (int i = 0; i < matrix_1.GetRows(); ++i) {
    for (int j = 0; j < matrix_2.GetCols(); ++j) {
      EXPECT_NEAR(matrix_1(i, j), matrix_2(i, j), 1E-5);
    }
  }
}

TEST(test_eq_matrix, test1) {
  S21Matrix matrix1(3, 3);
  matrix1(1, 1) = 2.000000010;
  S21Matrix matrix2(3, 3);
  matrix2(1, 1) = 2.000000011;
  EXPECT_EQ(matrix1.EqMatrix(matrix2), 1);
}

TEST(test_eq_matrix, test2) {
  S21Matrix matrix1(3, 3);
  S21Matrix matrix2(3, 3);
  matrix2(1, 1) = 2;
  EXPECT_EQ(matrix1.EqMatrix(matrix2), 0);
}

TEST(test_sum_matrix, test3) {
  S21Matrix matrix1(3, 3);
  matrix1(1, 1) = 2;
  S21Matrix matrix2(3, 3);
  matrix2(1, 1) = 2;
  S21Matrix matrix_result(3, 3);
  matrix_result(1, 1) = 4;
  matrix1.SumMatrix(matrix2);
  CompareMatrix(matrix1, matrix_result);
}

TEST(test_sum_matrix, test4) {
  S21Matrix matrix1(3, 3);
  S21Matrix matrix2;
  EXPECT_THROW(matrix1.SumMatrix(matrix2), std::out_of_range);
}

TEST(test_sub_matrix, test5) {
  S21Matrix matrix1(3, 3);
  S21Matrix matrix2(3, 3);
  S21Matrix matrix_result(3, 3);
  matrix1.SubMatrix(matrix2);
  CompareMatrix(matrix1, matrix_result);
}

TEST(test_sub_matrix, test6) {
  S21Matrix matrix1(3, 3);
  S21Matrix matrix2(2, 2);
  EXPECT_THROW(matrix1.SubMatrix(matrix2), std::out_of_range);
}

TEST(test_mul_number, test7) {
  S21Matrix matrix1(3, 3);
  S21Matrix matrix2(3, 3);
  double d = 3.0;
  matrix1.MulNumber(d);
  CompareMatrix(matrix1, matrix2);
}

TEST(test_mul_matrix, test8) {
  S21Matrix matrix1(3, 3);
  S21Matrix matrix2(3, 3);
  S21Matrix matrix_result(3, 3);
  matrix1.MulMatrix(matrix2);
  CompareMatrix(matrix1, matrix_result);
}

TEST(test_mul_matrix, test9) {
  S21Matrix matrix1(3, 3);
  S21Matrix matrix2(2, 2);
  EXPECT_THROW(matrix1.MulMatrix(matrix2), std::out_of_range);
}

TEST(test_transpose, test10) {
  S21Matrix matrix1(3, 3);
  matrix1(0, 1) = 2;
  S21Matrix matrix_result(3, 3);
  matrix_result(1, 0) = 2;
  CompareMatrix(matrix1.Transpose(), matrix_result);
}

TEST(test_calc_complements, test11) {
  S21Matrix matrix1(3, 3);
  matrix1(0, 0) = 1;
  matrix1(1, 1) = 2;
  matrix1(2, 2) = 3;
  S21Matrix matrix_result(3, 3);
  matrix_result(0, 0) = 6;
  matrix_result(1, 1) = 3;
  matrix_result(2, 2) = 2;
  CompareMatrix(matrix1.CalcComplements(), matrix_result);
}

TEST(test_calc_complements, test12) {
  S21Matrix matrix1(3, 2);
  EXPECT_THROW(matrix1.CalcComplements(), std::out_of_range);
}

TEST(test_determinant, test13) {
  S21Matrix matrix1(3, 3);
  matrix1(0, 0) = 1;
  matrix1(1, 1) = 2;
  matrix1(2, 2) = 3;
  double det = 6.0;
  EXPECT_EQ(matrix1.Determinant(), det);
}

TEST(test_determinant, test14) {
  S21Matrix matrix1(3, 2);
  EXPECT_THROW(matrix1.Determinant(), std::out_of_range);
}

TEST(test_inverse_matrix, test15) {
  S21Matrix matrix1(3, 3);
  matrix1(0, 0) = 1;
  matrix1(1, 1) = 2;
  matrix1(2, 2) = 3;
  S21Matrix matrix_result(3, 3);
  matrix_result(0, 0) = 1;
  matrix_result(1, 1) = 1.0 / 2.0;
  matrix_result(2, 2) = 1.0 / 3.0;
  CompareMatrix(matrix1.InverseMatrix(), matrix_result);
}

TEST(test_inverse_matrix, test16) {
  S21Matrix matrix1(1, 1);
  matrix1(0, 0) = 2;
  S21Matrix matrix_result(1, 1);
  matrix_result(0, 0) = 0.5;
  CompareMatrix(matrix1.InverseMatrix(), matrix_result);
}

TEST(test_inverse_matrix, test17) {
  S21Matrix matrix1(3, 2);
  EXPECT_THROW(matrix1.InverseMatrix(), std::out_of_range);
}

TEST(test_inverse_matrix, test18) {
  S21Matrix matrix1(3, 3);
  EXPECT_THROW(matrix1.InverseMatrix(), std::out_of_range);
}

TEST(test_operator_brackets, test19) {
  S21Matrix matrix1(3, 3);
  EXPECT_THROW(matrix1(4, 4), std::out_of_range);
}

TEST(test_operator_brackets, test20) {
  const S21Matrix matrix1(3, 3);
  double i = matrix1(0, 0);
  EXPECT_EQ(i, matrix1(0, 0));
}

TEST(test_operator_brackets, test21) {
  const S21Matrix matrix1(3, 3);
  EXPECT_THROW(matrix1(4, 4), std::out_of_range);
}

TEST(test_operator_add, test22) {
  S21Matrix matrix1(3, 3);
  matrix1(1, 1) = 2;
  S21Matrix matrix2(3, 3);
  matrix2(1, 1) = 2;
  S21Matrix matrix_sum(3, 3);
  S21Matrix matrix_result(3, 3);
  matrix_result(1, 1) = 4;
  matrix_sum = matrix1 + matrix2;
  CompareMatrix(matrix_sum, matrix_result);
}

TEST(test_operator_sub, test23) {
  S21Matrix matrix1(3, 3);
  matrix1(1, 1) = 1;
  S21Matrix matrix2(3, 3);
  matrix2(1, 1) = 1;
  S21Matrix matrix_sub(3, 3);
  S21Matrix matrix_result(3, 3);
  matrix_sub = matrix2 - matrix1;
  CompareMatrix(matrix_sub, matrix_result);
}

TEST(test_operator_mul, test24) {
  S21Matrix matrix1(3, 3);
  matrix1(1, 1) = 2;
  S21Matrix matrix2(3, 3);
  matrix2(1, 1) = 2;
  S21Matrix matrix_mul(3, 3);
  S21Matrix matrix_result(3, 3);
  matrix_result(1, 1) = 4;
  matrix_mul = matrix2 * matrix1;
  CompareMatrix(matrix_mul, matrix_result);
}

TEST(test_operator_mul, test25) {
  S21Matrix matrix1(3, 3);
  matrix1(1, 1) = 2;
  double m = 3;
  S21Matrix matrix_mul(3, 3);
  S21Matrix matrix_result(3, 3);
  matrix_result(1, 1) = 6;
  matrix_mul = matrix1 * m;
  CompareMatrix(matrix_mul, matrix_result);
}

TEST(test_operator_add_eq, test26) {
  S21Matrix matrix1(3, 3);
  matrix1(1, 1) = 2;
  S21Matrix matrix2(3, 3);
  matrix2(1, 1) = 2;
  S21Matrix matrix_result(3, 3);
  matrix_result(1, 1) = 4;
  matrix1 += matrix2;
  CompareMatrix(matrix1, matrix_result);
}

TEST(test_operator_sub_eq, test27) {
  S21Matrix matrix1(3, 3);
  matrix1(1, 1) = 2;
  S21Matrix matrix2(3, 3);
  matrix2(1, 1) = 1;
  S21Matrix matrix_result(3, 3);
  matrix_result(1, 1) = 1;
  matrix1 -= matrix2;
  CompareMatrix(matrix1, matrix_result);
}

TEST(test_operator_sub_eq, test28) {
  S21Matrix matrix1(3, 3);
  matrix1(1, 1) = 2;
  S21Matrix matrix2(3, 3);
  matrix2(1, 1) = 1;
  S21Matrix matrix_result(3, 3);
  matrix_result(1, 1) = 1;
  matrix1 -= matrix2;
  CompareMatrix(matrix1, matrix_result);
}

TEST(test_operator_mulmatr_eq, test29) {
  S21Matrix matrix1(3, 3);
  matrix1(1, 1) = 2;
  S21Matrix matrix2(3, 3);
  matrix2(1, 1) = 2;
  S21Matrix matrix_result(3, 3);
  matrix_result(1, 1) = 4;
  matrix1 *= matrix2;
  CompareMatrix(matrix1, matrix_result);
}

TEST(test_operator_mulnum_eq, test30) {
  S21Matrix matrix1(3, 3);
  matrix1(1, 1) = 2;
  double d = 2;
  S21Matrix matrix_result(3, 3);
  matrix_result(1, 1) = 4;
  matrix1 *= d;
  CompareMatrix(matrix1, matrix_result);
}

TEST(test_operator_eq_eq, test31) {
  S21Matrix matrix1(3, 3);
  EXPECT_EQ(matrix1 == matrix1, 1);
}

TEST(test_operator_eq_eq, test32) {
  S21Matrix matrix1(3, 3);
  S21Matrix matrix2(2, 3);
  EXPECT_EQ(matrix1 == matrix2, 0);
}

TEST(test_move, test33) {
  int rows = 5, cols = 5;
  S21Matrix matrix1(rows, cols);
  S21Matrix matrix2(std::move(matrix1));
  EXPECT_EQ(matrix1.GetRows(), 0);
  EXPECT_EQ(matrix1.GetCols(), 0);
  EXPECT_EQ(matrix2.GetRows(), rows);
  EXPECT_EQ(matrix2.GetCols(), cols);
}

TEST(test_set_rows, test34) {
  S21Matrix matrix1(3, 3);
  matrix1.SetRows(4);
  EXPECT_EQ(matrix1.GetRows(), 4);
}

TEST(test_set_cols, test35) {
  S21Matrix matrix1(3, 3);
  matrix1.SetCols(4);
  EXPECT_EQ(matrix1.GetCols(), 4);
}

int main(int argc, char *argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
