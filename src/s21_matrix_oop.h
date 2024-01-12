#ifndef S21_MATRIX_H_ON_CPP_MY_MATRIX_OOP_H_
#define S21_MATRIX_H_ON_CPP_MY_MATRIX_OOP_H_
#include <iostream>

class S21Matrix {
 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other) noexcept;
  ~S21Matrix();

  double& operator()(int row, int col);
  double operator()(int row, int col) const;
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix operator*(double num);
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator=(S21Matrix&& other) noexcept;
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator*=(const double num);
  bool operator==(const S21Matrix& other);

  bool EqMatrix(const S21Matrix& other) noexcept;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num) noexcept;
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  int GetRows() const;
  int GetCols() const;
  void SetRows(int rows);
  void SetCols(int cols);

 private:
  // Attributes
  int rows_, cols_;
  double** matrix_;
  // Other methods..
  void GetMinor(double** mat, double** temp, int skip_row, int skip_col, int n);
  double CalculateDeterminant(double** mat, int size);
  void CopyMatrix(double** matrix, int rows, int cols);
};

#endif  // S21_MATRIX_H_ON_CPP_MY_MATRIX_OOP_H_