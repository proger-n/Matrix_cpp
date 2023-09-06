#ifndef MATRIX_H_
#define MATRIX_H_
#include <iostream>

class S21Matrix {
 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other);
  ~S21Matrix();


  bool EqMatrix(const S21Matrix& other);
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  double& operator()(int row, int col);
  double operator()(int row, int col) const;
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator=(S21Matrix&& A) noexcept;
  S21Matrix& operator+(const S21Matrix& other);
  S21Matrix& operator-(const S21Matrix& other);
  S21Matrix& operator*(const S21Matrix& other);
  S21Matrix& operator*(double d);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator*=(const double d);
  bool operator==(const S21Matrix& other);

  // temp funcs, delete them
  void Print() {
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++)
        std::cout << this->matrix_[i][j] << " ";
      std::cout << std::endl;
    }
  }

 private:
  // Attributes
  int rows_, cols_;
  double** matrix_;
  // Other methods..
  void GetMinor(double** mat, double** temp, int skip_row, int skip_col, int n);
  double CalculateDeterminant(double** mat, int size);
};

#endif