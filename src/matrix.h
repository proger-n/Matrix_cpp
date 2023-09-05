#ifndef MATRIX_H_
#define MATRIX_H_
#include <iostream>
using namespace std;
// class Point {
// public:
//     Point();
//     ~Point();
// private:
//     int x, y;

// public:
//     void setCoord(int pt_x, int pt_y);
//     int getX();
//     int getY();
// };

class S21Matrix {
 public:
  S21Matrix();
  S21Matrix(int rows, int cols);  // Default constructor
  S21Matrix(const S21Matrix& other);
  // S21Matrix(S21Matrix&& other);
  ~S21Matrix();  // Destructor

  double getElement(int i, int j);
  void setElement(int i, int j, double value);

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
  // S21Matrix &operator=(S21Matrix&& A) noexcept;
  S21Matrix& operator+(const S21Matrix& other);
  S21Matrix& operator-(const S21Matrix& other);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);

  // temp funcs, delete them
  void Print() {
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++) cout << this->matrix_[i][j] << " ";
      cout << endl;
    }
  }

 private:
  // Attributes
  int rows_, cols_;  // Rows and columns
  double** matrix_;  // Pointer to the memory where the matrix is allocated
  // Other methods..
  void GetMinor(double** mat, double** temp, int skip_row, int skip_col, int n);
  double CalculateDeterminant(double** mat, int size);
};

#endif