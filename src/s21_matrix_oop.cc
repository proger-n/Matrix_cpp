#include "s21_matrix_oop.h"

#include <cmath>
#include <iostream>

S21Matrix::S21Matrix() : S21Matrix(1, 1) {}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows_ > 0 && cols_ > 0) {
    matrix_ = new double*[rows_]();
    for (int i = 0; i < rows_; i++) {
      matrix_[i] = new double[cols_]();
    }
  }
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : S21Matrix(other.rows_, other.cols_) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix&& other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
}

S21Matrix::~S21Matrix() {
  for (int i = 0; i < this->rows_; i++) {
    delete[] matrix_[i];
  }
  delete[] matrix_;
}

bool S21Matrix::EqMatrix(const S21Matrix& other) {
  bool ret = true;
  if (rows_ == other.rows_ && cols_ == other.cols_) {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        if (matrix_[i][j] * 1e7 != other.matrix_[i][j] * 1e7) ret = false;
      }
    }
  } else
    ret = false;
  return ret;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (cols_ != other.cols_) {
    throw std::out_of_range("ERROR: check matrix size!");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (cols_ != other.cols_) {
    throw std::out_of_range("ERROR: check matrix size!");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::out_of_range("ERROR: check matrix size!");
  }
  S21Matrix result(rows_, other.cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      for (int k = 0; k < rows_; k++) {
        result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  *this = result;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix matr_result(cols_, rows_);
  for (int i = 0; i < matr_result.rows_; i++) {
    for (int j = 0; j < matr_result.cols_; j++) {
      matr_result.matrix_[i][j] = matrix_[j][i];
    }
  }
  return matr_result;
}

S21Matrix S21Matrix::CalcComplements() {
  S21Matrix matr_result(rows_, cols_);

  if (rows_ != cols_) {
    throw std::out_of_range("ERROR: check matrix size!");
  }
  int size = (rows_)-1;
  for (int i = 0; i < matr_result.rows_; i++) {
    for (int j = 0; j < matr_result.cols_; j++) {
      S21Matrix temp(size, size);
      GetMinor(matrix_, temp.matrix_, i, j, size + 1);
      matr_result.matrix_[i][j] =
          pow(-1.0, i + j) * CalculateDeterminant(temp.matrix_, size);
    }
  }

  return matr_result;
}

double S21Matrix::Determinant() {
  double result = 0.0;
  if (rows_ != cols_) {
    throw std::out_of_range("ERROR: check matrix size!");
  }
  result = CalculateDeterminant(matrix_, rows_);
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  S21Matrix matr_result(rows_, cols_);
  if (rows_ != cols_) {
    throw std::out_of_range("ERROR: check matrix size!");
  }
  double det = 0.0;
  det = this->Determinant();
  if (det == 0) {
    throw std::out_of_range("ERROR: det=0!");
  }
  S21Matrix trans(rows_, cols_);
  S21Matrix comp(rows_, cols_);
  trans = this->Transpose();
  comp = trans.CalcComplements();
  for (int i = 0; i < matr_result.rows_; i++) {
    for (int j = 0; j < matr_result.cols_; j++) {
      matr_result.matrix_[i][j] = comp.matrix_[i][j] / det;
    }
  }
  return matr_result;
}

double& S21Matrix::operator()(int row, int col) {
  if (row > rows_ || col > cols_ || row < 0 || col < 0) {
    throw std::out_of_range("ERROR: check matrix size!");
  }
  return matrix_[row][col];
}

double S21Matrix::operator()(int row, int col) const {
  if (row > rows_ || col > cols_ || row < 0 || col < 0) {
    throw std::out_of_range("ERROR: check matrix size!");
  }
  return matrix_[row][col];
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if(this != &other) {
    this->~S21Matrix();
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = new double*[rows_]();
      for (int i = 0; i < rows_; i++) {
        matrix_[i] = new double[cols_]();
      }
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] = other.matrix_[i][j];
      }
    }
  }
  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& other) noexcept {
  this->~S21Matrix();
  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_ = other.matrix_;
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix result(*this);
  result.SubMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix result(*this);
  result.MulMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(double d) {
  S21Matrix result(*this);
  result.MulNumber(d);
  return result;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  this->SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  this->SubMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  this->MulMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(double d) {
  this->MulNumber(d);
  return *this;
}

bool S21Matrix::operator==(const S21Matrix& other) {
  return EqMatrix(other);
}

void S21Matrix::GetMinor(double** mat, double** temp, int skip_row,
                         int skip_col, int n) {
  int i = 0, j = 0;
  for (int row = 0; row < n; row++) {
    for (int col = 0; col < n; col++) {
      if (row != skip_row && col != skip_col) {
        temp[i][j++] = mat[row][col];
        if (j == n - 1) {
          j = 0;
          i++;
        }
      }
    }
  }
}

double S21Matrix::CalculateDeterminant(double** mat, int size) {
  double D = 0;
  if (size == 1)
    D = mat[0][0];
  else {
    S21Matrix temp(size, size);
    int sign = 1;
    for (int cur_col = 0; cur_col < size; cur_col++) {
      GetMinor(mat, temp.matrix_, 0, cur_col, size);
      D +=
          sign * mat[0][cur_col] * CalculateDeterminant(temp.matrix_, size - 1);
      sign = -sign;
    }
  }
  return D;
}

int S21Matrix::get_row_size() const { return rows_; }

int S21Matrix::get_col_size() const { return cols_; }

// void S21Matrix::set_rows(int rows) {
//   if (rows < 1) throw std::out_of_range("Incorrect rows size");
//   if (rows != rows_) {
//     S21Matrix result(rows, cols_);
//     PassValues(result);
//     *this = result;
//   }
// }

// void S21Matrix::set_cols(int cols) {
//   if (cols < 1) throw std::out_of_range("Incorrect cols size");
//   if (cols != cols_) {
//     S21Matrix result(rows_, cols);
//     PassValues(result);
//     *this = result;
//   }
// }

// int main() {
//   setlocale(LC_ALL, "rus");

//   S21Matrix matr1(3, 3);
//   S21Matrix matr2(2, 2);
//   S21Matrix matr3(2, 2);
//   matr1(0, 0) = 1.1;
//   matr1(0, 1) = 1;
//   matr1(0, 2) = 1;
//   matr1(1, 0) = 1;
//   matr1(1, 1) = 1;
//   matr1(1, 2) = 1;

//   matr2(0, 0) = 1;
//   matr2(0, 1) = 2;
//   matr2(1, 0) = 3;
//   matr2(1, 1) = 4;

//   // matr1.Print();
//   // matr1.MulNumber(3);
//   // matr1.Print();

//   // matr2.SumMatrix(matr1);
//   matr2.Print();
//   S21Matrix matr4(2, 2);
//   // matr4 = matr2.InverseMatrix();  // оператор копирования
//   matr4 = std::move(matr1);  //     оператор перемещения
//   // matr2.Print();
//   matr4.Print();
//   matr1.Print();

//   // cout << "eq 1 and 2=" << matr1.EqMatrix(matr2) << endl;
//   // cout << "eq 1 and 4=" << matr1.EqMatrix(matr4) << endl;
//   // cout << "det=" << matr2.Determinant() << endl;
//   // matr3=matr1+matr4;
//   // matr2 *= 3.3;
//   // matr2.Print();
//   return 0;
// }