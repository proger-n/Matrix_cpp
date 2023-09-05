#include "matrix.h"

#include <cmath>
#include <iostream>
using namespace std;

S21Matrix::S21Matrix() : S21Matrix(1, 1) {
  // cout << "Вызов конструктора объекта" << this << endl;
  // this->rows_ = 1;
  // this->cols_ = 1;
  // this->matrix_ = new double*[rows_];
  // this->matrix_[0] = new double[cols_];
}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows_ > 0 && cols_ > 0) {
    // cout << "Вызов конструктора объекта" << this << endl;
    matrix_ = new double*[rows_]();
    for (int i = 0; i < rows_; i++) {
      matrix_[i] = new double[cols_]();
    }
  }
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : S21Matrix(other.rows_, other.cols_) {
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
}

// S21Matrix::S21Matrix(S21Matrix&& other) {
//   this->~S21Matrix();
//   rows_=other.rows_;
//   cols_=other.cols_;
//   matrix_=other.matrix_;
//   other.matrix_ = nullptr;
//   other.rows_ = 0;
//   other.cols_ = 0;
//   return *this;
// }

S21Matrix::~S21Matrix() {
  // cout << "Вызов деструктора объекта" << this << endl;
  for (int i = 0; i < this->rows_; i++) delete[] this->matrix_[i];
  delete[] this->matrix_;
}

double S21Matrix::getElement(int i, int j) { return matrix_[i][j]; }

void S21Matrix::setElement(int i, int j, double value) {
  matrix_[i][j] = value;
}

bool S21Matrix::EqMatrix(const S21Matrix& other) {
  bool ret = true;
  if (this->rows_ == other.rows_ && this->cols_ == other.cols_) {
    for (int i = 0; i < this->rows_; i++)
      for (int j = 0; j < this->cols_; j++)
        if (this->matrix_[i][j] * 1e7 != other.matrix_[i][j] * 1e7) ret = false;
  } else
    ret = false;
  return ret;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (this->cols_ == other.cols_) {
    for (int i = 0; i < this->rows_; i++)
      for (int j = 0; j < this->cols_; j++)
        this->matrix_[i][j] += other.matrix_[i][j];
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (this->rows_ == other.rows_ && this->cols_ == other.cols_) {
    for (int i = 0; i < this->rows_; i++)
      for (int j = 0; j < this->cols_; j++)
        this->matrix_[i][j] -= other.matrix_[i][j];
  }
}

void S21Matrix::MulNumber(const double num) {
  if (this->rows_ > 0 && this->cols_ > 0) {
    for (int i = 0; i < this->rows_; i++)
      for (int j = 0; j < this->cols_; j++) this->matrix_[i][j] *= num;
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (this->cols_ == other.rows_) {
    S21Matrix result(this->rows_, other.cols_);
    for (int i = 0; i < this->rows_; i++)
      for (int j = 0; j < other.cols_; j++)
        for (int k = 0; k < this->rows_; k++)
          result.matrix_[i][j] += this->matrix_[i][k] * other.matrix_[k][j];
    *this = result;
  }
}

S21Matrix S21Matrix::Transpose() {
  if (rows_ > 0 && cols_ > 0) {
    S21Matrix matr_result(cols_, rows_);
    for (int i = 0; i < matr_result.rows_; i++)
      for (int j = 0; j < matr_result.cols_; j++)
        matr_result.matrix_[i][j] = matrix_[j][i];
    return matr_result;
  }
  throw "Bad matrix!";
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ == cols_) {
    S21Matrix matr_result(rows_, cols_);
    int size = (rows_)-1;
    for (int i = 0; i < matr_result.rows_; i++) {
      for (int j = 0; j < matr_result.cols_; j++) {
        S21Matrix temp(size, size);
        GetMinor(matrix_, temp.matrix_, i, j, size + 1);
        matr_result.matrix_[i][j] =
            pow(-1.0, i + j) * CalculateDeterminant(temp.matrix_, size);
        //   s21_remove_matrix(&temp);
        // temp.~S21Matrix();
      }
    }
    return matr_result;
  }
  throw "Bad matrix!";
}

double S21Matrix::Determinant() {
  // if (this->rows_ == this->cols_) {
  return CalculateDeterminant(matrix_, rows_);
  // } throw "Bad matrix!";
}

S21Matrix S21Matrix::InverseMatrix() {
  if (rows_ == cols_) {
    double det = 0.0;
    det = this->Determinant();
    if (det != 0) {
      S21Matrix matr_result(rows_, cols_);
      S21Matrix trans(rows_, cols_);
      S21Matrix comp(rows_, cols_);
      trans = this->Transpose();
      comp = trans.CalcComplements();
      for (int i = 0; i < matr_result.rows_; i++) {
        for (int j = 0; j < matr_result.cols_; j++) {
          matr_result.matrix_[i][j] = comp.matrix_[i][j] / det;
          //   s21_remove_matrix(&temp);
          // temp.~S21Matrix();
        }
      }
      return matr_result;
    }
  }
  throw "Bad matrix!";
}

double& S21Matrix::operator()(int row, int col) {
  if (this->rows_ > 0 && this->cols_ > 0) {
    return this->matrix_[row][col];
  }
  throw "Out of range!";
}

double S21Matrix::operator()(int row, int col) const {
  return this->matrix_[row][col];
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  for (int i = 0; i < this->rows_; i++)
    for (int j = 0; j < this->cols_; j++)
      this->matrix_[i][j] = other.matrix_[i][j];
  return *this;
}

// S21Matrix &S21Matrix::operator=(S21Matrix&& other) {
//   this->~S21Matrix();
//   rows_=other.rows_;
//   cols_=other.cols_;
//   matrix_=other.matrix_;
//   other.matrix_ = nullptr;
//   other.rows_ = 0;
//   other.cols_ = 0;
//   return *this;
// }

S21Matrix& S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

S21Matrix& S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix result(*this);
  result.SubMatrix(other);
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
    // s21_remove_matrix(&temp);
    // temp.~S21Matrix();
  }
  return D;
}

int main() {
  setlocale(LC_ALL, "rus");

  S21Matrix matr1(3, 3);
  S21Matrix matr2(2, 2);
  S21Matrix matr3;
  matr1(0, 0) = 1.1;
  matr1(0, 1) = 1;
  matr1(0, 2) = 1;
  matr1(1, 0) = 1;
  matr1(1, 1) = 1;
  matr1(1, 2) = 1;

  matr2(0, 0) = 1;
  matr2(0, 1) = 2;
  matr2(1, 0) = 3;
  matr2(1, 1) = 4;

  matr1.Print();
  matr1.MulNumber(3);
  matr1.Print();

  S21Matrix matr4(2, 2);
  matr4 = matr2.InverseMatrix();  // оператор копирования
  // matr4=std::move(matr1); //     оператор перемещения
  matr2.Print();
  matr4.Print();

  cout << "eq 1 and 2=" << matr1.EqMatrix(matr2) << endl;
  cout << "eq 1 and 4=" << matr1.EqMatrix(matr4) << endl;
  cout << "det=" << matr2.Determinant() << endl;
  matr2-=matr4;
  matr2.Print();
  return 0;
}