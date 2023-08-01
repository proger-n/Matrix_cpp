#include <iostream>
#include <cmath>
#include "matrix.h"
using namespace std;

// Point::Point() {
//     cout << "Вызов конструктора объекта" << this << endl;
// }

// Point::~Point() {
//     cout << "Вызов деструктора объекта" << this << endl;
// }

// void Point::setCoord(int x, int y) {
//     this->x = x;
//     this->y = y;
// }

// int Point::getX() {return x;}

// int Point::getY() {return y;}
S21Matrix::S21Matrix() : S21Matrix(1, 1) {
    // cout << "Вызов конструктора объекта" << this << endl;
    // this->rows_ = 1;
    // this->cols_ = 1;
    // this->matrix_ = new double*[rows_];
    // this->matrix_[0] = new double[cols_];
}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if(this->rows_ > 0 && this->cols_ > 0) {
    cout << "Вызов конструктора объекта" << this << endl;
    this->matrix_ = new double*[rows_];
    for(int i = 0; i < this->rows_; i++)
        this->matrix_[i] = new double[cols_];
  } 
}

S21Matrix::~S21Matrix() {
    cout << "Вызов деструктора объекта" << this << endl;
    for(int i = 0; i < this->rows_; i++)
        delete[] this->matrix_[i];
    delete[] this->matrix_;
}

bool S21Matrix::EqMatrix(const S21Matrix& other) {
    bool ret = true;
  if (this->rows_ == other.rows_ && this->cols_ == other.cols_) {
      for(int i = 0; i < this->rows_; i++)
            for(int j = 0; j < this->cols_; j++)
                if(this->matrix_[i][j] * 1e7 != other.matrix_[i][j] * 1e7)
                    ret = false;
  } else ret = false;
  return ret;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
    if(this->cols_ == other.cols_) {
        for(int i = 0; i < this->rows_; i++)
            for(int j = 0; j < this->cols_; j++)
                this->matrix_[i][j] = this->matrix_[i][j] + other.matrix_[i][j];
    }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (this->rows_ == other.rows_ && this->cols_ == other.cols_) {
      for(int i = 0; i < this->rows_; i++)
            for(int j = 0; j < this->cols_; j++)
                this->matrix_[i][j] += other.matrix_[i][j];
  }
}

void S21Matrix::MulNumber(const double num) {
  if (this->rows_ > 0 && this->cols_ > 0) {
      for(int i = 0; i < this->rows_; i++)
            for(int j = 0; j < this->cols_; j++)
                this->matrix_[i][j] *= num;
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (this->cols_ == other.rows_) {
    S21Matrix result(this->rows_, other.cols_);
      for(int i = 0; i < this->rows_; i++)
            for(int j = 0; j < other.cols_; j++)
                for(int k = 0; k < this->rows_; k++)
                    result.matrix_[i][j] += this->matrix_[i][k]*other.matrix_[k][j];
    *this = result;
  }
}

S21Matrix S21Matrix::Transpose() {
  if (this->rows_ > 0 && this->cols_ > 0) {
    S21Matrix matr_result(this->rows_, this->cols_);
      for(int i = 0; i < this->rows_; i++)
            for(int j = 0; j < this->cols_; j++)
                matr_result.matrix_[i][j] = this->matrix_[i][j];
    return matr_result;
  } throw "Bad matrix!";
    
}

S21Matrix S21Matrix::CalcComplements() {
  if (this->rows_ == this->cols_) {
    S21Matrix matr_result(this->rows_, this->cols_);
    int size = (this->rows_) - 1;
    for (int i = 0; i < matr_result.rows_; i++) {
      for (int j = 0; j < matr_result.cols_; j++) {
          S21Matrix temp(size, size);
          GetMinor(this->matrix_, temp.matrix_, i, j, size + 1);
          matr_result.matrix_[i][j] =
              pow(-1.0, i + j) * CalculateDeterminant(temp, size);
        //   s21_remove_matrix(&temp);
        // temp.~S21Matrix();
      }
    }
    return matr_result;
  } throw "Bad matrix!";
    
}

double S21Matrix::Determinant() {
  if (this->rows_ == this->cols_) {
    return CalculateDeterminant(*this, this->rows_);
  } throw "Bad matrix!";
  
}

double &S21Matrix::operator()(int row, int col){
    return this->matrix_[row][col];
}

S21Matrix &S21Matrix::operator=(const S21Matrix& other) {
    if (this->rows_ == other.rows_ && this->cols_ == other.cols_) {
      for(int i = 0; i < this->rows_; i++)
        for(int j = 0; j < this->cols_; j++)
          this->matrix_[i][j] = other.matrix_[i][j];
      return *this;
    }
    throw "matrixes not equal!";
}

void S21Matrix::GetMinor(double **mat, double **temp, int skip_row, int skip_col,
                   int n) {
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

double S21Matrix::CalculateDeterminant(S21Matrix mat, int size) {
    double D = 0;
  if (size == 1)
    D = mat.matrix_[0][0];
  else {
    S21Matrix temp(size, size);

    int sign = 1;

    for (int cur_col = 0; cur_col < size; cur_col++) {
      GetMinor(mat.matrix_, temp.matrix_, 0, cur_col, size);
      D += sign * mat.matrix_[0][cur_col] *
           CalculateDeterminant(temp, size - 1);

      sign = -sign;
    }
    // s21_remove_matrix(&temp);
    temp.~S21Matrix();
  }
  return D;
}

 int main() {
    setlocale(LC_ALL, "rus");

    S21Matrix matr1(2,2);
    S21Matrix matr2(2,2);
    S21Matrix matr3;
    matr1(0, 0) = 1.1;
    matr1(0, 1) = 1;
    matr1(1, 0) = 1;
    matr1(1, 1) = 1;

    matr2(0, 0) = 1;
    matr2(0, 1) = 1;
    matr2(1, 0) = 1;
    matr2(1, 1) = 1;
    matr3.SumMatrix(matr2);
    matr1.Print();
    // matr2.CalcComplements();
    S21Matrix matr4(2,2);
    matr4=matr1;
    matr4.Print();
    // cout << "eq=" << matr2.Determinant() << endl;
    
    //  Point pt;
    //  pt.setCoord(2, 3);

    //  Point *ptr = new Point();
    //  ptr->setCoord(4, 5);

    //  cout << pt.getX() << " " << pt.getY() << endl;
    //  cout << ptr->getX() << " " << ptr->getY() << endl;
     return 0;
 }