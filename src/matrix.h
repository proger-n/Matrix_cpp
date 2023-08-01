#ifndef MATRIX_H_
#define MATRIX_H_

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
    private:
        // Attributes
        int rows_, cols_;         // Rows and columns
        double **matrix_;         // Pointer to the memory where the matrix is allocated

    public:
        S21Matrix();
        S21Matrix(int rows, int cols);              // Default constructor
        ~S21Matrix();             // Destructor

        bool EqMatrix(const S21Matrix& other);
        void SumMatrix(const S21Matrix& other); 
        void SubMatrix(const S21Matrix& other);
        void MulNumber(const double num);
        void MulMatrix(const S21Matrix& other);
        S21Matrix Transpose();
        S21Matrix CalcComplements();
        double Determinant();
        double &operator()(int row, int col);
        S21Matrix &operator=(const S21Matrix& A);




        void Print(){
            for(int i = 0; i < this->rows_; i++){
                for(int j = 0; j < this->cols_; j++)
                    cout << this->matrix_[i][j] << " ";
                cout << endl;    
                }
        }
        // Other methods..
        void GetMinor(double **mat, double **temp, int skip_row, int skip_col,
                   int n);
        double CalculateDeterminant(S21Matrix mat, int size);
};

#endif