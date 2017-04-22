# Matrix Public Functions

## Matrix Constructors and Destructor  

~Matrix() = default;  

Matrix() = default;  

Matrix(const std::initializer_list\<std::initializer_list\<double\>\>& matrix);  

Matrix(const Matrix& matrix);  

Matrix(Matrix&& matrix);  

explicit Matrix(const Vector\<Vector\<double\>\>& matrix);  

explicit Matrix(size_t rows, size_t columns);  

## Matrix Public Static Functions  

static Matrix rVector(Vector\<double\> vector);  

static Matrix cVector(Vector\<double\> vector);  

static Matrix diagonal(Vector\<double\> vector);  

static Matrix identity(size_t power);  

## Matrix Print Funtion  

void print();  

## Matrix Elementary Transform Functions  

Matrix& rSwap(size_t ri, size_t rj);  

Matrix& rMultiply(size_t ri, double n);  

Matrix& rTransform(size_t ri, size_t rj, double n);  

Matrix& cSwap(size_t ci, size_t cj);  

Matrix& cMultiply(size_t ci, double n);  

Matrix& cTransform(size_t ci, size_t cj, double n);  

## Determinant of Square Matrix  

double determinant();  

## Inverse Matrix of Square Matrix  

Matrix inverse();  

## Adjoint Matrix of Square Matrix   

Matrix adjoint();  

## Remove(、Insert or Replace) Row(or Column) of Matrix  

Matrix& rRemove(size_t row);  

Matrix& rInsert(size_t row,const Matrix& rVector);  

Matrix& rInsert(size_t row,Vector\<double\> rVector);  

Matrix& rReplace(size_t row,const Matrix& rVector);  

Matrix& rReplace(size_t row,Vector\<double\> rVector);   

Matrix& cRemove(size_t column);  

Matrix& cInsert(size_t column,const Matrix& cVector);  

Matrix& cInsert(size_t column,Vector\<double\> cVector);  

Matrix& cReplace(size_t column,const Matrix& cVector);  

Matrix& cReplace(size_t column,Vector\<double\> cVector);  

## Get a Column From Matrix  

Matrix column(size_t column);  

## Trust Current as Extended Matrix to Calulate Solve Matrix by Cramer's Rule  

Matrix cramer();  

## To Set Matrix's Name  

Matrix& setName(std::string name);  

## Get Matrix's Name

std::string name();  

## Get Number of Rows  

size_t rows();  

## Get Number of Columns  

size_t columns();  

## Check the Matrix if it is Available  

bool isAvailable();  

## Check the Operation if it is Successful   

bool isSuccessful();  

## Override of Operation

Matrix& operator = (const Matrix& matrix);  

Matrix& operator = (Matrix&& matrix);  

Vector\<double\>& operator [] (size_t i);  

const Vector\<double\>& operator [] (size_t i) const;  

friend bool operator == (const Matrix& x, const Matrix& y);  

friend bool operator != (const Matrix& x, const Matrix& y);  

friend Matrix operator * (const Matrix& x, double n);  

friend Matrix operator * (double n, const Matrix& x);  

friend Matrix operator + (const Matrix& x, const Matrix& y);  

friend Matrix operator - (const Matrix& x, const Matrix& y);  

friend Matrix operator * (const Matrix& x, const Matrix& y);  
