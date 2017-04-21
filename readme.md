# Matrix 
~Matrix() = default;  

Matrix() = default;  

Matrix(const std::initializer_list\<std::initializer_list\<double\>\>& matrix);  

Matrix(const Matrix& matrix);  

Matrix(Matrix&& matrix);  

explicit Matrix(const Vector\<Vector\<double\>\>& matrix);  

explicit Matrix(size_t rows, size_t columns);  

static Matrix rVector(Vector\<double\> vector);  

static Matrix cVector(Vector\<double\> vector);  

static Matrix diagonal(Vector\<double\> vector);  

static Matrix identity(size_t power);  

void print();  

Matrix& rSwap(size_t ri, size_t rj);  

Matrix& rMultiply(size_t ri, double n);  

Matrix& rTransform(size_t ri, size_t rj, double n);  

Matrix& cSwap(size_t ci, size_t cj);  

Matrix& cMultiply(size_t ci, double n);  

Matrix& cTransform(size_t ci, size_t cj, double n);  

double determinant();  

Matrix inverse();  

Matrix adjoint();  

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

Matrix column(size_t column);  

Matrix cramer();  

Matrix& setName(std::string name)  

std::string name()  

size_t rows()  

size_t columns()  

bool isAvailable()  

bool isSuccessful()  

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
