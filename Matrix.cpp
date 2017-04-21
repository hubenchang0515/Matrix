/*
	Author : PlanC
	E-mail : hubenchang0515@outlook.com
	Blog   : www.kurukurumi.com

	File   : Matrix.c
	Date   : 2017-4-22
*/

#include <Matrix.h>
#include <iostream>
#include <iomanip>
#include <stdexcept>

#define ERROR_MESSAGE(MESSAGE) errorString( __FUNCTION__ ,MESSAGE)

Matrix::Matrix(size_t rows_, size_t columns_)
{
	this->matrix = Vector<Vector<double>>(rows_);
	this->rows_ = rows_;
	this->columns_ = columns_;

	for(auto& row : this->matrix)
		row.resize(columns_,0);
}

Matrix::Matrix(const std::initializer_list<std::initializer_list<double>>& matrix)
{
	this->rows_ = matrix.size();
	this->columns_ = matrix.begin()->size();
	for(auto row : matrix)
		this->matrix.push_back(row);
}

Matrix::Matrix(const Matrix& matrix)
{
	this->rows_ = matrix.rows_;
	this->columns_ = matrix.columns_;
	this->matrix = Vector<Vector<double>>(matrix.rows_);
	for(size_t row = 0; row < matrix.rows_; row++)
		this->matrix[row] = matrix.matrix[row];
}

Matrix::Matrix(Matrix&& matrix)
{
	this->rows_ = matrix.rows_;
	this->columns_ = matrix.columns_;
	this->matrix = Vector<Vector<double>>(matrix.rows_);
	for(size_t row = 0; row < matrix.rows_; row++)
		this->matrix[row] = matrix.matrix[row];
}

Matrix::Matrix(const Vector<Vector<double>>& matrix)
{
	this->matrix = matrix;
}


Matrix Matrix::rVector(Vector<double> vector)
{
	Matrix ans(1,vector.size());
	ans[0] = vector;
	return ans;
}

Matrix Matrix::cVector(Vector<double> vector)
{
	size_t rows = vector.size();
	Matrix ans(rows,1);
	for(size_t row = 0; row < rows; row++)
	{
		ans[row][0] = vector[row];
	}

	return ans;
}

Matrix Matrix::diagonal(Vector<double> vector)
{
	size_t rows = vector.size();
	Matrix ans(rows,rows);
	for(size_t row = 0; row < rows; row++)
	{
		ans[row][row] = vector[row];
	}

	return ans;
}

Matrix Matrix::identity(size_t power)
{
	Matrix ans(power,power);
	for(size_t element = 0; element < power; element++)
	{
		ans[element][element] = 1;
	}

	return ans;
}

void Matrix::print()
{
	std::cout << std::endl << name_ << " = " << std::endl;
	if(available_)
	{
		for(auto& row : this->matrix)
		{
			for(auto element : row)
			{
				std::cout << std::setw(10) << element << " ";
			}
			std::cout << ";" << std::endl;
		}
	}
	else
	{
		std::cout << "This Matrix is not available.\n";
	}
}


Matrix& Matrix::rSwap(size_t ri, size_t rj)
{
	if(ri >= rows_ || rj >= rows_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Row number out of range."));
	}
	else
	{
		std::swap(this->matrix[ri], this->matrix[rj]);
		successful_ = true;
	}

	return *this;
}

Matrix& Matrix::rMultiply(size_t ri, double n)
{
	if (ri >= rows_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Row number out of range."));
	}
	else
	{
		for(auto& element : this->matrix[ri])
		{
			element *= n;
		}
		successful_ = true;
	}

	return *this;
}

Matrix& Matrix::rTransform(size_t ri, size_t rj, double n)
{
	if(ri >= rows_ || rj >= rows_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Row number out of range."));
	}
	else
	{
		for(size_t element = 0; element < columns_; element++)
		{
			this->matrix[ri][element] += n * this->matrix[rj][element];
		}
		successful_ = true;
	}

	return *this;
}


Matrix& Matrix::cSwap(size_t ci, size_t cj)
{
	if(ci >= columns_ || cj >= columns_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Column number out of range."));
	}
	else
	{
		for(size_t row = 0; row < rows_; row++)
		{
			std::swap(this->matrix[row][ci],this->matrix[row][cj]);
		}
		successful_ = true;
	}

	return *this;
}

Matrix& Matrix::cMultiply(size_t ci, double n)
{
	if(ci >= columns_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Column number out of range."));
	}
	else
	{
		for(auto& row : this->matrix)
		{
			row[ci] *= n;
		}
		successful_ = true;
	}

	return *this;
}

Matrix& Matrix::cTransform(size_t ci, size_t cj, double n)
{
	if(ci >= columns_ || cj >= columns_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Column number out of range."));
	}
	else
	{
		for(size_t row = 0; row < rows_; row++)
		{
			this->matrix[row][ci] += n * this->matrix[row][cj];
		}
		successful_ = true;
	}

	return *this;
}

double Matrix::determinant()
{
	double result = 1;
	if(rows_ != columns_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Not a square matrix."));
	}
	else
	{
		successful_ =  true;
		Matrix temp = *this;
		/* Transform temp Matrix to be a lower triangular matrix. */
		for(size_t column = 0; column < temp.columns_-1; column++)
		{
			/* Ensure temp[column][column] isn't 0. */
			for(size_t i = column+1; temp[column][column] == 0 && i < temp.columns_; i++)
			{
				temp.cSwap(column,i);
				result *= -1;
			}
			if(temp[column][column] == 0)
			{
				return 0;
			}

			/* change current column lower elements to be 0 */
			for(size_t row = column+1; row < temp.rows_; row++)
			{
				temp.rTransform(row,column,-temp[row][column]/temp[column][column]);
			}
		}

		for(size_t i = 0; i < temp.rows_; i++)
		{
			result *= temp[i][i];
		}
		
	}

	return result;
}

Matrix Matrix::inverse()
{
	Matrix ans;
	if(rows_ != columns_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Not a square matrix."));
	}
	else
	{
		Matrix temp = *this;
		ans = Matrix::identity(this->rows_);
		/* Transform this Matrix to be an identity matrix. */
		for(size_t column = 0; column < temp.columns_; column++)
		{
			/* Ensure temp[column][column] isn't 0. */
			for(size_t row = column + 1; temp[column][column] == 0 && row < temp.rows_; row++)
			{
				temp.rSwap(column,row);
				ans.rSwap(column,row);
			}
			if(temp[column][column] == 0)
			{
				successful_ = false;
				throw std::logic_error(ERROR_MESSAGE("Don't have inverse matrix."));
				return ans;
			}

			/* change diagonal element to be 1 */
			double k = 1/temp[column][column];
			temp.rMultiply(column,k);
			ans.rMultiply(column,k);

			/* change current column upper elements to be 0 */
			for(size_t row = 0; row < column; row++)
			{
				double k = -temp[row][column]/temp[column][column];
				temp.rTransform(row,column,k);
				ans.rTransform(row,column,k);
			}
			/* change current column lower elements to be 0 */
			for(size_t row = column + 1; row < temp.rows_; row++)
			{
				double k = -temp[row][column]/temp[column][column];
				temp.rTransform(row,column,k);
				ans.rTransform(row,column,k);
			}
			
		}
	}

	return ans;
}

Matrix Matrix::adjoint()
{
	if(rows_ != columns_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Not a square matrix."));
		return Matrix();
	}
	return this->determinant() * this->inverse();
}

Matrix& Matrix::rRemove(size_t row)
{
	if(row >= this->rows_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Row number out of range."));
	}
	else
	{
		successful_ = true;
		Vector<Vector<double>>::iterator r = this->matrix.begin() + row;
		this->matrix.erase(r);
		this->rows_ -= 1;
	}

	return *this;
}
	
Matrix& Matrix::rInsert(size_t row,const Matrix& rVector)
{
	if(row > this->rows_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Row number out of range."));
	}
	else if(rVector.rows_ != 1)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Parameter isn't a row vector matrix."));
	}
	else if(rVector.columns_ != this->columns_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Parameter length is't correct."));
	}
	else
	{
		successful_ = true;
		Vector<Vector<double>>::iterator r = this->matrix.begin() + row;
		this->matrix.insert(r,rVector[0]);
		this->rows_ += 1;
	}

	return *this;
}
	
Matrix& Matrix::rInsert(size_t row,Vector<double> rVector)
{
	if(row > this->rows_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Row number out of range."));
	}
	else if(rVector.size() != this->columns_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Parameter length is't correct."));
	}
	else
	{
		successful_ = true;
		Vector<Vector<double>>::iterator r = this->matrix.begin() + row;
		this->matrix.insert(r,rVector);
		this->rows_ += 1;
	}

	return *this;
}

Matrix& Matrix::rReplace(size_t row,const Matrix& rVector)
{
	if(row >= this->rows_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Row number out of range."));
	}
	else if(rVector.rows_ != 1)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Parameter isn't a row vector matrix."));
	}
	else if(rVector.columns_ != this->columns_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Parameter length is't correct."));
	}
	else
	{
		rRemove(row);
		rInsert(row,rVector);
	}

	return *this;
}

Matrix& Matrix::rReplace(size_t row,Vector<double> rVector)
{
	if(row >= this->rows_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Row number out of range."));
		return *this;
	}
	else if(rVector.size() != this->columns_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("parameter length is't correct."));
	}
	else
	{
		rRemove(row);
		rInsert(row,rVector);	
	}
	
	return *this;
}

Matrix& Matrix::cRemove(size_t column)
{
	if(column >= this->columns_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Column number out of range."));
	}
	else
	{
		successful_ = true;
		for(auto& row : this->matrix)
		{
			Vector<double>::iterator c = row.begin() + column;
			row.erase(c);
		}

		this->columns_ -= 1;
	}

	return *this;
}

Matrix& Matrix::cInsert(size_t column,const Matrix& cVector)
{
	if(column > this->columns_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Column number out of range."));
	}
	else if(cVector.columns_ != 1)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Parameter isn't a column vector matrix."));
	}
	else if(cVector.rows_ != this->rows_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Parameter length is't correct."));
	}
	else
	{
		successful_ = true;
		for(size_t row = 0; row < this->rows_; row++)
		{
			Vector<double>::iterator c = (this->matrix)[row].begin() + column;
			(this->matrix)[row].insert(c,cVector[row][0]);
		}
		this->columns_ += 1;
	}

	return *this;
}

Matrix& Matrix::cInsert(size_t column,Vector<double> cVector)
{
	if(column > this->columns_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Column number out of range."));
	}
	else if(cVector.size() != this->rows_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Parameter length is't correct."));
	}
	else
	{
		successful_ = true;
		for(size_t row = 0; row < this->rows_; row++)
		{
			Vector<double>::iterator c = (this->matrix)[row].begin() + column;
			(this->matrix)[row].insert(c,cVector[row]);
		}
		this->columns_ += 1;
	}

	return *this;
}

Matrix& Matrix::cReplace(size_t column,const Matrix& cVector)
{
	if(column >= this->columns_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Column number out of range."));
	}
	else if(cVector.columns_ != 1)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Parameter isn't a column vector matrix."));
	}
	else if(cVector.rows_ != this->rows_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Parameter length is't correct."));
	}
	else
	{
		cRemove(column);
		cInsert(column,cVector);
	}

	return *this;
}

Matrix& Matrix::cReplace(size_t column,Vector<double> cVector)
{
	if(column >= this->columns_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Column number out of range."));
	}
	else if(cVector.size() != this->rows_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Parameter length is't correct."));
	}
	else
	{
		cRemove(column);
		cInsert(column,cVector);
	}

	return *this;
}

Matrix Matrix::column(size_t column)
{
	Matrix ans(this->rows_,1);
	if(column >= this->columns_)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Column number out of range."));
	}
	else
	{
		for(size_t i = 0; i < rows_; i++)
		{
			ans[i][0] = this->matrix[i][column];
		}
	}

	return ans;
}

Matrix Matrix::cramer()
{
	Matrix result(rows_,1);
	Matrix A = *this;
	Matrix B = this->column(this->columns_-1);
	A.cRemove(A.columns_-1);
	double denominator = A.determinant();
	if(denominator == 0)
	{
		successful_ = false;
		throw std::logic_error(ERROR_MESSAGE("Don't have singular solve."));
		return result;
	}
	else
	{
		for(size_t i = 0; i < A.columns_; i++)
		{
			A.cReplace(i,B);
			result[i][0] = A.determinant()/denominator;
			A.cReplace(i,this->column(i));
		}
	}

	return result;
}

Vector<double>& Matrix::operator [] (size_t i)
{
	return this->matrix[i];
}

const Vector<double>& Matrix::operator [] (size_t i) const
{
	return this->matrix[i];
}

Matrix& Matrix::operator = (const Matrix& matrix)
{
	
	this->available_ = matrix.available_;
	this->rows_ = matrix.rows_;
	this->columns_ = matrix.columns_;
	this->matrix = matrix.matrix;
	return *this;
}

Matrix& Matrix::operator = (Matrix&& matrix)
{
	this->available_ = matrix.available_;
	this->rows_ = matrix.rows_;
	this->columns_ = matrix.columns_;
	this->matrix = matrix.matrix;
	return *this;
}


bool operator == (const Matrix& x, const Matrix& y)
{
	return x.matrix == y.matrix;
}

bool operator != (const Matrix& x, const Matrix& y)
{
	return !(x.matrix == y.matrix);
}

Matrix operator * (const Matrix& x, double n)
{
	return n * x;
}


Matrix operator * (double n,const Matrix& x)
{
	Matrix ans = x;
	for(auto& row : ans.matrix)
		for(auto& element : row)
		{
			element *= n;
		}

	return ans;
}

Matrix operator + (const Matrix& x,const Matrix& y)
{
	Matrix ans;
	if(x.rows_ != y.rows_ || x.columns_ != y.columns_)
	{
		ans.available_ = false;
	}
	else
	{
		ans = x;
		for(size_t row = 0;  row < ans.rows_; row++)
			for(size_t column = 0; column < ans.columns_; column++)
			{
				ans.matrix[row][column] = x.matrix[row][column] + y.matrix[row][column];
			}
	}

	return ans;
}

Matrix operator - (const Matrix& x,const Matrix& y)
{
	Matrix ans;
	if(x.rows_ != y.rows_ || x.columns_ != y.columns_)
	{
		ans.available_ = false;
	}
	else
	{
		ans = x;
		for(size_t row = 0;  row < ans.rows_; row++)
			for(size_t column = 0; column < ans.columns_; column++)
			{
				ans.matrix[row][column] = x.matrix[row][column] - y.matrix[row][column];
			}
	}

	return ans;
}

Matrix operator * (const Matrix& x,const Matrix& y)
{
	Matrix ans;
	if(x.columns_ != y.rows_)
	{
		ans.available_ = false;
	}
	else
	{
		ans = Matrix(x.rows_,y.columns_);
		for(size_t xrow = 0; xrow < ans.rows_; xrow++)
			for(size_t ycolumn = 0; ycolumn < ans.columns_; ycolumn++)
			{
				double sum = 0;
				for(size_t element = 0; element < y.rows_; element++)
				{
					sum += x.matrix[xrow][element] * y.matrix[element][ycolumn];
				}
				ans.matrix[xrow][ycolumn] = sum;
			}
	}

	return ans;
}

