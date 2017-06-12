/*
	Author : PlanC
	E-mail : hubenchang0515@outlook.com
	Blog   : www.kurukurumi.com

	File   : Matrix.h
	Date   : 2017-4-22
*/

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <cstddef>
#include <string>
#include <iostream>
#include <cstdio>

template<typename T> using Vector = std::vector<T>;

class Matrix
{	
public:

	~Matrix() = default;
	Matrix() = default;

	Matrix(const std::initializer_list<std::initializer_list<double>>& matrix);
	Matrix(const Matrix& matrix);
	Matrix(Matrix&& matrix);

	explicit Matrix(const Vector<Vector<double>>& matrix);
	explicit Matrix(size_t rows, size_t columns);
	
	/* Create a row-vector Matrix. */
	static Matrix rVector(Vector<double> vector);
	/* Create a column-vector Matrix. */
	static Matrix cVector(Vector<double> vector);
	/* Create a diagonal Matrix. */
	static Matrix diagonal(Vector<double> vector);
	/* Create a identity Matrix.  */
	static Matrix identity(size_t power);

	/* Print Matrix's Elements */
	void print();

	/* Elementary transform of Matrix */
	/* swap (ri) and (rj) row */
	Matrix& rSwap(size_t ri, size_t rj);
	/* elements of (ri) row multiply (n) */
	Matrix& rMultiply(size_t ri, double n);
	/* elements of (ri) row add  */
	Matrix& rTransform(size_t ri, size_t rj, double n);
	/* column transform */
	Matrix& cSwap(size_t ci, size_t cj);
	Matrix& cMultiply(size_t ci, double n);
	Matrix& cTransform(size_t ci, size_t cj, double n);

	/* calculate deternimant of Matrix. */
	double determinant();
	/* calculate the inverse Matrix. */
	Matrix inverse();
	/* calculate the adjoint Matrix. */
	// Matrix adjoint();

	/* Change row or column to be a new Matrix. */
	/* remove a row */
	Matrix& rRemove(size_t row);
	/* insert a row */
	Matrix& rInsert(size_t row,const Matrix& rVector);
	Matrix& rInsert(size_t row,Vector<double> rVector);
	/* replace a row */
	Matrix& rReplace(size_t row,const Matrix& rVector);
	Matrix& rReplace(size_t row,Vector<double> rVector);
	/* remove a column */
	Matrix& cRemove(size_t column);
	/* insert a column */
	Matrix& cInsert(size_t column,const Matrix& cVector);
	Matrix& cInsert(size_t column,Vector<double> cVector);
	/* replace a column */
	Matrix& cReplace(size_t column,const Matrix& cVector);
	Matrix& cReplace(size_t column,Vector<double> cVector);

	/* get a column vector */
	Matrix column(size_t column);

	/* treat this as a extended matrix to calulate the solve matrix */
	Matrix cramer();

	inline Matrix& setName(std::string name)
	{
		this->name_ = name;
		return *this;
	}

	inline std::string name()
	{
		return name_;
	}

	inline size_t rows()
	{
		return rows_;
	}

	inline size_t columns()
	{
		return columns_;
	}

	inline bool isAvailable()
	{
		return available_;
	}

	inline bool isSuccessful()
	{
		return successful_;
	}

	Matrix& operator = (const Matrix& matrix);
	Matrix& operator = (Matrix&& matrix);
	Vector<double>& operator [] (size_t i);
	const Vector<double>& operator [] (size_t i) const;
	
	friend bool operator == (const Matrix& x, const Matrix& y);
	friend bool operator != (const Matrix& x, const Matrix& y);
	friend Matrix operator * (const Matrix& x, double n);
	friend Matrix operator * (double n, const Matrix& x);
	friend Matrix operator + (const Matrix& x, const Matrix& y);
	friend Matrix operator - (const Matrix& x, const Matrix& y);
	friend Matrix operator * (const Matrix& x, const Matrix& y);

private:
	bool available_ = true;
	bool successful_ = true;
	size_t rows_ = 0;
	size_t columns_ = 0;
	std::string name_ = "Matrix";
	Vector<Vector<double>> matrix;

	std::string errorString(std::string function,std::string message)
	{
		char address[64];
		sprintf(address,"%p",this);
		std::string result = "[" + name_ + "(" + address + ")::";
		result += function + "] : " + message;
		return result;
	}
};

#endif
