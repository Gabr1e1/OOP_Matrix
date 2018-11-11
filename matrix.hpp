#ifndef SJTU_MATRIX_HPP
#define SJTU_MATRIX_HPP

#include <cstddef>
#include <initializer_list>
#include <utility>
#include <iterator>

using std::size_t;

namespace sjtu
{
template<class T>
class Matrix
{
private:
	T *mat;
	size_t row, col;

	T *operator[](size_t p)
	{
		return mat + p * m;
	}

public:
	Matrix() = default;

	Matrix(size_t n, size_t m, T _init = T()) : mat(new T[n * m]), row(n), col(m)
	{
		for (size_t i = 0; i < row; i++)
			for (size_t j = 0; j < col; j++)
			{
				mat[i][j] = _init;
			}
	}

	explicit Matrix(std::pair <size_t, size_t> sz, T _init = T()) : mat(new T[sz.first * sz.second]), row(sz.first),
																	col(sz.second)
	{
		for (size_t i = 0; i < row; i++)
			for (size_t j = 0; j < col; j++)
			{
				mat[i][j] = _init;
			}
	}

	Matrix(std::initializer_list <std::initializer_list<T>> il) : mat(new T[il[0].size() * il.size()]),
																  row(il[0].size()), col(il.size())
	{
		for (size_t i = 0; i < row; i++)
			for (size_t j = 0; j < col; j++)
			{
				mat[i][j] = il[i][j];
			}
	}

	Matrix(const Matrix &o) : mat(new T[o.rowLength() * o.colLength()]), row(o.rowLength()), col(o.colLength())
	{
		for (size_t i = 0; i < row; i++)
			for (size_t j = 0; j < col; j++)
			{
				mat[i][j] = o.mat[i][j];
			}
	}

	template<class U>
	Matrix(const Matrix<U> &o)
	{

	}

	Matrix &operator=(const Matrix &o)
	{

	}

	template<class U>
	Matrix &operator=(const Matrix<U> &o)
	{

	}

	Matrix(Matrix &&o) noexcept
	{

	}

	Matrix &operator=(Matrix &&o) noexcept
	{

	}

	~Matrix()
	{}

public:
	size_t rowLength() const
	{
		return row;
	}

	size_t columnLength() const
	{
		return col;
	}

	void resize(size_t _n, size_t _m, T _init = T())
	{

	}

	void resize(std::pair <size_t, size_t> sz, T _init = T())
	{

	}

	std::pair <size_t, size_t> size() const
	{

	}

	void clear()
	{

	}

public:
	const T &operator()(size_t i, size_t j) const
	{

	}

	T &operator()(size_t i, size_t j)
	{

	}

	Matrix<T> row(size_t i) const
	{

	}

	Matrix<T> column(size_t i) const
	{

	}


public:
	template<class U>
	bool operator==(const Matrix<U> &o) const
	{

	}

	template<class U>
	bool operator!=(const Matrix<U> &o) const
	{

	}

	Matrix operator-() const
	{

	}

	template<class U>
	Matrix &operator+=(const Matrix<U> &o)
	{

	}

	template<class U>
	Matrix &operator-=(const Matrix<U> &o)
	{

	}

	template<class U>
	Matrix &operator*=(const U &x)
	{

	}

	Matrix tran() const
	{

	}

public: // iterator
	class iterator
	{
	public:
		using iterator_category = std::random_access_iterator_tag;
		using value_type        = T;
		using pointer           = T *;
		using reference         = T &;
		using size_type         = size_t;
		using difference_type   = std::ptrdiff_t;

		iterator() = default;

		iterator(const iterator &) = default;

		iterator &operator=(const iterator &) = default;

		iterator(const int curn, const int curm, Matrix *mat,
				 std::pair <size_t, size_t> _subL = {0, 0}, std::pair <size_t, size_t> _subR = {0, 0}, pointer p)
				: curRow(curn), curCol(curm), subL(_subL), subR(_subR), corresMatrix(mat), ptr(p)
		{
			if (subR == {0, 0}) subR = {n - 1, m - 1};
		}

	private:
		pointer ptr;
		size_t curRow, curCol;
		std::pair <size_t, size_t> subL, subR;
		Matrix *corresMatrix;

		difference_type index()
		{
			return curRow * (subR - subL + 1) + curCol;
		}

	public:
		difference_type operator-(const iterator &o)
		{
			return index() - o.index();
		}

		iterator operator+(difference_type offset) const
		{
			auto tPtr = p - curRow * corresMatrix.col + curCol;
			curCol += offset;
			if (curCol <= subR.second) return iterator(curRow, curCol, mat, subL, subR, ptr + offset);
			curRow += curCol / (subR.second - subL.second + 1);
			curCol %= subR.second - subL.second + 1;
			return iterator(curRow, curCol, mat, subL, subR, tPtr + curRow * corresMatrix.col + curCol)
		}

		iterator operator-(difference_type offset) const
		{
			auto tPtr = p - curRow * corresMatrx.col + curCol;
			curCol -= offset;
			if (curCol >= subL.second) return iterator(curRow, curCol, mat, subL, subR, ptr - offset);
			curRow -= curCol
		}

		iterator &operator+=(difference_type offset)
		{
			*this = (*this) + offset;
			return *this;
		}

		iterator &operator-=(difference_type offset)
		{
			*this = (*this) - offset;
			return *this;
		}

		iterator &operator++() //++i
		{
			(*this) += 1;
			return *this;
		}

		iterator operator++(int) //i++
		{
			iterator tmp = (*this);
			++(*this);
			return tmp;
		}

		iterator &operator--()
		{
			(*this) -= 1;
			return *this;
		}

		iterator operator--(int)
		{
			iterator tmp = (*this);
			--(*this);
			return tmp;
		}

		reference operator*() const
		{
			return *p;
		}

		pointer operator->() const
		{
			return p;
		}

		bool operator==(const iterator &o) const
		{
			return ptr == o.ptr && curRow == o.curRow && curCol == o.curCol
				   && subL == o.subL && subR == o.subR && corresMatrix == o.corresMatrix;
		}

		bool operator!=(const iterator &o) const
		{
			return !((*this) == o);
		}
	};

	iterator begin()
	{

	}

	iterator end()
	{

	}

	std::pair <iterator, iterator> subMatrix(std::pair <size_t, size_t> l, std::pair <size_t, size_t> r)
	{

	}
};

}

//
namespace sjtu
{
template<class T, class U>
auto operator*(const Matrix<T> &mat, const U &x)
{

}

template<class T, class U>
auto operator*(const U &x, const Matrix<T> &mat)
{

}

template<class U, class V>
auto operator*(const Matrix<U> &a, const Matrix<V> &b)
{

}

template<class U, class V>
auto operator+(const Matrix<U> &a, const Matrix<V> &b)
{

}

template<class U, class V>
auto operator-(const Matrix<U> &a, const Matrix<V> &b)
{

}

}

#endif //SJTU_MATRIX_HPP

