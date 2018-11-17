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
	size_t rowNum, colNum;

	T *operator[](size_t p)
	{
		return mat + p * colNum;
	}

public:
	Matrix() = default;

	Matrix(size_t n, size_t m, T _init = T()) : rowNum(n), colNum(m)
	{
        mat = new T[n * m];
		for (size_t i = 0; i < rowNum; i++)
			for (size_t j = 0; j < colNum; j++)
			{
				mat[i][j] = _init;
			}
	}

	explicit Matrix(std::pair <size_t, size_t> sz, T _init = T()) : Matrix(sz.first, sz.second, _init) { }

	Matrix(std::initializer_list <std::initializer_list<T>> il) : rowNum(il[0].size()), colNum(il.size())
	{
        mat = new T[il[0].size * il.size()];
		for (size_t i = 0; i < rowNum; i++)
			for (size_t j = 0; j < colNum; j++)
			{
				mat[i][j] = il[i][j];
			}
	}

	Matrix(const Matrix &o) : rowNum(o.rowNumLength()), colNum(o.colLength())
	{
        mat = new T[o.rowNumLength() * o.colLength()];
		for (size_t i = 0; i < rowNum; i++)
			for (size_t j = 0; j < colNum; j++)
			{
				mat[i][j] = o.mat[i][j];
			}
	}

	template<class U>
	Matrix(const Matrix<U> &o) : rowNum(o.rowNumLength()), colNum(o.colLength())
	{
        mat = new T[o.rowNum.Length() * o.colLength()];
        for (size_t i = 0; i < rowNum; i++)
            for (size_t j = 0; j < colNum; j++)
            {
                mat[i][j] = static_cast<T>(o.mat[i][j]);
            }
	}

	Matrix &operator=(const Matrix &o)
	{
        (*this) = Matrix(o);
        return (*this);
	}

	template<class U>
	Matrix &operator=(const Matrix<U> &o)
	{
        (*this) = Matrix(o);
        return (*this);
	}

	Matrix(Matrix &&o) noexcept : rowNum(o.rowLength()), colNum(o.colLength()), mat(o.mat)
	{
        o.mat = nullptr;
	}

	Matrix &operator=(Matrix &&o) noexcept
	{
        (*this) = Matrix(o);
        return *this;
	}

	~Matrix()
	{
        if (rowNum || colNum) delete[] mat;
    }

public:
	size_t rowLength() const
	{
		return rowNum;
	}

	size_t columnLength() const
	{
		return colNum;
	}

	void resize(size_t n, size_t m, T _init = T())
	{
        if (rowNum * colNum < n * m)
        {
            T *newmat = new T[n * m];
            for (int i = 0; i < rowNum * colNum; i++) newmat[i] = mat[i];
            for (int i = rowNum * colNum; i < n * m; i++) mat[i] = _init;
            delete[] mat;
            mat = newmat;
        }
        rowNum = n;
        colNum = m;
	}

	void resize(std::pair <size_t, size_t> sz, T _init = T())
	{
        resize(sz.first, sz.second, _init);
	}

	std::pair <size_t, size_t> size() const
	{
        return std::make_pair(rowNum, colNum);
	}

	void clear()
	{
        if (rowNum || colNum) delete[] mat;
        rowNum = colNum = 0;
	}

public:
	const T &operator()(size_t i, size_t j) const
	{
        return mat[i][j];
	}

	T &operator()(size_t i, size_t j)
	{
        return mat[i][j];
	}

	Matrix row(size_t i) const
	{
        Matrix ret = Matrix(1,colNum);
        for (int i = 0; i < colNum; i++) ret.mat[0][i] = mat[0][i];
        return ret;
	}

	Matrix column(size_t i) const
	{
        Matrix ret = Matrix(rowNum,1);
        for (int i = 0; i < rowNum; i++) ret.mat[i][0] = mat[i][0];
        return ret;
	}


public:
	template<class U>
	bool operator==(const Matrix<U> &o) const
	{
        if (rowNum != o.rowLength() || colNum != o.colLength()) return false;
        for (int i = 0; i < rowNum; i++)
            for (int j = 0; j < colNum; j++)
                if (mat[i][j] != o(i,j)) return false;
        return true;
	}

	template<class U>
	bool operator!=(const Matrix<U> &o) const
	{
        return !((*this) == o);
	}

	Matrix operator-() const
	{
        Matrix ret = Matrix(rowNum, colNum);
        for (int i = 0; i < rowNum; i++)
            for (int j = 0; j < colNum; j++) ret.mat[i][j] = -mat[i][j];
        return ret;
	}

	template<class U>
	Matrix &operator+=(const Matrix<U> &o)
	{
        for (int i = 0; i < rowNum; i++)
            for (int j = 0; j < colNum; j++) mat[i][j] += o(i,j);
        return (*this);
	}

	template<class U>
	Matrix &operator-=(const Matrix<U> &o)
	{
        for (int i = 0; i < rowNum; i++)
            for (int j = 0; j < colNum; j++) mat[i][j] -= o(i,j);
        return (*this);
	}

	template<class U>
	Matrix &operator*=(const U &x)
	{
        for (int i = 0; i < rowNum; i++)
            for (int j = 0; j < colNum; j++) mat[i][j] *= static_cast<T>(x);
        return (*this);
	}

	Matrix tran() const
	{
        Matrix ret = Matrix(colNum, rowNum);
        for (int i = 0; i < rowNum; i++)
            for (int j = 0; j < colNum; j++) ret.mat[j][i] = mat[i][j];
        return ret;
	}

public: // iterator
	class iterator
	{
	public:
		using iterator_category = std::random_access_iterator_tag;
		using value_type        = T;
		using pointer           = T*;
		using reference         = T&;
		using size_type         = size_t;
		using difference_type   = std::ptrdiff_t;

		iterator() = default;

		iterator(const iterator &) = default;

		iterator &operator=(const iterator &) = default;

		iterator(int _curRow, int _curCol, Matrix *mat, pointer p, std::pair <size_t, size_t> _subL = {0, 0}, std::pair <size_t, size_t> _subR = {0, 0})
				: curRow(_curRow), curCol(_curCol), subL(_subL), subR(_subR), corresMatrix(mat), ptr(p)
		{
			if (subR == std::make_pair((size_type)0, (size_type)0)) subR = std::make_pair(corresMatrix -> rowLength() - 1, corresMatrix -> colLength() - 1);
		}

	private:
		size_type curRow, curCol;
		Matrix *corresMatrix;
        pointer ptr;
		std::pair <size_t, size_t> subL, subR;

		difference_type index()
		{
			return curRow * (subR.second - subL.second + 1) + curCol;
		}

	public:
		difference_type operator-(const iterator &o)
		{
			return index() - o.index();
		}

		iterator operator+(difference_type offset) const
		{
            size_type colLen = subR.second - subL.second + 1;
            size_type tmp = (curCol + offset - subL.second) / colLen + ((offset < 0) && (curCol + offset - subL.second) % colLen != 0);
            auto p = ptr;
            p -= curRow * mat.colLength() + curCol;
            p += (curRow + tmp) * mat.colLength() + (curCol + offset - tmp * colLen);
            return iterator(curRow + tmp, curCol + offset - tmp * colLen, mat, subL, subR, ptr); 
		}

		iterator operator-(difference_type offset) const
		{
            return (*this) + (-offset);
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
			return *ptr;
		}

		pointer operator->() const
		{
			return ptr;
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
        return iterator(0, 0, this, &mat[0][0]);
	}

	iterator end()
	{
        return iterator(rowNum -  1, colNum - 1, this, &mat[rowNum - 1][colNum - 1]);
	}

	std::pair <iterator, iterator> subMatrix(std::pair <size_t, size_t> l, std::pair <size_t, size_t> r)
	{
        return std::make_pair(iterator(0, 0, this, &mat[l.first][l.second], l, r), iterator(r.first - l.first, r.second - l.second, this, &mat[r.first][r.second], l, r));
	}
};

}

//
namespace sjtu
{
template<class T, class U>
auto operator*(const Matrix<T> &mat, const U &x)
{
    Matrix<decltype(T() * U())> ret = mat;
    for (int i = 0; i < ret.rowLength(); i++)
        for (int j = 0; j < ret.colLength(); j++) ret(i,j) *= x;
    return ret;
}

template<class T, class U>
auto operator*(const U &x, const Matrix<T> &mat)
{
    return mat * x;
}

template<class U, class V>
auto operator*(const Matrix<U> &a, const Matrix<V> &b)
{
    auto ret = Matrix<decltype(U() * V())>(a.rowLength(), b.colLength());
    for (int i = 0; i < ret.rowLength(); i++)
    {
        for (int j = 0; j < ret.colLength(); j++)
        {
            for (int k = 0; k < a.colLength(); k++)
            {
                ret(i,j) += a(i,k) * b(k,j);
            }
        }
    }
    return ret;
}

template<class U, class V>
auto operator+(const Matrix<U> &a, const Matrix<V> &b)
{
    auto ret = Matrix<decltype(U() + V())>(a.rowLength(), a.colLength());
    for (int i = 0; i < ret.rowLength(); i++)
    {
        for (int j = 0; j < ret.colLength(); j++) ret(i,j) = a(i,j) + b(i,j);
    }
    return ret;
}

template<class U, class V>
auto operator-(const Matrix<U> &a, const Matrix<V> &b)
{
    return a + (-b);
}

}

#endif //SJTU_MATRIX_HPP

