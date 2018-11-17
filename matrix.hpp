#ifndef SJTU_MATRIX_HPP
#define SJTU_MATRIX_HPP

#include <cstddef>
#include <initializer_list>
#include <utility>
#include <iterator>
#include <stdexcept>

using std::size_t;

namespace sjtu
{
template<class T>
class Matrix
{
private:
	size_t rowNum, colNum;
	T* mat;

	T* operator[](size_t p)
	{
		return mat + p * colNum;
	}

    const T* operator[] (size_t p) const
    {
        return mat + p * colNum;
    }
    
public:
	Matrix() : rowNum(0), colNum(0), mat(nullptr) {}

	Matrix(size_t n, size_t m, T _init = T()) : rowNum(n), colNum(m), mat(nullptr)
	{
        if (rowNum * colNum) mat = new T[rowNum * colNum];
		for (size_t i = 0; i < rowNum; i++)
			for (size_t j = 0; j < colNum; j++)
			{
				(*this)[i][j] = _init;
			}
	}

	explicit Matrix(std::pair <size_t, size_t> sz, T _init = T()) : Matrix(sz.first, sz.second, _init) {}

	Matrix(std::initializer_list <std::initializer_list<T>> il) : rowNum(il.size()), colNum(il.begin() -> size()), mat(nullptr)
	{
        if (rowNum * colNum) mat = new T[rowNum * colNum];
        auto curRow = il.begin();
		for (size_t i = 0; i < rowNum; i++, curRow++)
        {
        	if (curRow -> size() != colNum)
			{
        		delete[] mat;
        		throw(std::invalid_argument("initializer list has invalid size"));
			}
            auto curElement = curRow -> begin();
			for (size_t j = 0; j < colNum; j++, curElement++)
			{
				(*this)[i][j] = *curElement;
			}
        }
	}

	Matrix(const Matrix &o) : rowNum(o.rowLength()), colNum(o.columnLength()), mat(nullptr)
	{
        if (rowNum * colNum) mat = new T[rowNum * colNum];
		for (size_t i = 0; i < rowNum; i++)
			for (size_t j = 0; j < colNum; j++)
			{
				(*this)[i][j] = o[i][j];
			}
	}

	template<class U>
	Matrix(const Matrix<U> &o) : rowNum(o.rowLength()), colNum(o.columnLength()), mat(nullptr)
	{
        if (rowNum * colNum) mat = new T[rowNum * colNum];
        for (size_t i = 0; i < rowNum; i++)
            for (size_t j = 0; j < colNum; j++)
            {
				(*this)[i][j] = static_cast<T>(o(i,j));
            }
	}

	Matrix(Matrix &&o) noexcept : rowNum(o.rowLength()), colNum(o.columnLength()), mat(o.mat)
	{
		o.mat = nullptr;
	}


	Matrix &operator=(const Matrix &o)
	{
		if (mat == o.mat) return (*this);
        if (mat != nullptr) delete[] mat;

		rowNum = o.rowLength();
		colNum = o.columnLength();
		mat = new T[rowNum * colNum];
		for (size_t i = 0; i < rowNum; i++)
			for (size_t j = 0; j < colNum; j++)
			{
				(*this)[i][j] = o(i,j);
			}
		return (*this);
	}

	template<class U>
	Matrix &operator=(const Matrix<U> &o)
	{
        if (mat != nullptr) delete[] mat;
        rowNum = o.rowLength();
	    colNum = o.columnLength();
        mat = (rowNum * colNum) ? (new T[rowNum * colNum]) : nullptr;

        for (size_t i = 0; i < rowNum; i++)
            for (size_t j = 0; j < colNum; j++)
            {
                (*this)[i][j] = static_cast<T>(o(i, j));
            }
        return (*this);
	}

	Matrix &operator=(Matrix &&o) noexcept
	{
	    if (mat == o.mat) return *this;
		if (mat != nullptr) delete[] mat;
		rowNum = o.rowNum; colNum = o.colNum;
		mat = o.mat;
		o.mat = nullptr;
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

	std::pair<size_t, size_t> size() const
	{
		return std::make_pair(rowNum, colNum);
	}

	void resize(size_t n, size_t m, T _init = T())
	{
        if (rowNum * colNum != n * m)
        {
            T *newmat = (n * m) ? (new T[n * m]) : nullptr;
            for (size_t i = 0; i < std::min(n * m, rowNum * colNum); i++) newmat[i] = mat[i];
            for (size_t i = rowNum * colNum; i < n * m; i++) newmat[i] = _init;
            if (mat != nullptr) delete[] mat;
            mat = newmat;
        }
        rowNum = n;
        colNum = m;
	}

	void resize(std::pair <size_t, size_t> sz, T _init = T())
	{
        resize(sz.first, sz.second, _init);
	}

	void clear()
	{
        if (mat != nullptr) delete[] mat;
        rowNum = colNum = 0;
	}

public:
	const T &operator()(size_t i, size_t j) const
	{
		if (i < 0 || i >= rowNum || j < 0 || j >= colNum) throw(std::invalid_argument("operator() : invalid index"));
        return (*this)[i][j];
	}

	T &operator()(size_t i, size_t j)
	{
		if (i < 0 || i >= rowNum || j < 0 || j >= colNum) throw(std::invalid_argument("operator() : invalid index"));
		return (*this)[i][j];
	}

	Matrix row(size_t i) const
	{
		if (i < 0 || i >= rowNum) throw(std::invalid_argument("row : invalid index"));
		Matrix ret(1,colNum);
        for (size_t j = 0; j < colNum; j++) ret[0][j] = (*this)[i][j];
        return ret;
	}

	Matrix column(size_t i) const
	{
		if (i < 0 || i >= colNum) throw(std::invalid_argument("column : invalid index"));
		Matrix ret(rowNum,1);
        for (size_t j = 0; j < rowNum; j++) ret[j][0] = (*this)[j][i];
        return ret;
	}


public:
	template<class U>
	bool operator==(const Matrix<U> &o) const
	{
        if (rowNum != o.rowLength() || colNum != o.columnLength()) return false;
        for (size_t i = 0; i < rowNum; i++)
            for (size_t j = 0; j < colNum; j++)
                if ((*this)[i][j] != o(i,j)) return false;
        return true;
	}

	template<class U>
	bool operator!=(const Matrix<U> &o) const
	{
        return !((*this) == o);
	}

	Matrix operator-() const
	{
        Matrix ret(rowNum, colNum);
        for (size_t i = 0; i < rowNum; i++)
            for (size_t j = 0; j < colNum; j++) ret[i][j] = -(*this)[i][j];
        return ret;
	}

	template<class U>
	Matrix &operator+=(const Matrix<U> &o)
	{
		if (rowNum != o.rowLength() || colNum != o.columnLength()) throw(std::invalid_argument("operator += : invalid addend"));
		for (size_t i = 0; i < rowNum; i++)
            for (size_t j = 0; j < colNum; j++) (*this)[i][j] = static_cast<T>((*this)[i][j] + o(i,j));
        return (*this);
	}

	template<class U>
	Matrix &operator-=(const Matrix<U> &o)
	{
		if (rowNum != o.rowLength() || colNum != o.columnLength()) throw(std::invalid_argument("operator += : invalid subtrahend"));
		for (size_t i = 0; i < rowNum; i++)
            for (size_t j = 0; j < colNum; j++) (*this)[i][j] = static_cast<T>((*this)[i][j] - o(i,j));
        return (*this);
	}

	template<class U>
	Matrix &operator*=(const U &x)
	{
        for (size_t i = 0; i < rowNum; i++)
            for (size_t j = 0; j < colNum; j++) (*this)[i][j] = static_cast<T>((*this)[i][j] * x);
        return (*this);
	}

	Matrix tran() const
	{
        Matrix ret(colNum, rowNum);
        for (size_t i = 0; i < rowNum; i++)
            for (size_t j = 0; j < colNum; j++) ret[j][i] = (*this)[i][j];
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

		iterator(const iterator &o) = default;

		iterator &operator=(const iterator &o) = default;

		iterator(size_type _curRow, size_type _curCol, Matrix* _corresMatrix, pointer _ptr, std::pair <size_type, size_type> _subL = {0, 0}, std::pair <size_type, size_type> _subR = {0, 0})
            : curRow(_curRow), curCol(_curCol), corresMatrix(_corresMatrix), ptr(_ptr), subL(_subL), subR(_subR)
		{
			if (subR == std::make_pair((size_type)0, (size_type)0)) subR = std::make_pair(corresMatrix -> rowLength() - 1, corresMatrix -> columnLength() - 1);
		}

	private:
		size_type curRow, curCol;
		Matrix* corresMatrix;
        pointer ptr;
		std::pair <size_type, size_type> subL, subR;

		difference_type index() const
		{
			return (curRow - subL.first) * (subR.second - subL.second + 1) + (curCol - subL.second);
		}

	public:
		difference_type operator-(const iterator &o)
		{
			return index() - o.index();
		}

		iterator operator+(difference_type offset) const
		{
            size_type colLen = subR.second - subL.second + 1;
            size_type tmp = (curCol + offset) / colLen + ((offset < 0) && (curCol + offset) % colLen != 0);
            auto p = ptr;
            p -= curRow  * (corresMatrix -> columnLength()) + curCol + (subL.first * corresMatrix -> columnLength() + subL.second);
            p += (curRow + tmp) * (corresMatrix -> columnLength()) + (curCol + offset - tmp * colLen) + (subL.first * corresMatrix -> columnLength() + subL.second);
            return iterator(curRow + tmp, curCol + offset - tmp * colLen, corresMatrix, p, subL, subR);
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
        return iterator(0, 0, this, &(*this)[0][0]);
	}

	iterator end()
	{
        return iterator(rowNum -  1, colNum - 1, this, &(*this)[rowNum - 1][colNum - 1]);
	}

	std::pair <iterator, iterator> subMatrix(std::pair <size_t, size_t> l, std::pair <size_t, size_t> r)
	{
        return std::make_pair(iterator(0, 0, this, &(*this)[l.first][l.second], l, r), iterator(r.first - l.first, r.second - l.second, this, &(*this)[r.first][r.second], l, r));
	}
};

}

//
namespace sjtu
{
template<class T, class U>
auto operator*(const Matrix<T> &mat, const U &x)
{
    Matrix<decltype(T() * U())> ret(mat.rowLength(), mat.columnLength());
    for (size_t i = 0; i < ret.rowLength(); i++)
        for (size_t j = 0; j < ret.columnLength(); j++)
        	ret(i,j) = static_cast<decltype(T() * U())>(mat(i,j)) * static_cast<decltype(T() * U())>(x);
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
	if (a.columnLength() != b.rowLength()) throw(std::invalid_argument("operator * : a.col != b.row"));
    auto ret = Matrix<decltype(U() * V())>(a.rowLength(), b.columnLength());
    for (size_t i = 0; i < ret.rowLength(); i++)
    {
        for (size_t j = 0; j < ret.columnLength(); j++)
        {
            for (size_t k = 0; k < a.columnLength(); k++)
            {
                ret(i,j) += static_cast<decltype(U() * V())>(a(i,k)) * static_cast<decltype(U() * V())>(b(k,j));

			}
        }
    }
    return ret;
}

template<class U, class V>
auto operator+(const Matrix<U> &a, const Matrix<V> &b)
{
	if (a.rowLength() != b.rowLength() || a.columnLength() != b.columnLength()) throw(std::invalid_argument("operator + : invalid size"));
    auto ret = Matrix<decltype(U() + V())>(a.rowLength(), a.columnLength());
    for (size_t i = 0; i < ret.rowLength(); i++)
    {
        for (size_t j = 0; j < ret.columnLength(); j++)
			ret(i,j) = static_cast<decltype(U() * V())>(a(i,j)) + static_cast<decltype(U() * V())>(b(i,j));
    }
    return ret;
}

template<class U, class V>
auto operator-(const Matrix<U> &a, const Matrix<V> &b)
{
	if (a.rowLength() != b.rowLength() || a.columnLength() != b.columnLength()) throw(std::invalid_argument("operator - : invalid size"));
	return a + (-b);
}

}

#endif //SJTU_MATRIX_HPP

