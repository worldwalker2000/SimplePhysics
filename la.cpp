#include "la.hpp"

#include <cstdlib>
#include <cassert>
#include <iostream>

Vec::Vec(int len)
    : buf(len)
{
}

void Vec::Debug(const char* name)
{
    std::printf("%s (%ld)\n", name, Size());

    for (int i = 0; i < Size(); i++)
        std::printf("%0.15f, ", At(i));

    std::printf("\n\n");
}

double Vec::At(int i) const
{
    return buf[i];
}

double& Vec::At(int i)
{
    return buf[i];
}

Vec& Vec::operator+=(const Vec& rhs)
{
    assert(Size() == rhs.Size());

    #pragma omp simd
    for (int i = 0; i < Size(); i++)
        buf[i] += rhs.buf[i];

    return *this;
}

Vec operator+(Vec lhs, const Vec& rhs)
{
    assert(lhs.Size() == rhs.Size());

    #pragma omp simd
    for (int i = 0; i < lhs.Size(); i++)
        lhs.buf[i] += rhs.buf[i];

    return lhs;
}

Vec& Vec::operator*=(const double c)
{
    #pragma omp simd
    for (int i = 0; i < Size(); i++)
        buf[i] += c;

    return *this;
}

Vec operator*(Vec lhs, const double c)
{
    #pragma omp simd
    for (int i = 0; i < lhs.Size(); i++)
        lhs.buf[i] *= c;

    return lhs;
}

Vec operator*(const double c, Vec rhs)
{
    #pragma omp simd
    for (int i = 0; i < rhs.Size(); i++)
        rhs.buf[i] *= c;

    return rhs;
}

Vec& Vec::operator*=(const Vec& rhs)
{
    assert(Size() == rhs.Size());

    #pragma omp simd
    for (int i = 0; i < Size(); i++)
        buf[i] *= rhs.buf[i];

    return *this;

}

Vec operator*(Vec lhs, const Vec& rhs)
{
    assert(lhs.Size() == rhs.Size());

    #pragma omp simd
    for (int i = 0; i < lhs.Size(); i++)
        lhs.buf[i] *= rhs.buf[i];

    return lhs;

}

void Vec::Zero()
{
    #pragma omp simd
    for (int i = 0; i < Size(); i++)
        At(i) = 0.0;
}

Mat::Mat(int r_, int c_)
    : r(r_), c(c_), buf(r_*c_)
{
}

void Mat::Debug(const char* name)
{
    std::printf("%s (%ld x %ld)\n", name, Rows(), Cols());

    for (int i = 0; i < Rows(); i++)
    {
        for (int j = 0; j < Cols(); j++)
            std::printf("%0.2f ", At(i, j));
        std::printf("\n");
    }

    std::printf("\n");
}


double Mat::At(int row, int col) const
{
    return buf[col + c * row];
}

double& Mat::At(int row, int col)
{
    return buf[col + c * row];
}

Vec operator*(Mat mat, const Vec& vec)
{
    assert(mat.c == vec.Size());

    Vec res(mat.r);

    for (int i = 0; i < mat.r; i++) // y
    {
        double sum = 0.0;
        for (int j = 0; j < mat.c; j++) // x
        {
            sum += mat.At(i, j) * vec.At(j);
        }

        res.At(i) = sum;
    }

    return res;
}

// TODO: parallel and rhs by ref
Mat operator*(const Mat& lhs, Mat rhs)
{
    assert(lhs.c == rhs.r);

    Mat mat(lhs.r, rhs.c);

    rhs.Transpose();

    int inner = lhs.c;

    for (int i = 0; i < mat.r; i++) // y
    {
        for (int j = 0; j < mat.c; j++) // x
        {
            double sum = 0.0;

            for (int k = 0; k < inner; k++)
                sum += lhs.At(i, k) * rhs.At(j, k); // rhs is transposed

            mat.At(i, j) = sum;
        }
    }

    return mat;
}

void Mat::Transpose()
{
    Mat mat(c, r);

    for (int i = 0; i < r; i++) // y
        for (int j = 0; j < c; j++) // x
            mat.At(j, i) = At(i, j);

    buf = mat.buf;
    r = mat.r;
    c = mat.c;
}

void Mat::Zero()
{
    for (int i = 0; i < Rows()*Cols(); i++)
        buf[i] = 0.0;
}

static void SwapRows(std::vector<double>& A, int n, int i, int j)
{
    for (int k = 0; k < n; k++)
    {
        double tmp = A[k + n * i];
        A[k + n * i] = A[k + n * j];
        A[k + n * j] = tmp;
    }
}

// TODO: optimize, use the bufs directly instread of copyies
Vec Mat::Solve(Mat mat, Vec bvec)
{
    assert(mat.Rows() == mat.Cols());
    assert(mat.Rows() == bvec.Size());

    std::vector<double>& A = mat.buf;
    std::vector<double>& b = bvec.buf;

    int n = mat.Rows();

    std::vector<std::pair<int, int>> swaps;

    for (int row = 0; row < n; row++)
    {
#if 1
        int bestRow = row;
        float value = A[row + n * bestRow];
        for (int r = row + 1; r < n; r++)
        {
            if (std::abs(A[row + n * r]) > value)
            {
                bestRow = r;
                value = A[row + n * r];
            }
        }

        // if (A[row + n * bestRow] == 0) assert(false && "Matrix is singular");
        // Just solve as 0 if non singular
        if (A[row + n * bestRow] == 0)
        {
            Vec xvec(n);
            // TODO: I think it is already zeroed
            xvec.Zero();
            std::printf("WARNING: Matrix is singular\n");

            return xvec;
        }

        if (bestRow != row)
        {
            SwapRows(A, n, row, bestRow);
            SwapRows(b, 1, row, bestRow);
            swaps.push_back(std::make_pair(row, bestRow));
        }
#endif

        for (int r = row + 1; r < n; r++)
        {
            float factor = A[row + n * r] / A[row + n * row];

            for (int c = 0; c < n; c++)
            {
                A[c + n * r] -= A[c + n * row] * factor;
            }

            b[r] -= b[row] * factor;
        }
    }

    std::vector<double> x(n);

    for (int r = n - 1; r >= 0; r--)
    {
        float v = b[r];
        for (int i = r + 1; i < n; i++)
        {
            v -= A[i + n * r] * x[i];
        }

        x[r] = v / A[r * n + r];
    }

    for (int i = swaps.size() - 1; i >= 0; i--)
    {
        SwapRows(x, 1, swaps[i].first, swaps[i].second);
    }

    Vec xvec(n);
    for (int i = 0; i < n; i++)
        xvec.At(i) = x[i];

    return xvec;
}
