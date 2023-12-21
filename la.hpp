#pragma once

#include <vector>

struct Vec
{
    std::vector<double> buf;

    explicit Vec(int len);

    void Debug(const char* name);

    double At(int i) const;
    double& At(int i);

    Vec& operator+=(const Vec& rhs);
    friend Vec operator+(Vec lhs, const Vec& rhs);

    Vec& operator*=(const double c);
    friend Vec operator*(Vec lhs, const double c);
    friend Vec operator*(const double c, Vec rhs);

    Vec& operator*=(const Vec& rhs);
    friend Vec operator*(Vec lhs, const Vec& rhs);

    void Zero();

    inline std::size_t Size() const { return buf.size(); }
};

struct Mat
{
    int r, c;
    std::vector<double> buf;

    explicit Mat(int r_, int c_);

    void Debug(const char* name);

    double At(int row, int col) const;
    double& At(int row, int col);

    friend Vec operator*(Mat mat, const Vec& vec);

    friend Mat operator*(const Mat& lhs, Mat rhs);

    void Transpose();
    void Zero();

    inline std::size_t Rows() const { return r; }
    inline std::size_t Cols() const { return c; }

    static Vec Solve(Mat mat, Vec bvec);
};
