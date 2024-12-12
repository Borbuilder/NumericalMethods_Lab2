#include "SolveMatrix.h"

//std::vector<double> SolveMatrix(const std::vector<double>& A,
//    const std::vector<double>& C,
//    const std::vector<double>& B,
//    const std::vector<double>& Fi,
//    const std::vector<double>& ae,
//    const std::vector<double>& mu)
//{
//    const unsigned long long n = Fi.size() + 2;
//    if (A.size() != n - 2 || C.size() != n - 2 || B.size() != n - 2 || ae.size() != 2 || mu.size() != 2)
//        throw std::logic_error("Wrong size of vectors!");
//
//    std::vector<double> results(n);
//    std::vector<double> alpha(n - 1);
//    std::vector<double> beta(n - 1);
//
//    alpha[0] = ae[0];
//    beta[0] = mu[0];
//
//    for (unsigned long long i = 0; i < n - 2; ++i)
//    {
//        alpha[i + 1] = B[i] / (C[i] - alpha[i] * A[i]);
//        beta[i + 1] = (Fi[i] + A[i] * beta[i]) / (C[i] - alpha[i] * A[i]);
//    }
//
//    results[n - 1] = (mu[1] + beta[n - 2] * ae[1]) / (1 - alpha[n - 2] * ae[1]);
//
//    for (long long i = n - 2; i >= 0; i--)
//    {
//        results[i] = alpha[i] * results[i + 1] + beta[i];
//    }
//
//    return results;
//}

std::vector<double> SolveMatrix(const std::vector<std::vector<double>>& matrix_A, const std::vector<double>& vector_b)
{
    int size = vector_b.size();
    std::vector<double> y(size);
    std::vector<double> alpha(size, 0);
    std::vector<double> betta(size, 0);
    std::vector<double> A(size, 0);
    std::vector<double> B(size, 0);
    std::vector<double> C(size, 0);

    // Заполняем матрицу A и вектор B
    for (int i = 1; i < size - 1; ++i) {
        A[i] = matrix_A[i][0];
        C[i] = -matrix_A[i][1];
        B[i] = matrix_A[i][2];
    }

    alpha[1] = 0;
    betta[1] = vector_b[0];
    for (int i = 1; i < size - 1; ++i) {
        alpha[i + 1] = B[i] / (C[i] - A[i] * alpha[i]);
        betta[i + 1] = (-vector_b[i] + A[i] * betta[i]) / (C[i] - A[i] * alpha[i]);
    }

    // Начальный элемент вектора y
    y[size - 1] = vector_b[size - 1];

    // Обратный ход
    for (int i = size - 2; i >= 0; --i) {
        y[i] = alpha[i + 1] * y[i + 1] + betta[i + 1];
    }

    return y;
}
