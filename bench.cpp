#include <iostream>
#include <chrono>
#include <stdio.h>
#include <random>
#include <string>
#include <fstream>

// #define CHECK 1
#ifndef CHECK
#define CHECK 0
#endif 

#ifndef MODEL
#define MODEL 1
#endif

#if MODEL == 1 // MKL
#define MODE 0
#include <mkl.h>
#endif

#if MODEL == 2 // OPENBLAS
#define MODE 0
#include <cblas.h>
#include <lapacke.h>
#endif

#if MODEL == 3 // EIGEN
#define MODE 1
#define EIGEN_NO_DEBUG
#include <cblas.h>
#include <lapacke.h>
#endif

#if MODEL == 4 // EIGEN + MKL
#define MODE 1
#define EIGEN_USE_MKL_ALL
#define EIGEN_NO_DEBUG
#include <mkl.h>
#endif

#if MODEL == 5 // EIGEN + OPENBLAS
#define MODE 1
#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#define EIGEN_NO_DEBUG
#include <cblas.h>
#include <lapacke.h>
#endif

#include <Eigen/Dense>

void F_GEMM(double* const&, double* const&, double* const&, const int&, const int&, const int&);
void F_GESV(double* const&, double* const&, double* const&, const int&, const int&, const int&);
void F_POTRF(double* const&, const int&, const int&, const int&);
void check(const int&);

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cout << "No. of input: " << argc << std::endl;
        puts("./a.out <dim> <nrep>");
        exit(0);
    }

    int n = std::stoi(argv[1]);
    int nrep = std::stoi(argv[2]);

    double* x = new double[n*n];
    double* y = new double[n*n];
    double* z = new double[n*n];
    double* v = new double[n];
    double* b = new double[n];
    int p = 0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            p = i * n + j;
            x[p] = 0.0;
            y[p] = 0.0;
            z[p] = 0.0;
        }
        v[i] = 0;
        b[i] = 0;
    }   

    if (CHECK == 1) {
        check(n);
    }

    // Perform Martix Multiplication
    F_GEMM(x, y, z, n, nrep, MODE);

    // Perform Matrix Inversion
    F_GESV(v, y, b, n, nrep, MODE);

    // Perform Cholesky Decomposition
    F_POTRF(y, n, nrep, MODE);

    delete[] x;
    delete[] y;
    delete[] z;
    delete[] v;
    delete[] b;
    return 0;
}

void F_GEMM(double* const& x, double* const& y, double* const& z, const int& n, const int& nrep, const int& mode) {
    std::ofstream ot;
    ot.open("result.txt", std::ios_base::app);
    double t = 0.0;
    int p = 0;
    auto t0 = std::chrono::system_clock::now();
    auto t1 = std::chrono::system_clock::now();
    std::random_device device;
    std::mt19937 generator(device());
    std::normal_distribution<double> normal(0.0, 1.0);
    for (size_t i = 0; i < nrep; i++)
    {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                p = i * n + j;
                y[p] = normal(generator);
                z[p] = normal(generator);
            }
        }
        if (mode == 1) {
            Eigen::MatrixXd X(n, n), Y(n, n), C(n, n);
            int p = 0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    p = i * n + j;
                    X(i, j) = y[p];
                    Y(i, j) = z[p];
                    C(i, j) = 0;
                }
            }
            t0 = std::chrono::system_clock::now();
            C.noalias() = X * Y;
            t1 = std::chrono::system_clock::now();
        }
        else {
            t0 = std::chrono::system_clock::now();
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, y, n, z, n, 0.0, x, n);
            t1 = std::chrono::system_clock::now();
        }
        std::chrono::duration<double> elapsed = t1 - t0;
        std::cout << "DGEMM elapsed time: " << elapsed.count() << "s\n";
        t += static_cast<double>(elapsed.count());
    }
    std::cout << "DGEMM average time: " << t / nrep << "s\n";
    ot << "DGEMM," << MODEL << "," << n << "," << t / nrep << "\n";
    ot.close();
}

void F_GESV(double* const& v, double* const& y, double* const& b, const int& n, const int& nrep, const int& mode) {
    std::ofstream ot;
    ot.open("result.txt", std::ios_base::app);
    double t = 0.0;
    double* y1 = new double[n * n];
    double* y2 = new double[n * n];
    int p = 0, q = 0;
    lapack_int ret = 0;
    auto t0 = std::chrono::system_clock::now();
    auto t1 = std::chrono::system_clock::now();
    std::random_device device;
    std::mt19937 generator(device());
    std::normal_distribution<double> normal(0.0, 1.0);
    for (size_t i = 0; i < nrep; i++)
    {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                p = i * n + j;
                q = j * n + i;
                y1[p] = normal(generator);
                y2[q] = y1[p];
            }
            b[i] = normal(generator);
        }
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, y1, n, y2, n, 0.0, y, n);

        if (mode == 1) {
            Eigen::MatrixXd X(n, n);
            Eigen::VectorXd C(n), Y(n);
            int p = 0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    p = i * n + j;
                    X(i, j) = y[p];
                }
                Y(i) = b[i];
                C(i) = 0;
            }
            t0 = std::chrono::system_clock::now();
            C = X.ldlt().solve(Y);
            t1 = std::chrono::system_clock::now();
        }
        else {
            int* ipiv = new int[n];
            t0 = std::chrono::system_clock::now();
            ret = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, 1, y, n, ipiv, b, 1);
            t1 = std::chrono::system_clock::now();
            delete[] ipiv;
        }
        std::chrono::duration<double> elapsed = t1 - t0;
        std::cout << "DGESV elapsed time: " << elapsed.count() << "s\n";
        t += static_cast<double>(elapsed.count());
    }
    std::cout << "DGESV average time: " << t / nrep << "s\n";
    ot << "DGESV," << MODEL << "," << n << "," << t / nrep << "\n";
    ot.close();
    delete[] y1;
    delete[] y2;
}

void F_POTRF(double* const& y, const int& n, const int& nrep, const int& mode) {
    std::ofstream ot;
    ot.open("result.txt", std::ios_base::app);
    double t = 0.0;
    int p = 0, q = 0;
    lapack_int ret = 0;
    double* y1 = new double[n * n];
    double* y2 = new double[n * n];
    auto t0 = std::chrono::system_clock::now();
    auto t1 = std::chrono::system_clock::now();
    std::random_device device;
    std::mt19937 generator(device());
    std::normal_distribution<double> normal(0.0, 1.0);
    for (size_t i = 0; i < nrep; i++)
    {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                p = i * n + j;
                q = j * n + i;
                y1[p] = normal(generator);
                y2[q] = y1[p];
            }
        }
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, y1, n, y2, n, 0.0, y, n);

        if (mode == 1) {
            Eigen::MatrixXd X(n, n);
            int p = 0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    p = i * n + j;
                    X(i, j) = y[p];
                }
            }
            t0 = std::chrono::system_clock::now();
            Eigen::MatrixXd T = X.llt().matrixL();
            t1 = std::chrono::system_clock::now();
        }
        else {
            int* ipiv = new int[n];
            t0 = std::chrono::system_clock::now();
            ret = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', n, y, n);
            t1 = std::chrono::system_clock::now();
            delete[] ipiv;
        }
        std::chrono::duration<double> elapsed = t1 - t0;
        std::cout << "DPOTRF elapsed time: " << elapsed.count() << "s\n";
        t += static_cast<double>(elapsed.count());
    }
    std::cout << "DPOTRF average time: " << t / nrep << "s\n";
    ot << "DPOTRF," << MODEL << "," << n << "," << t / nrep << "\n";
    ot.close();
    delete[] y1;
    delete[] y2;
}

void check(const int& n) {
    double* x = new double[n * n];
    double* y = new double[n * n];
    double* z = new double[n * n];
    double* y1 = new double[n * n];
    double* y2 = new double[n * n];
    double* y3 = new double[n * n];
    double* a = new double[n];
    double* b = new double[n];
    int* ipiv = new int[n];
    int p = 0, q = 0;
    std::random_device device;
    std::mt19937 generator(device());
    std::normal_distribution<double> normal(0.0, 1.0);
    double c1 = 0.0, c2 = 0.0, c3 = 0.0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            p = i * n + j;
            q = j * n + i;
            y1[p] = normal(generator);
            y2[q] = y1[p];
            z[p] = normal(generator);
            z[q] = z[p];
            x[p] = 0.0;
            x[q] = 0.0;
        }
        a[i] = normal(generator);
        b[i] = 0.0;
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, y1, n, y2, n, 0.0, y, n);

    Eigen::MatrixXd X(n, n), Y(n, n), Z(n, n);
    Eigen::VectorXd A(n), B(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            p = i * n + j;
            Y(i, j) = y[p];
            Z(i, j) = z[p];
            X(i, j) = 0.0;
        }
        A(i) = a[i];
        B(i) = b[i];
    }

    // DGEMM
    X.noalias() = Y * Z;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, y, n, z, n, 0.0, x, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            c1 = c1 + abs(x[i * n + j] - X(i, j));
        }
    }

    // DGESV
    B = Y * A;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, 1, n, 1.0, y, n, a, 1, 0.0, b, 1);
    
    A = Y.ldlt().solve(B);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            y3[i * n + j] = y[i * n + j];
        }
    }
    LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, 1, y, n, ipiv, b, 1);

    for (int i = 0; i < n; i++) {
        c2 = c2 + abs(b[i] - A(i));
    }

    //DPOTRF
    Eigen::MatrixXd T = Y.llt().matrixL();
    LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', n, y3, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < (i + 1); j++) {
           c3 = c3 + abs(T(i, j) - y3[i * n + j]);
        }
    }
    
    std::cout << c1 << "," << c2 << "," << c3 << '\n';

    delete[] x;
    delete[] y;
    delete[] z;
    delete[] y1;
    delete[] y2;
    delete[] y3;
    delete[] a;
    delete[] b;
    delete[] ipiv;
}