#include <iostream>
#include <chrono>
#include <stdio.h>
#include <random>
#include <string>
#include <fstream>

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
            for (size_t j = i; j < n; j++) {
                p = i * n + j;
                q = j * n + i;
                y[p] = normal(generator);
                y[q] = y[p];
            }
            b[i] = normal(generator);
        }
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
            ret = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, 1, y, n, ipiv, b, n);
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
}

void F_POTRF(double* const& y, const int& n, const int& nrep, const int& mode) {
    std::ofstream ot;
    ot.open("result.txt", std::ios_base::app);
    double t = 0.0;
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
            for (size_t j = i; j < n; j++) {
                p = i * n + j;
                q = j * n + i;
                y[p] = normal(generator);
                y[q] = y[p];
            }
        }
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
            auto T = X.llt().matrixL();
            t1 = std::chrono::system_clock::now();
        }
        else {
            int* ipiv = new int[n];
            t0 = std::chrono::system_clock::now();
            ret = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'U', n, y, n);
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
}