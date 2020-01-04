# Efficient use of linear algebra library for high performance computing

## Background and Introduction

I'm working on researches in the field of Statistics. Modern days, because of the evolve of applications with big data, the statistical models have to cater the problems with more and more covariates included. Many people use R to build and test their model. However, running simulations with a large amount of covariates may take several weeks or months. Therefore, we want to know whether the simple codes using base R is appropriate for such projects. If not, we want to code in an efficient way so that we can spend more time on the theory or model itself rather than switching between servers and test the model through intensive but less efficient simulations. Since many people who work in this field haven't studied in computer science before, the following article will guide you to install the math library accordingly and the codes can be served as a basic benchmarking tool to test whether you have done the correct configuration so that when you code your own project using C++ and use such libraries, it will run in the most efficient way.

## List of Contents

* Installation of required packages and the benchmarking tool
* Linear algebra operations for different libraries in c++. MKL v.s. OpenBLAS v.s. Eigen
* Linear algebra operations for different languages. c++ (MKL) v.s. Matlab (MKL) v.s. R (LAPACK)
* Linear algebra operations for different architechures. Intel KabyLake v.s. AMD Zen+
* Conclusion and suggestion

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

The MKL, OpenBLAS and Eigen libraries should be installed (linux version) in order to run the following benchmark. 

#### OpenBLAS installation

OpenBLAS can be cloned into your local directory by using `git`.
```
git clone "https://github.com/xianyi/OpenBLAS.git"
```
It can be installed by invoking `make` in the directory which contains `Makefile`. Make sure `gfortran` is installed and `gcc` is updated which is compatible for `c++14`. Use `make install` to copy your `include` and `lib` folder your preferred installation directory.
```
make install PREFIX=your_installation_directory
```

#### MKL installation

MKL can be downloaded from intel's official webset for free. After dowloading the whole package, the installation can be done by invoking the program `install_GUI`. The default installation directory is `/home/usr/intel`, the `include` and `lib` folder can be found in the subfolder `/intel64` for 64 bit processors.

#### Eigen installation

Similar to MKL, Eigen can also be downloaded from the official webset. Eigen can be used driectly without installation. After untar the archive, the `include` libraries are in the subforder `/Eigen`. Usually, the dense matrix manipulation libraries can be used by including `<Eigen/Dense>` as the header file in your C++ code.

### Installation of Benchmarking tool

Clone the files into your local directory.
```
git clone "https://github.com/heilokchow/math_lib_benchmark.git"
```
Change the value of variables for include and library directories `MKL_INCDIR`, `MKL_LIBDIR`, `OMP_INCDIR`, `EIGEN_INCDIR`, `OPENBLAS_INCDIR` and `OPENBLAS_LIBDIR` in the `Makefile` according to your installation path of MKL, OpenBLAS and Eigen libraries.

By invoking `make MODEL=x`, the benchmark can be done for different libraries.

 * **MODEL=1:** MKL
 * **MODEL=2:** OpenBLAS
 * **MODEL=3:** Eigen
 * **MODEL=4:** Eigen + MKL
 * **MODEL=5:** Eigen + OpenBLAS

Then, the benchmark can be run by `./a.out <n> <nrep>`, where `n` is the dimension of the matrix and `nrep` is the number of replications. You can use the command `export OMP_NUM_THREADS=k` before running the program to force MKL and OpenBLAS use k threads. An important note is that if your own program use parrallel programming techniques already, you should use single thread for these math libraries. Also, Eigen can use `EIGEN_DONT_PARALLELIZE` to disable parallelization.

## Benchmarking for different libraries (i5-8250u @ 3.2GHz)

Here, the benchmark is run on i5-8250u @ 3.2GHz. Only single threaded performance is considered. The comparison results is ploted using `summary.R` and the results are shown below. The y-axis is the real time (not cpu time) spent and x-axis is the dimension of the matrix. Most paper uses the term "GFlops" to measure the performance of different implementations or algorithms. However, in real-life applications from other fields such as Statistics, the actual time spent is more crucial. Therefore, the elapsed time is applied for benchmarking. The `gcc` flags are `-O2`, `-mfma`. The FMA instruction is more efficient than SSE or AVX2, if the code is not compiled with FMA instrution sets, the Eigen library will have a dramastic performance loss.

### DGEMM subroutine

The DGEMM subroutine from BLAS library which calculates the general matrix and matrix product. It is widely used in all kinds of scientific computation and the efficiency of matrix multiplication is of vital importance. Both five methods have similar efficiency.

![Imgur](https://i.imgur.com/DKLtvFT.jpg)

### DGESV subroutine

The DGESV subroutine solves the linear equations which is widly used in mathematical optimization problems. Here, Eigen library performs better than MKL and OpenBLAS. However, the `llt()` function is used here for Eigen and only positive definite matrix is tested here. If more general martix is evaluated in your application or you want more accuarcy, you may consider function `partialPivLu()` or `fullPivLu()` which will be slower. 

![Imgur](https://i.imgur.com/MdRJkpR.jpg)

### DPOTRF subroutine

The DPOTRF subroutine does the Cholesky factorization. The Cholesky factorization is the essential part for dependent multivariate normal distribution construction. Both five methods have similar efficiency.

![Imgur](https://i.imgur.com/QltmWGA.jpg)

## Benchmarking for different languages (i5-8250u @ 3.2GHz)

Here, we still consider the single threaded performance when running benchmarks. The benchmarking code is contained in `summary.R` and `matlab_r2019.mlx`. Note that the default run on Matlab uses multi-threaded MKL library. To force Matlab using single thread, parrallel tool box could be installed and one worker thread should be used during the benchmarking run. The x-axis is the dimension of the matrix and the y-axis is `t^(1/3)` in term of seconds since all subroutines have an `O(n^3)` computation time. For example, from the figure below, the actual elapsed time for R to perform a 2000 x 2000 matrix multiplication is 1.55^3 which is around 4 seconds.

![Imgur](https://i.imgur.com/ZSdCmRd.jpg)

 Matlab and Eigen performs similar and are about 10-50 times faster than R(base) for large matrix related operations. Note that the default BLAS and LAPACK library is used for R. To boost the matrix manipulations, RcppEigen packages or relinking the library to OpenBLAS could be considerred. But both two methods are not easy to be implemented, the previous one still requires coding in c++ and the later one requires a reinstallation of R which is inapplicable on most HPC servers. Therefore, R may not be a suitable choice for big data manipulation if such linear operations is needed.

## Benchmarking for different architectures 

### i5-8250u @ 3.2GHz v.s. r7-2700x @ 4.0GHz

Due to the advantages of Eigen from benchmark results above (efficient and portable), the Eigen library is adopt for benchmarking on different architectures. Here, Intel's KabyLake is compared against with AMD's Zen+. Since the actual real application's performance in researches is of vital important, both computers are run on auto boost with different speeds. The single thread performance is recorded below. y-axis is the actual time (in seconds) takes, the lower, the better.

![Imgur](https://i.imgur.com/jrKWfmI.jpg)

Running at a lower clock, Intel's cpu is 10-25% times faster than AMD's cpu when running the Eigen library. Therefore, Zen or Zen+ architechture is not a good choice for double precision floating point computation which is common in most researches. If dense matrix operations is an essentail part of your program, Intel's most architectures or Zen2 architecture is more suitable for these tasks.

## Conclusion

If only 100 covariates is included in your model, R is still a good choice since the coding process is much simplier than C++ and their are many existing packages for you to use. However, if you works on models with more than 1000 covariates, you may consider other choices. If a controllable parallelization is not needed, Matlab is a good choice. Otherwise, using Eigen library with correct configuration will help you do simulations more efficiently.
