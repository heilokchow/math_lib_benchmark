# MKL v.s. OpenBLAS v.s. Eigen

## List of Contents

* Installation of required packages and the benchmarking tool
* Linear algebra operations for different libraries in c++. MKL v.s. OpenBLAS v.s. Eigen
* Linear algebra operations for different languages. c++ (MKL) v.s. Matlab (MKL) v.s. R (LAPACK)
* Linear algebra operations for different architechures. Intel KabyLake v.s. AMD Zen+
* Analysis of the numercial errors
* Conclusion and suggestion

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

The MKL, OpenBLAS and Eigen libraries should be installed (linux version) in order to run the following benchmark. 

#### OpenBLAS installation

OpenBLAS can be cloned into your local directory and installed by invoking `make`. Make sure `gfortran` is installed and `gcc` is updated which is compatible for `c++14`.
```
git clone "https://github.com/heilokchow/OpenBLAS.git"
make
```
Use `make install` to copy your `include` and `lib` folder your preferred installation directory.
```
make install PREFIX=your_installation_directory
```

#### MKL installation

MKL can be downloaded from intel's official webset for free. After dowloading the whole package, the installation can be done by invoking the program `install_GUI`. The default installation directory is `/home/usr/intel`, the `include` and `lib` folder can be found in the subfolder `/intel64` for 64 bit processors.

#### Eigen installation

Similar to MKL, Eigen can also be downloaded from the official webset. After untar the archive, the `include` libraries are in the subforder `/Eigen`. Usually, the dense matrix manipulation libraries can be used by including `<Eigen/Dense>` as the header file.

### Installation

Clone the files into your local directory.
```
git clone "https://github.com/heilokchow/math_lib_benchmark.git"
```
Change the value of variables for include and library directories `MKL_INCDIR`, `MKL_LIBDIR`, `OMP_INCDIR`, `EIGEN_INCDIR`, `OPENBLAS_INCDIR` and `OPENBLAS_LIBDIR` in the Makefile according to your installation path of MKL, OpenBLAS and Eigen libraries.

By invoking `make MODEL=x`, the benchmark can be done for different libraries.

 * **MODEL=1:** MKL
 * **MODEL=2:** OpenBLAS
 * **MODEL=3:** Eigen
 * **MODEL=4:** Eigen + MKL
 * **MODEL=5:** Eigen + OpenBLAS

Then, the benchmark can be run by `./a.out <n> <nrep>`, where `n` is the dimension of the matrix and `nrep` is the number of replications. 

## Benchmarking for different libraries (i5-8250u @ 3.2GHz)

Here, the benchmark is run on i5-8250u @ 3.2GHz. Only single threaded performance is considered. The comparison results is ploted using `summary.R` and the results are shown below. The y-axis is the real time (not cpu time) spent and x-axis is the dimension of the matrix. Most paper uses the term "GFlops" to measure the performance of different implementations or algorithms. However, in real-life applications from other fields such as Physics and Statistics, the actual time spent is more crucial. Therefore, the elapsed time is applied for benchmarking.

### DGEMM subroutine

The DGEMM subroutine from BLAS library which calculates the general matrix and matrix product. It is widely used in all kinds of scientific computation and the efficiency of matrix multiplication is of vital importance. The following results reveals that other than the pure Eigen library, other libaries perform similarly.

![Imgur](https://i.imgur.com/DKLtvFT.jpg)

### DGESV subroutine

The DGESV subroutine solves the linear equations which is widly used in mathematical optimization problems. Since Eigen doesn't use the LAPACK's DGESV subroutine as its backend, the performance of Eigen+MKL and Eigen+OpenBLAS is not quite good and it is similar to the pure Eigen library's result.

![Imgur](https://i.imgur.com/MdRJkpR.jpg)

### DPOTRF subroutine

The DPOTRF subroutine does the Cholesky factorization. The Cholesky factorization is the essential part for positive definite matrix's inversion. The covariance matricies from statistical analysis belongs to such matrix class. The result of Cholesky factorization is similar to the result of matrix multiplication where 4 methods which applies LAPACK subroutines perform similar. 

![Imgur](https://i.imgur.com/QltmWGA.jpg)

