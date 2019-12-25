# MKL v.s. OpenBLAS v.s. Eigen

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

Similar to MKL, Eigen can be downloaded from the official webset. After untar the archive, the `include` libraries are in the subforder `/Eigen`. Usually, the dense matrix manipulation libraries can be used by including `<Eigen/Dense>` as the header file.

### Installation

Clone the files into your local directory.
```
git clone "https://github.com/heilokchow/math_lib_benchmark.git"
```
Change the value of variable for include and library directories `MKL_INCDIR`, `MKL_LIBDIR`, `OMP_INCDIR`, `EIGEN_INCDIR`, `OPENBLAS_INCDIR` and `OPENBLAS_LIBDIR` in the Makefile according to your installation path of MKL, OpenBLAS and Eigen libraries.

By invoking `make MODEL=x`, the benchmark can be done for different libraries.

 * **MODEL=1:** MKL
 * **MODEL=2:** OpenBLAS
 * **MODEL=3:** Eigen
 * **MODEL=4:** Eigen + MKL
 * **MODEL=5:** Eigen + OpenBLAS

Then, the benchmark can be run by `./a.out <n> <nrep>`, where `n` is the dimension of the matrix and `nrep` is the number of replications. 

## Benchmarking

### DGEMM subroutine

Here, the benchmark is run on i5-8250u @ 3.2GHz. 4 cores are registered for parallizing purpose. The DGEMM subroutine is applied for the benchmark. The comparison results is ploted using `summary.R` and the result is shown below. 

#### Multithreading Performance (8 threads)

![Imgur](https://i.imgur.com/C1GSGMg.jpg)

## Conclusion
From the result above, the MKL performs best on intel's KabyLake CPU. However, the result from OpenBLAS is quite close to MKL's performance. Also. OpenBLAS is designed to be optimized on different CPU architectures. Therefore, the performance might be in favour of OpenBLAS when running on Zen2 architecture. The performance of other subroutines such as matrix inversion and LU factorization on different kinds of architectures will be updated later. Also, as for Eigen library, the speed is similar to other libraries if the backends use OpenBLAS or MKL. Therefore, considering the fact that Eigen is much earier to use comparing to other libraries with only a little trade off on performance, Eigen might be the best linear algebra library when writing C++ codes for real applications.
