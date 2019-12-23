# MKL v.s. OpenBLAS v.s. Eigen

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

The MKL, OpenBLAS and Eigen libraries should be installed (linux version) in order to run the following benchmark

### Installing

Change the include and library directories in the Makefile according to your installation path of MKL, OpenBLAS and Eigen libraries.

### Benchmarking

Here, the benchmark is run on i5-8250u @ 3.2GHz. 4 cores are registered for parallizing purpose. The DGEMM subroutine is applied for the benchmark. The comparison results is ploted using `summary.R` and the result is shown below. 

![Imgur](https://i.imgur.com/C1GSGMg.jpg)

## Conclusion
From the result above, the MKL performs best on intel's KabyLake CPU. However, the result from OpenBLAS is quite close to MKL's performance. Also. OpenBLAS is designed to be optimized on different CPU architectures. Therefore, the performance might be in favour of OpenBLAS when running on Zen2 architecture. The performance of other subroutines such as matrix inversion and LU factorization on different kinds of architectures will be updated later. Also, as for Eigen library, the speed is similar to other libraries if the backends use OpenBLAS or MKL. Therefore, considering the fact that Eigen is much earier to use comparing to other libraries with only a little trade off on performance, Eigen might be the best linear algebra library when writing C++ codes for real applications.
