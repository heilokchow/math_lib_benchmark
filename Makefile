MKL_INCDIR := -I/home/heilokchow/intel/mkl/include/
MKL_LIBDIR := /home/heilokchow/intel/mkl/lib/intel64/
OMP_LIBDIR := /home/heilokchow/intel/lib/intel64/
EIGEN_INCDIR := -I/home/heilokchow/Documents/external_package/EIGENDIR/
OPENBLAS_INCDIR := -I/home/heilokchow/Documents/Course/OpenBLAS/include/
OPENBLAS_LIBDIR := -L/home/heilokchow/Documents/Course/OpenBLAS/lib/
MKL_LIB := -Wl,--start-group $(MKL_LIBDIR)libmkl_intel_lp64.a $(MKL_LIBDIR)libmkl_gnu_thread.a $(MKL_LIBDIR)libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
LIB := -l:libopenblas.a -pthread
CCFLAG := -O2 -std=c++14
TARGETS := a.out
OBJECTS := bench.o
CC := g++

$(TARGETS) : $(OBJECTS)
	$(CC) -o $@ $< $(MKL_LIB)

bench.o : bench.cpp
	$(CC) -c $(DEPENDENCE) $< $(CCFLAG) $(EIGEN_INCDIR) $(MKL_INCDIR)

.PHONY : clean
clean :
	-rm $(TARGETS) $(OBJECTS)
