MKL_INCDIR := -I/home/heilokchow/intel/mkl/include/
MKL_LIBDIR := /home/heilokchow/intel/mkl/lib/intel64/
OMP_LIBDIR := /home/heilokchow/intel/lib/intel64/
EIGEN_INCDIR := -I/home/heilokchow/Documents/external_package/EIGENDIR/
OPENBLAS_INCDIR := -I/home/heilokchow/Documents/Course/OpenBLAS/include/
OPENBLAS_LIBDIR := /home/heilokchow/Documents/Course/OpenBLAS/lib/
MKL_LIB := -Wl,--start-group $(MKL_LIBDIR)libmkl_intel_lp64.a $(MKL_LIBDIR)libmkl_gnu_thread.a $(MKL_LIBDIR)libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
OPENBLAS_LIB := -L$(OPENBLAS_LIBDIR) -l:libopenblas.a -pthread
CCFLAG := -O2 -std=c++14
TARGETS := a.out
OBJECTS := bench.o
CC := g++

ifndef MODEL
MODEL = 1
endif

ifeq ($(MODEL), 1)
INC = $(MKL_INCDIR)
LIB = $(MKL_LIB)
else ifeq ($(MODEL), 2)
INC = $(OPENBLAS_INCDIR)
LIB = $(OPENBLAS_LIB)
else ifeq ($(MODEL), 3)
INC = $(EIGEN_INCDIR)
else ifeq ($(MODEL), 4)
INC = $(EIGEN_INCDIR) $(MKL_INCDIR)
LIB = $(MKL_LIB)
else
INC = $(EIGEN_INCDIR) $(OPENBLAS_INCDIR)
LIB = $(OPENBLAS_LIB)
endif

DEPENDENCE = -DMODEL=$(MODEL)

$(TARGETS) : $(OBJECTS)
	$(CC) -o $@ $< $(LIB)

bench.o : bench.cpp
	$(CC) -c $(DEPENDENCE) $< $(CCFLAG) $(INC)

.PHONY : clean
clean :
	-rm $(TARGETS) $(OBJECTS)
