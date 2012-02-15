# on Franklin and Hopper, we will benchmark you against Cray LibSci, the default vendor-tuned BLAS. The Cray compiler wrappers handle all the linking. If you wish to compare with other BLAS implementations, check the NERSC documentation.
# This makefile is intended for the GNU C compiler. On Franklin and Hopper, the Portland compilers are default, so you must instruct the Cray compiler wrappers to switch to GNU: type "module swap PrgEnv-pgi PrgEnv-gnu"
# Your code must compile (with GCC) with the given CFLAGS. You may experiment with the OPT variable to invoke additional compiler options.

CC = cc 
#OPT = -O3 -ffast-math -funroll-loops -fstrict-aliasing -funsafe-loop-optimizations -msse -msse2 -msse3 -mfpmath=both -march=amdfam10 -fprofile-generate
#OPT = -O3 -ffast-math -funroll-loops -fstrict-aliasing -funsafe-loop-optimizations -msse -msse2 -msse3 -mfpmath=both -march=amdfam10 -fprofile-use
OPT = -O3 -ffast-math -funroll-loops -fstrict-aliasing -funsafe-loop-optimizations -msse -msse2 -msse3 -mfpmath=both -march=amdfam10 -ftree-vectorize 
#OPT = -O0 -g -ggdb
#OPT = -O3 -ffast-math -funroll-loops -fstrict-aliasing -funsafe-loop-optimizations -msse -msse2 -msse3 -mfpmath=both -march=amdfam10 -ftree-vectorize -S -fno-inline


CFLAGS = -Wall -std=gnu99 $(OPT) #-fprofile-use
LDFLAGS = -Wall #-g -ggdb
#LDFLAGS = -Wall -fprofile-use
#LDFLAGS = -Wall -fprofile-generate
# librt is needed for clock_gettime
LDLIBS = -lrt

sources = $(wildcard dgemm-*.c)
targets = $(sources:dgemm-%.c=benchmark-%)
objects = benchmark.o $(sources:.c=.o)

.PHONY : default
default : all

.PHONY : all
all : clean_ $(targets)
#all : $(targets)

.PHONY :  clean_
clean_:
	rm -f dgemm-11.o

benchmark-naive : benchmark.o dgemm-naive.o 
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)
benchmark-recursive : benchmark.o dgemm-recursive.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)
benchmark-11 : benchmark.o dgemm-11.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)
benchmark-blocked : benchmark.o dgemm-blocked.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)
benchmark-blas : benchmark.o dgemm-blas.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

%.o : %.c
	$(CC) -c $(CFLAGS) $<

.PHONY : clean
clean:
	rm -f $(targets) $(objects)
