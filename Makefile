# on Franklin and Hopper, we will benchmark you against Cray LibSci, the default vendor-tuned BLAS. The Cray compiler wrappers handle all the linking. If you wish to compare with other BLAS implementations, check the NERSC documentation.
# This makefile is intended for the GNU C compiler. On Franklin and Hopper, the Portland compilers are default, so you must instruct the Cray compiler wrappers to switch to GNU: type "module swap PrgEnv-pgi PrgEnv-gnu"
# Your code must compile (with GCC) with the given CFLAGS. You may experiment with the OPT variable to invoke additional compiler options.

CC = cc 
#OPT = -O3 -ffast-math -funroll-loops -fstrict-aliasing -funsafe-loop-optimizations -march=native -msse -msse2 -msse3 -mfpmath=both
#OPT = -O3 -ffast-math -funroll-loops -fstrict-aliasing -funsafe-loop-optimizations -march=native -msse3 -mfpmath=both
#OPT = -O3 -ffast-math -funroll-loops -fstrict-aliasing -funsafe-loop-optimizations -msse -msse2 -msse3 -mfpmath=sse -march=amdfam10 -fprofile-arcs
#OPT = -O3 -ffast-math -funroll-loops -fstrict-aliasing -funsafe-loop-optimizations -msse -msse2 -msse3 -mfpmath=both -march=amdfam10 -fprofile-use
OPT = -O3 -ffast-math -funroll-loops -fstrict-aliasing -funsafe-loop-optimizations -msse -msse2 -msse3 -mfpmath=both -march=amdfam10 -ftree-vectorize
#OPT = -g -ggdb
#OPT = -O3 -ffast-math -funroll-loops -fno-inline -fstrict-aliasing -funsafe-loop-optimizations -march=native -msse3 -mfpmath=both -S
#PGI
#OPT = -fast -fastsse -O4 -Mvect -Munroll
#Pathscale
#OPT = -Ofast -std=gnu99 -msse -msse2 -msse3 -msse4a -ffast-math
#OPT = -std=gnu99 -msse -msse2 -msse3 -msse4a -O4 -ffast-math


CFLAGS = -Wall -std=gnu99 $(OPT) #-fprofile-use
#PGI
#CFLAGS = -Minform=inform -c99 $(OPT) #-fprofile-use
#CFLAGS = $(OPT) -std=gnu99
LDFLAGS = -Wall #-g -ggdb
#LDFLAGS = -Wall -fprofile-use
#LDFLAGS = -Wall -fprofile-arcs
#PGI
#LDFLAGS = -Minform=inform
#LDFLAGS = 
# librt is needed for clock_gettime
LDLIBS = -lrt

sources = $(wildcard dgemm-*.c)
targets = $(sources:dgemm-%.c=benchmark-%)
objects = benchmark.o $(sources:.c=.o)

.PHONY : default
default : all

.PHONY : all
#all : clean $(targets)
all : clean_ $(targets)
#all : $(targets)

.PHONY :  clean_
clean_:
	rm -f dgemm-jornada-auto.o

benchmark-naive : benchmark.o dgemm-naive.o 
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)
benchmark-recursive : benchmark.o dgemm-recursive.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)
benchmark-register-blocked : benchmark.o dgemm-register-blocked.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)
benchmark-jornada-sse : benchmark.o dgemm-jornada-sse.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)
benchmark-jornada-auto : benchmark.o dgemm-jornada-auto.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)
benchmark-jornada : benchmark.o dgemm-jornada.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)
benchmark-blocked : benchmark.o dgemm-blocked.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)
benchmark-blas : benchmark.o dgemm-blas.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)
benchmark-sse-t : benchmark.o dgemm-see-t.c
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)
benchmark-sse-no-t : benchmark.o dgemm-sse-no-t.c
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

%.o : %.c
	$(CC) -c $(CFLAGS) $<

.PHONY : clean
clean:
	rm -f $(targets) $(objects)
