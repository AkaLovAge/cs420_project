SHELL=/bin/bash

#
M=1000
N=1000

CC=mpicc

OPT_LEVEL=-o3

all: householder_2D_parallel:

householder_2D_parallel: householder_2D_parallel.c
	$(CC) src/householder_2D_parallel.c -o h2d_$(N)_$(M)

clean:
	rm -r bin
