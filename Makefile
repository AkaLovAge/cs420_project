SHELL=/bin/bash

# Matrix dimensions
M=64
N=64

CC=mpicc

OPT_LEVEL=-O3

OPENMP_FLAG=-fopenmp

#basic setting
SRC_DIR=./src
INCLUDE_DIR=./include
BIN_DIR=./bin
TEST_DIR=./test
OBJ_DIR=./obj

#get the name of file
SOURCE=$(wildcard $(SRC_DIR)/*.c)
OBJS=$(SOURCE:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)


#openmp configuration 
CC_VERSION=$(shell $(CC) --version)
ifeq ($(findstring Intel,$(CC_VERSION)),Intel)
  OPENMP_FLAG=-qopenmp
endif

# Used for checking the results.
ERROR_THRESHOLD=1e-4

# The papi library location.

COMMON_PROG_ARGS=-DN=$(N) -DM=$(M) -std=c99 $(OPT_LEVEL) -DERROR_THRESHOLD=$(ERROR_THRESHOLD)

#Program arguments 
cc_PROG_ARGS=$(COMMON_PROG_ARGS)

#Include tage 
ACT_INCLUDE_TAG=-I$(INCLUDE_DIR)

all: $(BIN_DIR)/mpi_householder \
	$(BIN_DIR)/householder \
	$(BIN_DIR)/2d_omp

$(OBJ_DIR)/%.o:$(SRC_DIR)/%.c
	$(CC) -c $< -o $@ $(ACT_INCLUDE_TAG)

$(BIN_DIR)/householder: $(TEST_DIR)/householder.c $(OBJS)
	$(CC) $^ -o $@ $(cc_PROG_ARGS) $(ACT_INCLUDE_TAG)

$(BIN_DIR)/mpi_householder: $(TEST_DIR)/mpi_householder.c $(OBJS)
	$(CC) $^ -o $@ $(cc_PROG_ARGS) $(ACT_INCLUDE_TAG)

$(BIN_DIR)/2d_omp: $(TEST_DIR)/2d_omp.c $(OBJS)
	$(CC) $^ -o $@ $(cc_PROG_ARGS) $(ACT_INCLUDE_TAG) $(OPENMP_FLAG)

clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(BIN_DIR)/*
