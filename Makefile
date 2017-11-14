SHELL=/bin/bash

# Matrix dimensions
M=1000
N=1000

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
#PAPI_LIB_DIR=/usr/local/apps/papi/5.4.1/lib
#PAPI_INC_DIR=/usr/local/apps/papi/5.4.1/include
COMMON_PROG_ARGS=-DN=$(N) -DM=$(M) -std=c99 $(OPT_LEVEL) -DERROR_THRESHOLD=$(ERROR_THRESHOLD)

#Program arguments 
cc_PROG_ARGS=$(COMMON_PROG_ARGS)

#Include tage 
ACT_INCLUDE_TAG=-I$(INCLUDE_DIR)

all: $(BIN_DIR)/2D_householder \

$(OBJ_DIR)/%.o:$(SRC_DIR)/%.c
	$(CC) -c $< -o $@ $(ACT_INCLUDE_TAG)

$(BIN_DIR)/2D_householder: $(TEST_DIR)/householder.c $(OBJS)
	$(CC) $^ -o $@ $(cc_PROG_ARGS) $(ACT_INCLUDE_TAG)

clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(BIN_DIR)/*
