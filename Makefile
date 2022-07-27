#=========================== Makefile for Poseidon ============================#

## Select Number of Processes (Mostly for MacBook)

NAME := KnotAFish

NPROCS=1

## Select Machine

MACHINE_NAME    =NicksMacBook
#MACHINE_NAME = sjdunham
#MACHINE_NAME	=MacBook
#MACHINE_NAME	=Rhea
#MACHINE_NAME	=BlueWaters


## Select Mode

CMODE	=DEBUG
#CMODE	=OPTIMIZE


## Compile with Openmp

#OPENMP_MODE	=ON
OPENMP_MODE	=OFF



## Compile with PETSc

#PETSC_MODE	=ON
PETSC_MODE	=OFF


## Compile with HDF5

HDF5_MODE       =ON
#HDF5_MODE      =OFF


## Compile with MPI

MPI_MODE        =ON
#MPI_MODE       =OFF

#AMREX_MODE      =ON
AMREX_MODE      =OFF



include ./Build/Makefile_Core

VPATH += ./Obj

#---------------------------- Compilation Rules ------------------------------------#

#objb = $(POSEIDON).o
#$(foreach Poseidon, $(POSEIDON), $(eval $(POSEIDON) := $(objb)))

#OBJS := $(addsuffix .o, $(POSEIDON))
#DEPS := $(addsuffix .d, $(POSEIDON))
#LIB  := $(patsubst %, lib.%a, $(NAME))


#$(info $$OBJS is [${OBJS}])

$(LIB): $(OBJS)
	ar crv $@ $^


.PHONY: PoseidonLib clean clean_output

PoseidonLib: $(POSEIDON_o)
	ar crv $(OBJ)/poseidon.a $(OBJ)/*.o




#------------------------------- Clean Up Rule  ------------------------------------#

clean:
	rm -f $(OBJ)/*.o
	rm -f $(OBJ)/*.mod
	rm -f $(BIN)/*.x
	rm -f $(OBJ)/*.a

clean_output:
	rm -f $(OUT)/*.out
	rm -f $(OUT)/Results/*.out
	rm -f $(OUT)/Reports/*.out
	rm -f $(OUT)/Reports/Iteration_Reports/*.out
	rm -f $(OUT)/Objects/*.out
	rm -f $(OUT)/Objects/Sources/*.out
	rm -f $(OUT)/Objects/Linear_System/*.out
	rm -f $(OUT)/Objects/Residual/*.out
