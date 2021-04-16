#=========================== Makefile for Poseidon ============================#

## Select Number of Processes (Mostly for MacBook)

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

PETSC_MODE	=ON
#PETSC_MODE	=OFF


## Compile with HDF5

HDF5_MODE       =ON
#HDF5_MODE      =OFF


## Compile with MPI

MPI_MODE        =ON
#MPI_MODE       =OFF


Poseidon_Dir = /Users/nickroberts/Poseidon/Restructure

include ./Build/Makefile_Core

VPATH += ./Obj

#---------------------------- Compilation Rules ------------------------------------#


main  : $(POSEIDON)
	@echo "         compiling with $(COMP_$(MACHINE_NAME)) :"
	$(FORT) -c $(STD) $(OUTPUT_LINKER) $(OBJ) $(INCLUDE_LINKER) $(OBJ) $(DRV)/$(CODE_drv:%.o=%.$(EXT)) -o $(DRIVER_OBJ)
	$(FORT) $(STD) $(OBJ)/*.o -o $(BIN)/main.x
	@echo ">>> compiled on `hostname -s` with $(FORT_$(MACHINE_NAME)) <<<"



PoseidonLib: $(POSEIDON)
	ar crv $(OBJ)/Poseidon.a $(OBJ)/*.o





#-------------------------- Execution Rule ------------------------------------#


run:
	$(BIN)/main.x

run_yahil:
	mpirun  -np $(NPROCS) ./$(BIN)/main_$(DIMENSION).x

run_mpi:
	mpirun  -np $(NPROCS) ./$(BIN)/main_$(DIMENSION).x

run_replot:
	mpirun ./$(BIN)/replot.x


#------------------------------- Clean Up Rule  ------------------------------------#

clean:
	rm -f $(OBJ)/*.o
	rm -f $(OBJ)/*.mod
	rm -f $(BIN)/*.x


clean_output:
	rm -f $(OUT)/*.out
	rm -f $(OUT)/Results/*.out
	rm -f $(OUT)/Reports/*.out
	rm -f $(OUT)/Reports/Iteration_Reports/*.out
	rm -f $(OUT)/Objects/*.out
	rm -f $(OUT)/Objects/Sources/*.out
	rm -f $(OUT)/Objects/Linear_System/*.out
	rm -f $(OUT)/Objects/Residual/*.out
