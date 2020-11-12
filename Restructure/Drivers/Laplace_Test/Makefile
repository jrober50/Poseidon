include ../../Build/Makefile_Core

TEST_NAME    =Laplace_Test


## Select Number of Processes (Mostly for MacBook)

NPROCS=1

## Select Machine

MACHINE_NAME    =NicksMacBook
#MACHINE_NAME = sjdunham
#MACHINE_NAME   =MacBook
#MACHINE_NAME   =Rhea
#MACHINE_NAME   =BlueWaters


## Select Mode

CMODE   =DEBUG
#CMODE  =OPTIMIZE


## Compile with Openmp

#OPENMP_MODE    =ON
OPENMP_MODE     =OFF



## Compile with PETSc

PETSC_MODE      =ON
#PETSC_MODE     =OFF


## Compile with HDF5

HDF5_MODE       =ON
#HDF5_MODE      =OFF


## Compile with MPI

MPI_MODE        =ON
#MPI_MODE       =OFF


Main : $(POSEIDON)
	@echo "         compiling with $(COMP_$(MACHINE_NAME)) :"
	$(FORT) -c $(STD) $(OUTPUT_LINKER) $(OBJ) $(INCLUDE_LINKER) $(OBJ) Laplace_Test.f90  -o $(OBJ)/Laplace_Test.o
	$(FORT) $(STD) $(OBJ)/*.o -o $(BIN)/$(TEST_NAME).x
	@echo ">>> compiled on `hostname -s` with $(FORT_$(MACHINE_NAME)) <<<"

Mapping : $(POSEIDON)
	@echo "         compiling with $(COMP_$(MACHINE_NAME)) :"
	$(FORT) -c $(STD) $(OUTPUT_LINKER) $(OBJ) $(INCLUDE_LINKER) $(OBJ) Laplace_Mapping.f90  -o $(OBJ)/Laplace_Mapping.o
	$(FORT) $(STD) $(OBJ)/*.o -o $(BIN)/Mapping.x
	@echo ">>> compiled on `hostname -s` with $(FORT_$(MACHINE_NAME)) <<<"




run : 
	./$(BIN)/$(TEST_NAME).x

run_Mapping : 
	./$(BIN)/Laplace_Mapping.x

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



