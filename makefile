#=========================== Makefile for Poseidon ============================#

# Select Dimension

#DIMENSION=1D
#DIMENSION=2D
DIMENSION=3D

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





#-------------------------------- Code Objects ---------------------------------#

#
#	Contains the list of object files that define the operations of the
#	Poseidon code.
#
include makefile_objects






#-------------------------------- Code Dictionary -----------------------------#

#
#	Contains the definitions of terms used below in the compilation of
#	the Poseidon code.
#
include makefile_dictionary





#---------------------------- Compilation Rules ------------------------------------#


yahil : $(CODE_com) $(CODE_par) $(CODE_o) $(CODE_ext) $(CODE_itf)
	@echo "         compiling with $(COMP_$(MACHINE_NAME)) :"
	$(FORT) -c $(STD) $(OUTPUT_LINKER) $(OBJ) $(INCLUDE_LINKER) $(OBJ) $(SRC)/$(CODE_drv:%.o=%.$(EXT)) -o $(DRIVER_OBJ)
	$(FORT) $(STD) $(OBJ)/*.o -o $(BIN)/main_$(DIMENSION).x
	@echo ">>> compiled on `hostname -s` with $(FORT_$(MACHINE_NAME)) <<<"


replot: $(CODE_com) $(CODE_par) $(CODE_o) $(CODE_ext) $(CODE_itf)
	@echo "         compiling with $(COMP_$(MACHINE_NAME)) :"
	$(FORT) -c $(STD) $(OUTPUT_LINKER) $(OBJ) $(INCLUDE_LINKER) $(OBJ) $(SRC)/d.replot.F90 -o $(OBJ)/d.replot.o
	$(FORT) $(STD) $(OBJ)/*.o -o $(BIN)/replot.x
	@echo ">>> compiled on `hostname -s` with $(FORT_$(MACHINE_NAME)) <<<"



PoseidonLib: $(CODE_com) $(CODE_par) $(CODE_o) $(CODE_ext) $(CODE_itf)
	ar crv $(OBJ)/Poseidon.a $(OBJ)/*.o



$(CODE_com):%.o: $(SRC)/%.$(EXT)
	$(FORT) -c $(STD) $(OUTPUT_LINKER) $(OBJ) $(INCLUDE_LINKER) $(OBJ) $< -o $(OBJ)/$@

$(CODE_par):%.o: $(SRC)/%.$(EXT)
	$(FORT) -c $(STD) $(OUTPUT_LINKER) $(OBJ) $(INCLUDE_LINKER) $(OBJ) $< -o $(OBJ)/$@

$(CODE_o):%.o:  $(SRC)/%.$(EXT)
	$(FORT) -c $(STD) $(OUTPUT_LINKER) $(OBJ) $(INCLUDE_LINKER) $(OBJ) $< -o $(OBJ)/$@

$(CODE_ext):%.o: $(SRC)/%.$(EXT)
	$(FORT) -c $(STD) $(OUTPUT_LINKER) $(OBJ) $(INCLUDE_LINKER) $(OBJ) $< -o $(OBJ)/$@

$(CODE_itf):%.o: $(SRC)/%.$(EXT)
	$(FORT) -c $(STD) $(OUTPUT_LINKER) $(OBJ) $(INCLUDE_LINKER) $(OBJ) $< -o $(OBJ)/$@




#-------------------------- Execution Rule ------------------------------------#


run:
	$(BIN)/main_$(DIMENSION).x

run_yahil:
	mpirun  -np $(NPROCS) ./$(BIN)/main_$(DIMENSION).x

run_mpi:
	mpirun  -np $(NPROCS) ./$(BIN)/main_$(DIMENSION).x

run_replot:
	mpirun ./$(BIN)/replot.x


#------------------------------- Clean Up Rule  ------------------------------------#

clean:
	rm -f $(OBJ)/com.*.o
	rm -f $(OBJ)/z.*.o
	rm -f $(OBJ)/d.*.o
	rm -f $(OBJ)/p.*.o
	rm -f $(OBJ)/i.*.o
	rm -f $(OBJ)/e.*.o
	rm -f $(OBJ)/*.mod DONE
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
