HDF5_DIR       = /usr/local/hdf5

AMREX_HOME    ?= /Users/nickroberts/amrex
#AMREX_HOME    ?= $(AMREX_DIR)
HDF5_HOME     ?= $(HDF5_DIR)
LAPACK_HOME   ?= $(LAPACK_DIR)


HDF5_INC       = $(HDF5_HOME)/include
HDF5_LIB       = $(HDF5_HOME)/lib

#LAPACK_INC     = $(LAPACK_HOME)/include
#LAPACK_LIB     = $(LAPACK_HOME)/lib

USE_F_INTERFACES = TRUE

DEBUG     = FALSE
#DEBUG     = TRUE

USE_MPI   = TRUE
USE_OMP   = FALSE

#MEM_PROFILE  = FALSE
#TINY_PROFILE = FALSE

DIM       = 1

EBASE     = Yahil_Profile_XCFC



#INCLUDE_LOCATIONS += $(LAPACK_INC)
#LIBRARY_LOCATIONS += $(LAPACK_LIB)
LIBRARIES += -framework Accelerate


# Set pre-processor directives
# This is also where, e.g., USE_HDF5 and USE_PETSC are used
include $(AMREX_HOME)/Tools/GNUMake/Make.defs

# Include basic AMReX framework
include $(AMREX_HOME)/Src/Base/Make.package

# Include directories and source files needed for amrex_parallel_module
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/F_Interfaces/Base/Make.package
include $(AMREX_HOME)/Src/F_Interfaces/AmrCore/Make.package



# Include directories and source files needed from Poseidon
POSEIDON_MACHINE_NAME ?= jrober50
include /Users/nickroberts/Poseidon/Build/Make.Poseidon




# Include dicrectories and source files needed for the Driver
include Makefile_Driver
include Makefile_Path
include AMReX_Routines/Makefile_AMReX_Routines

f90EXE_sources += $(DRIVER_f90)

F90EXE_sources += $(DRIVER_AMREX_ROUTINES)
F90EXE_sources += $(DRIVER_F90)




# Include rules for compilation
include $(AMREX_HOME)/Tools/GNUMake/Make.rules





.PHONY : clean_output
clean_output :
	rm -f $(OUT)/*.out
	rm -f $(OUT)/Results/*.out
	rm -f $(OUT)/Reports/*.out
	rm -f $(OUT)/Reports/Iteration_Reports/*.out
	rm -f $(OUT)/Objects/*.out
	rm -f $(OUT)/Objects/Sources/*.out
	rm -f $(OUT)/Objects/Linear_System/*.out
	rm -f $(OUT)/Objects/Residual/*.out
	rm -f $(OUT)/Objects/Mesh/*.out
	rm -f $(OUT)/Objects/Coefficients/*.out


