HDF5_DIR       = /usr/local/hdf5

AMREX_HOME    ?= $(AMREX_CG_DIR)
HDF5_HOME     ?= $(myHDF5_DIR)
LAPACK_HOME   ?= $(myLAPACK_DIR)



HDF5_INC       = $(HDF5_HOME)/include
HDF5_LIB       = $(HDF5_HOME)/lib

LAPACK_INC     = $(LAPACK_HOME)/include
LAPACK_LIB     = $(LAPACK_HOME)/lib

USE_F_INTERFACES = TRUE

#DEBUG     = FALSE
DEBUG     = TRUE

USE_MPI   = TRUE
USE_OMP   = FALSE

#MEM_PROFILE  = FALSE
#TINY_PROFILE = FALSE

DIM       = 1

EBASE     = Yahil_XCFC_AMReX



INCLUDE_LOCATIONS += $(LAPACK_INC)
LIBRARY_LOCATIONS += $(LAPACK_LIB)
#LIBRARIES += -framework Accelerate
LIBRARIES += -llapack -lblas

ifneq ("$(wildcard $(AMREX_HOME)/Src/DG/Make.package)","")
  # Requires amrex fork from github.com:dunhamsj/amrex
  # Branch: MeshRefinement_DG
  include $(AMREX_HOME)/Src/DG/Make.package
endif


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
POSEIDON_ACTIVE_DIR ?= $(POSEIDON_DIR)
include $(POSEIDON_ACTIVE_DIR)/Build/Make.Poseidon

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
	rm -f Poseidon_Output/*.out
	rm -f Poseidon_Output/Results/*.out
	rm -f Poseidon_Output/Reports/*.out
	rm -f Poseidon_Output/Reports/Iteration_Reports/*.out
	rm -f Poseidon_Output/Objects/*.out
	rm -f Poseidon_Output/Objects/Sources/*.out
	rm -f Poseidon_Output/Objects/Linear_System/*.out
	rm -f Poseidon_Output/Objects/Residual/*.out
	rm -f Poseidon_Output/Objects/Mesh/*.out
	rm -f Poseidon_Output/Objects/Coefficients/*.out


