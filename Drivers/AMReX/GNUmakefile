HDF5_DIR       = /usr/local/hdf5

AMREX_HOME    ?= /Users/nickroberts/amrex
HDF5_HOME     ?= $(HDF5_DIR)
LAPACK_HOME   ?= $(LAPACK_DIR)


HDF5_INC       = $(HDF5_HOME)/include
HDF5_LIB       = $(HDF5_HOME)/lib

#LAPACK_INC     = $(LAPACK_HOME)/include
#LAPACK_LIB     = $(LAPACK_HOME)/lib

USE_F_INTERFACES = TRUE

DEBUG     = TRUE

USE_MPI   = TRUE
USE_OMP   = FALSE

#MEM_PROFILE  = FALSE
#TINY_PROFILE = FALSE

DIM = 3

EBASE     = AMReX_Mapping

PETSC_MODE = OFF
AMREX_MODE = ON


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
include /Users/nickroberts/Poseidon/Build/Make.Poseidon
include Makefile_Driver
f90EXE_sources += $(DRIVER_f90)
F90EXE_sources += $(DRIVER_F90)


DEFINES += $(POSEIDON_DEFS)

# Include rules for compilation
include $(AMREX_HOME)/Tools/GNUMake/Make.rules






