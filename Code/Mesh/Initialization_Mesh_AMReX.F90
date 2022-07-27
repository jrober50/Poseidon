   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Initialization_Mesh_AMReX_Module                                      !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!                                                                         !##!
!###############################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !##########################################################################!


!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE Poseidon_Kinds_Module, &
            ONLY : idp

USE Poseidon_Message_Routines_Module, &
            ONLY :  Init_Message

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Num_Levels,   &
                    iNumLeafElements

USE Variables_Tables, &
            ONLY :  Level_DX

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    Num_T_Elements,         &
                    Num_P_Elements,         &
                    rlocs,                  &
                    tlocs,                  &
                    plocs,                  &
                    drlocs,                 &
                    dtlocs,                 &
                    dplocs,                 &
                    R_Inner

USE Variables_AMReX_Core, &
            ONLY :  Table_Offsets,           &
                    FEM_Elem_Table

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Mesh_Flags,    &
                    iPF_Init_Mesh_Init

IMPLICIT NONE


CONTAINS



 !+101+############################################################!
!                                                                   !
!          Determine_AMReX_Mesh                                     !
!                                                                   !
 !#################################################################!
SUBROUTINE Determine_AMReX_Mesh()

INTEGER                                         ::  lvl
INTEGER                                         ::  Here
INTEGER                                         ::  There
INTEGER                                         ::  re
INTEGER                                         ::  elem

IF ( Verbose_Flag ) CALL Init_Message('Initializing Mesh Variables.')

DO lvl = 0,AMReX_Num_Levels-1
    Here  = Table_Offsets(lvl)
    There = Table_Offsets(lvl+1)-1
    DO elem = Here,There
        drlocs(FEM_Elem_Table(elem)) = level_dx(lvl,1)
    END DO
END DO


rlocs(0) = R_Inner
DO re = 1,Num_R_Elements
    rlocs(re) = rlocs(re-1)+drlocs(re-1)
!    rlocs(re) = R_Inner + sum(drlocs(0:re-1)
END DO


lPF_Init_Mesh_Flags(iPF_Init_Mesh_Init) = .TRUE.

END SUBROUTINE Determine_AMReX_Mesh






END MODULE Initialization_Mesh_AMReX_Module
