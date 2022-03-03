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


END SUBROUTINE Determine_AMReX_Mesh






END MODULE Initialization_Mesh_AMReX_Module
