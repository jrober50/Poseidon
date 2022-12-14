   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Allocation_Maps                                                       !##!
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
            ONLY :  idp

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Num_Levels,       &
                    iNumLeafElements,       &
                    iLeafElementsPerLvl,    &
                    Findloc_Table,          &
                    FEM_Elem_Table,         &
                    Table_Offsets


IMPLICIT NONE


CONTAINS

 !+101+############################################################!
!                                                                   !
!          Allocate_AMReX_Maps                                      !
!                                                                   !
 !#################################################################!
SUBROUTINE Allocate_AMReX_Maps()

#ifdef POSEIDON_AMREX_FLAG
ALLOCATE( FindLoc_Table(0:iNumLeafElements-1)  )
ALLOCATE( FEM_Elem_Table(0:iNumLeafElements-1) )
ALLOCATE( Table_Offsets(0:AMReX_Num_Levels)    )
#endif

END SUBROUTINE Allocate_AMReX_Maps






 !+101+############################################################!
!                                                                   !
!          Deallocate_AMReX_Maps                                    !
!                                                                   !
 !#################################################################!
SUBROUTINE Deallocate_AMReX_Maps()

#ifdef POSEIDON_AMREX_FLAG
DEALLOCATE( FindLoc_Table  )
DEALLOCATE( FEM_Elem_Table )
DEALLOCATE( Table_Offsets  )
#endif

END SUBROUTINE Deallocate_AMReX_Maps




 !+101+############################################################!
!                                                                   !
!          Reallocate_AMReX_Maps                                      !
!                                                                   !
 !#################################################################!
SUBROUTINE Reallocate_AMReX_Maps()

#ifdef POSEIDON_AMREX_FLAG
CALL Deallocate_AMReX_Maps()
CALL Allocate_AMReX_Maps
#endif

END SUBROUTINE Reallocate_AMReX_Maps






END MODULE Allocation_Maps

