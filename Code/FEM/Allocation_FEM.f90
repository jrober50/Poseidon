   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Allocation_FEM                                                  	     !##!
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
USE Variables_FEM_Module, &
            ONLY :  FEM_Node_xlocs

USE Poseidon_Parameters, &
            ONLY :  Degree

IMPLICIT NONE


CONTAINS



 !+101+####################################################!
!                                                           !
!          Allocate_FEM                                	    !
!                                                           !
 !#########################################################!
SUBROUTINE Allocate_FEM()

ALLOCATE( FEM_Node_xlocs(0:Degree) )

END SUBROUTINE Allocate_FEM



 !+101+####################################################!
!                                                           !
!          Deallocate_FEM                                   !
!                                                           !
 !#########################################################!
SUBROUTINE Deallocate_FEM()

DEALLOCATE( FEM_Node_xlocs )

END SUBROUTINE Deallocate_FEM




END MODULE Allocation_FEM
