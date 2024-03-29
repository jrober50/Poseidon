   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Initialization_FEM_Module                                       	     !##!
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

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations

USE Allocation_FEM, &
            ONLY :  Allocate_FEM

IMPLICIT NONE


CONTAINS



 !+101+####################################################!
!                                                           !
!          Initialization_FEM                         	    !
!                                                           !
 !#########################################################!
SUBROUTINE Initialization_FEM

CALL ALlocate_FEM()
FEM_Node_xlocs = Initialize_LGL_Quadrature_Locations(DEGREE)

END SUBROUTINE Initialization_FEM


END MODULE Initialization_FEM_Module
