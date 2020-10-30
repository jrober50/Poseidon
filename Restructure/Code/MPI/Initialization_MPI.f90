   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Initialization_MPI                                                           !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!


!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!

USE Poseidon_Kinds_Module, &
                ONLY : idp

USE Variables_Mesh, &
                ONLY :  Num_R_Elements,             &
                        Num_T_Elements,             &
                        Num_P_Elements
USE Variables_MPI, &
                ONLY :  Num_R_Elems_Per_Shell,      &
                        Num_T_Elems_Per_Block,      &
                        Num_P_Elems_Per_Block,      &
                        Ratio_BNDLperBLCK,          &
                        Ratio_T_BNDLperBLCK,        &
                        Ratio_P_BNDLperBLCK

USE Functions_MPI, &
                ONLY : Create_Poseidon_Communicators

IMPLICIT NONE

CONTAINS

 !+101+################################################################################!
!                                                                                       !
!       Initialize_Quad                                                                 !
!                                                                                       !
 !#####################################################################################!
SUBROUTINE Initialize_MPI()


Num_R_Elems_Per_Shell = Num_R_Elements
Num_T_Elems_Per_Block = Num_T_Elements
Num_P_Elems_Per_Block = Num_P_Elements

Ratio_T_BNDLperBLCK = 1
Ratio_P_BNDLperBLCK = 1
Ratio_BNDLperBLCK = Ratio_T_BNDLperBLCK * Ratio_P_BNDLperBLCK

CALL Create_Poseidon_Communicators()


END SUBROUTINE Initialize_MPI





END MODULE Initialization_MPI

