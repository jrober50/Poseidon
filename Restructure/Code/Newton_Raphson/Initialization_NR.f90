   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Initialization_NR                                                            !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the top level subroutines needed to inialize, run, and close       !##!
!##!        Poseidon.                                                               !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!    +101+   Poseidon_Initialize                                                 !##!
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
USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag


USE Variables_Functions, &
            ONLY :  Calc_3D_Values_At_Location,    &
                    Calc_1D_CFA_Values

USE Allocation_NR, &
            ONLY :  Allocate_NR

USE Functions_NR,   &
            ONLY :  Calc_NR_Values_At_Location,  &
                    Calc_1D_CFA_Values_NR

USE mpi




IMPLICIT NONE




                    !*F&S*==========================================!
                    !                                               !
                    !           Functions & Subroutines             !
                    !                                               !
                    !===============================================!
CONTAINS










 !+101+####################################################################################!
!                                                                                           !
!       Initialize_FP                                                          !
!                                                                                           !
!===========================================================================================!
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Initialize_NR( )



IF ( Verbose_Flag ) THEN
    PRINT*,"-Initializing Newton-Raphson Method variables. "
END IF

CALL Allocate_NR()

Calc_3D_Values_At_Location  => Calc_NR_Values_At_Location
Calc_1D_CFA_Values          => Calc_1D_CFA_Values_NR



END SUBROUTINE Initialize_NR










END MODULE Initialization_NR


