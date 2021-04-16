   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Initialization_Quadrature                                                    !##!
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

USE Poseidon_Numbers_Module, &
                ONLY :  pi

USE Poseidon_Parameters, &
                ONLY :  Verbose_Flag

USE Variables_Quadrature, &
                ONLY :  Num_R_Quad_Points,      &
                        Num_T_Quad_Points,      &
                        Num_P_Quad_Points,      &
                        Num_TP_Quad_Points,     &
                        Int_R_Locations,        &
                        Int_T_Locations,        &
                        Int_P_Locations,        &
                        Int_R_Weights,          &
                        Int_T_Weights,          &
                        Int_P_Weights,          &
                        Int_TP_Weights

USE Functions_Quadrature, &
                ONLY :  Initialize_LG_Quadrature,       &
                        Initialize_Trapezoid_Quadrature

USE Allocation_Quadrature, &
                ONLY :  Allocate_Quadrature

IMPLICIT NONE

CONTAINS

 !+101+################################################################################!
!                                                                                       !
!       Initialize_Quad                                                                 !
!                                                                                       !
 !#####################################################################################!
SUBROUTINE Initialize_Quadrature()

INTEGER                                 :: td, pd


IF ( Verbose_Flag ) THEN
    PRINT*,"-Initializing Quadrature variables. "
END IF

CALL Allocate_Quadrature


CALL Initialize_LG_Quadrature(NUM_R_QUAD_POINTS, INT_R_LOCATIONS, INT_R_WEIGHTS)
CALL Initialize_LG_Quadrature(NUM_T_QUAD_POINTS, INT_T_LOCATIONS, INT_T_WEIGHTS)
CALL Initialize_LG_Quadrature(NUM_P_QUAD_POINTS, INT_P_LOCATIONS, INT_P_WEIGHTS)
!CALL Initialize_Trapezoid_Quadrature(NUM_P_QUAD_POINTS, INT_P_LOCATIONS, INT_P_WEIGHTS)


IF ( NUM_T_QUAD_POINTS == 1 ) THEN
    ! 1D Cheat !
    INT_T_WEIGHTS = 4.0_idp/pi
END IF


DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS
        INT_TP_WEIGHTS( (td-1)*NUM_P_QUAD_POINTS+pd ) = INT_T_WEIGHTS(td)*INT_P_WEIGHTS(pd)
    END DO
END DO


END SUBROUTINE Initialize_Quadrature





END MODULE Initialization_Quadrature
