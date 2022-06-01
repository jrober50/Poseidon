   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE IG_Flat_Guess_Module                                                  !##!
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

USE Poseidon_Parameters, &
            ONLY :  Method_Flag

USE IG_Input_Native_Module, &
            ONLY :  IG_Input_Native

USE Variables_Mesh, &
        ONLY :  NUM_R_ELEMENTS,             &
                NUM_T_ELEMENTS,             &
                NUM_P_ELEMENTS

USE Variables_Quadrature, &
        ONLY :  NUM_R_QUAD_POINTS,          &
                NUM_T_QUAD_POINTS,          &
                NUM_P_QUAD_POINTS,         &
                INT_R_LOCATIONS,            &
                INT_T_LOCATIONS,            &
                INT_P_LOCATIONS,            &
                Local_Quad_DOF,             &
                xLeftLimit,                 &
                xRightLimit

IMPLICIT NONE


CONTAINS



!+101+###########################################################################!
!                                                                                !
!               Poseidon_Input_Guess                                             !
!                                                                                !
!################################################################################!
SUBROUTINE IG_Init_Flat_Guess()


REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Psi_Guess
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  AlphaPsi_Guess
REAL(idp), DIMENSION(:,:,:,:,:), ALLOCATABLE            ::  Beta_Guess


REAL(idp)                                               :: Left_Limit
REAL(idp)                                               :: Right_Limit
INTEGER                                                 :: Num_DOF



Num_DOF =NUM_R_QUAD_POINTS*NUM_T_QUAD_POINTS*NUM_P_QUAD_POINTS

ALLOCATE( Psi_Guess(    1:Local_Quad_DOF,       &
                        0:Num_R_Elements-1,     &
                        0:Num_T_Elements-1,     &
                        0:Num_P_Elements-1  )   )
ALLOCATE( AlphaPsi_Guess(   1:Local_Quad_DOF,       &
                            0:Num_R_Elements-1,     &
                            0:Num_T_Elements-1,     &
                            0:Num_P_Elements-1  )   )
ALLOCATE( Beta_Guess(   1:Local_Quad_DOF,       &
                        0:Num_R_Elements-1,     &
                        0:Num_T_Elements-1,     &
                        0:Num_P_Elements-1,     &
                        1:3                 )   )




Psi_Guess = 1.0_idp
AlphaPsi_Guess = 1.0_idp
Beta_Guess = 0.0_idp


CALL IG_Input_Native( Psi_Guess,                                                &
                      AlphaPsi_Guess,                                           &
                      Beta_Guess,                                               &
                      Num_R_Elements, Num_R_Elements, Num_R_Elements,           &
                      NUM_R_QUAD_POINTS, NUM_T_QUAD_POINTS, NUM_P_QUAD_POINTS,  &
                      INT_R_LOCATIONS, INT_T_LOCATIONS, INT_P_LOCATIONS,        &
                      xLeftLimit, xRightLimit                     )





END SUBROUTINE IG_Init_Flat_Guess



END MODULE IG_Flat_Guess_Module
