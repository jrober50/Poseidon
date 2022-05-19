   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_XCFC_Interface_Module                                               !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the top level subroutines needed to run the Poseidon in XCFC mode. !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+                                                                       !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

USE Poseidon_Kinds_Module, &
            ONLY : idp

USE Poseidon_Parameters, &
            ONLY :  DEGREE

USE Parameters_Variable_Indices, &
            ONLY :  iVB_X,                      &
                    iVB_S,                      &
                    iU_CF,                      &
                    iU_LF,                      &
                    iU_S1,                      &
                    iU_S2,                      &
                    iU_S3,                      &
                    iU_X1,                      &
                    iU_X2,                      &
                    iU_X3


USE Variables_Tables, &
            ONLY :  Ylm_Values,                 &
                    Ylm_dt_Values,              &
                    Ylm_dp_Values,              &
                    Lagrange_Poly_Table


USE Variables_Derived, &
            ONLY :  LM_LENGTH,          &
                    Beta_Prob_Dim

USE Variables_FP, &
            ONLY :  FP_Coeff_Vector_A,  &
                    FP_Coeff_Vector_B

USE Variables_Mesh, &
            ONLY :  rlocs,              &
                    tlocs,              &
                    plocs

USE Variables_Mesh, &
            ONLY :  rlocs,              &
                    tlocs,              &
                    plocs

USE Maps_Fixed_Point, &
            ONLY :  FP_Array_Map_TypeB

USE Maps_Quadrature, &
            ONLY :  Map_To_tpd

USE Maps_Domain, &
            ONLY :  Map_To_FEM_Node

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations

USE Functions_Math, &
            ONLY :  Lagrange_Poly


USE XCFC_Method_Module, &
            ONLY :  XCFC_Method_Part1,      &
                    XCFC_Method_Part2


CONTAINS


!+101+######################################################################################!
!                                                                                           !
!       Poseidon_XCFC_Run_Part1 - Performs the steps to calculate the XCFC Conformal        !
!                                 Factor.  Upon completion of this call, the routine        !
!                                 Poseidon_Return_ConFactor() can be used to acquire        !
!                                 values of the conformal factor at desired locations.      !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_XCFC_Run_Part1()


CALL XCFC_Method_Part1()


END SUBROUTINE Poseidon_XCFC_Run_Part1





!+102+######################################################################################!
!                                                                                           !
!       Poseidon_XCFC_Run_Part2 - Performs the steps to calculate the XCFC Lapse Function,  !
!                                 and Shift Vector.  Upon completion of this call, the      !
!                                 routine, Poseidon_Return_Lapse() and                      !
!                                 Poseidon_Return_Shift() can be used to acquire            !
!                                 values of the those variables at desired locations.       !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_XCFC_Run_Part2()


CALL XCFC_Method_Part2()


END SUBROUTINE Poseidon_XCFC_Run_Part2






















END MODULE Poseidon_XCFC_Interface_Module
