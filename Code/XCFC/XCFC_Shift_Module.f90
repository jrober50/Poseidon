   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE XCFC_Shift_Module                                                         !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!


USE Poseidon_Kinds_Module, &
        ONLY :  idp


USE Poseidon_Parameters, &
        ONLY :  Verbose_Flag,               &
                Degree

USE Variables_Mesh, &
        ONLY :  Num_R_Elements

USE Variables_Derived, &
        ONLY :  Beta_Prob_Dim,              &
                LM_Length


USE Variables_FP,  &
        ONLY :  FP_Coeff_Vector,            &
                FP_Coeff_Vector_Beta,       &
                FP_Source_Vector_Beta,      &
                Beta_Diagonals,             &
                Beta_MVL_Banded,            &
                Beta_IPIV,                  &
                Beta_Factorized_Flag


USE FP_Beta_Banded, &
        ONLY :  Factorize_Beta_Banded,      &
                Jacobi_PC_MVL_Banded_Vector

USE FP_Functions_BC,  &
        ONLY :  DIRICHLET_BC_Beta_Banded

USE FP_Functions_Mapping,   &
        ONLY :  FP_Beta_Array_Map,          &
                FP_FEM_Node_Map
                
USE XCFC_Source_Vector_Module, &
        ONLY :  XCFC_Calc_Shift_Source

IMPLICIT NONE

CONTAINS



!+101+##########################################################################!
!                                                                               !
!                       XCFC_Shift_Solve                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Shift_Solve()

INTEGER                                                 ::  lm_loc
INTEGER                                                 ::  here, there
INTEGER                                                 ::  ui, re, d



IF ( Verbose_Flag ) THEN
    PRINT*,"Begining Shift system. "
END IF



CALL XCFC_Calc_Shift_Source()
Call XCFC_Solve_Shift_System()


DO ui = 3,5
DO re = 0,Num_R_Elements -1
DO d = 0,Degree
DO lm_loc = 1,LM_Length
    here = FP_Beta_Array_Map(re,d,ui-2,lm_loc)
    There = FP_FEM_Node_Map(re,d)
    FP_Coeff_Vector(There,lm_loc,ui) = FP_Coeff_Vector_Beta(Here)
END DO ! lm_loc
END DO ! d
END DO ! re
END DO ! ui


END SUBROUTINE XCFC_Shift_Solve
















!+201+###########################################################################!
!                                                                                !
!           XCFC_Solve_Shift_System                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Solve_Shift_System()


INTEGER                                                                     ::  INFO
INTEGER, DIMENSION(1:Beta_Prob_Dim)                                         ::  IPIV


COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                              ::  WORK_VEC
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:)                            ::  WORK_MAT


INTEGER                                                                     ::  i, j, Col, Row
INTEGER                                                                     ::  ui, re, d, l
INTEGER                                                                     ::  uj, rep, dp, lp


INTEGER                                                                     ::  Here, There

REAL(idp), DIMENSION(1:4)                                                   ::  timer

REAL(idp)                                                                   ::  RCOND

IF ( Verbose_Flag ) THEN
    PRINT*,"--In XCFC Iteration, In XCFC_Solve_Shift_System."
END IF

IF ( .NOT. Beta_Factorized_Flag ) THEN
    CALL Factorize_Beta_Banded()
END IF



ALLOCATE( WORK_VEC( 1:Beta_Prob_Dim ) )

Work_Vec = FP_Source_Vector_Beta


CALL DIRICHLET_BC_Beta_Banded(Beta_Prob_Dim, Work_Vec )

CALL Jacobi_PC_MVL_Banded_Vector( Work_Vec )



CALL ZGBTRS( 'N',                   &
             Beta_Prob_Dim,         &
             Beta_Diagonals,        &
             Beta_Diagonals,        &
             1,                     &
             Beta_MVL_Banded,              &
             3*Beta_Diagonals+1,    &
             Beta_IPIV,             &
             -Work_Vec,              &
             Beta_Prob_Dim,         &
             INFO                   )

IF (INFO .NE. 0) THEN
    print*,"ZGBTRS has failed with INFO = ",INFO
END IF

FP_Coeff_Vector_Beta(:) = Work_Vec(:)



DEALLOCATE( Work_Vec )



END SUBROUTINE XCFC_Solve_Shift_System











END MODULE XCFC_Shift_Module

