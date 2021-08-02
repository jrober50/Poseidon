   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE XCFC_X_Module                                                                !##!
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
                FP_Coeff_Vector_X,          &
                FP_Source_Vector_X,         &
                Beta_Diagonals,             &
                Beta_MVL_Banded,            &
                Beta_IPIV,                  &
                Beta_Factorized_Flag


USE FP_Beta_Banded, &
        ONLY :  Factorize_Beta_Banded,          &
                Jacobi_PC_MVL_Banded_Vector

USE FP_Functions_BC,  &
        ONLY :  DIRICHLET_BC_Beta_Banded

USE FP_Functions_Mapping, &
        ONLY :  FP_Beta_Array_Map,          &
                FP_FEM_Node_Map

USE XCFC_Source_Vector_Module, &
        ONLY :  XCFC_Calc_X_Source

IMPLICIT NONE

CONTAINS

!+201+##########################################################################!
!                                                                               !
!                       Solve_X_System                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_X_Solve()


INTEGER                                                 ::  lm_loc
INTEGER                                                 ::  here, there
INTEGER                                                 ::  ui, re, d



IF ( Verbose_Flag ) THEN
    PRINT*,"Begining X system. "
END IF



CALL XCFC_Calc_X_Source()

Call XCFC_Solve_X_System()


DO ui = 6,8
DO re = 0,Num_R_Elements -1
DO d = 0,Degree
DO lm_loc = 1,LM_Length
    here = FP_Beta_Array_Map(re,d,ui-5,lm_loc)
    There = FP_FEM_Node_Map(re,d)
    FP_Coeff_Vector(There,lm_loc,ui) = FP_Coeff_Vector_X(Here)
END DO  ! l
END DO ! d
END DO ! re
END DO ! ui



END SUBROUTINE XCFC_X_Solve




!+301+###########################################################################!
!                                                                                !
!           Call Solve_FP_System                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Solve_X_System()


INTEGER                                                             ::  INFO

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                      ::  WORK_VEC

!REAL(idp), DIMENSION(1:4)                                           ::  timer





IF ( Verbose_Flag ) THEN
    PRINT*,"--In XCFC Iteration, In XCFC_Solve_X_System."
END IF






IF ( .NOT. Beta_Factorized_Flag ) THEN
    CALL Factorize_Beta_Banded()
END IF





ALLOCATE( WORK_VEC( 1:Beta_Prob_Dim ) )
Work_Vec = FP_Source_Vector_X


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



FP_Coeff_Vector_X(:) = Work_Vec(:)


DEALLOCATE( Work_Vec )




END SUBROUTINE XCFC_Solve_X_System











END MODULE XCFC_X_Module
