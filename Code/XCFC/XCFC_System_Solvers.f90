   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE XCFC_System_Solvers_Module                                            !##!
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
        ONLY :  DEGREE,                     &
                L_LIMIT,                    &
                Verbose_Flag

USE Parameters_Variable_Indices, &
        ONLY :  iVB_X,                      &
                iVB_S,                      &
                iU_X1,                      &
                iU_X2,                      &
                iU_X3,                      &
                iU_S1,                      &
                iU_S2,                      &
                iU_S3

USE Variables_Derived, &
        ONLY :  Beta_Prob_Dim,              &
                Num_R_Nodes

USE Variables_FP,  &
        ONLY :  FP_Coeff_Vector_A,            &
                FP_Coeff_Vector_B,          &
                FP_Source_Vector_A,         &
                FP_Source_Vector_B,         &
                Beta_Diagonals,             &
                Beta_MVL_Banded,            &
                Beta_IPIV,                  &
                Beta_Factorized_Flag,       &
                MCF_Flag,              &
                Factored_NNZ,               &
                Laplace_Factored_Val,       &
                Laplace_Factored_Col,       &
                Laplace_Factored_Row,       &
                FP_Update_Vector,           &
                CFA_Var_Map

USE Poseidon_Cholesky_Module,   &
        ONLY :  CCS_Back_Substitution,          &
                CCS_Forward_Substitution,       &
                Cholesky_Factorization


USE FP_Factorize_Beta_Banded, &
        ONLY :  Factorize_Beta_Banded,          &
                Jacobi_PC_MVL_Banded_Vector

USE FP_Functions_BC,  &
        ONLY :  DIRICHLET_BC_Beta_Banded,       &
                Dirichlet_BC_CHOL,              &
                Neumann_BC_CCS

USE FP_Functions_Mapping, &
        ONLY :  FP_LM_Map

IMPLICIT NONE


CONTAINS



!+301+###########################################################################!
!                                                                                !
!           Call Solve_FP_System                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Solve_System_TypeA(iU)

INTEGER, INTENT(IN)                                             ::  iU

COMPLEX(KIND = idp), DIMENSION(NUM_R_NODES)                     ::  WORK_VEC
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                  ::  WORK_ELEM_VAL

INTEGER                                                         ::  l, m
INTEGER                                                         ::  lm_loc





IF ( Verbose_Flag ) THEN
    IF( iU == 1 ) THEN
        WRITE(*,'(A)')"--In XCFC Iteration, Begining Conformal Factor System Solve."
    ELSE IF ( iU == 2 ) THEN
        WRITE(*,'(A)')"--In XCFC Iteration, Begining Lapse Function System Solve."
    ELSE
        WRITE(*,'(A)')"Incompatable iU value passed to XCFC_Solve_System_TypeA."
        WRITE(*,'(A,3I3.3)')"iU value received ",iU
        WRITE(*,'(A)')"iU must be 1 or 2."
    END IF
END IF



IF ( MCF_Flag == 0 ) THEN
    
    !
    !   This only needs to be done everytime the radial mesh is defined/redefined.
    !   This performs Cholesky factorization on the stiffness matrix and overwrites
    !   the stiffness matrix variables STF_ELEM_VAL, STF_ROW_IND, and STF_COL_PTR to
    !   represent the factorization matrix, L.  This matrix can then be reused to
    !   solve the linear system using forward and backward substitution.
    !

    CALL Cholesky_Factorization()
    MCF_Flag = 1

END IF


ALLOCATE( Work_Elem_Val(0:Factored_NNZ-1))




DO l = 0,L_LIMIT
DO m = -l,l

    lm_loc = FP_LM_Map(l,m)
    WORK_VEC = -FP_Source_Vector_A(:,lm_loc,iU)
    WORK_ELEM_VAL(:) = Laplace_Factored_VAL(:,l,CFA_Var_Map(iU))

    
    CALL DIRICHLET_BC_CHOL( NUM_R_NODES,                &
                            Factored_NNZ,               &
                            l,                          &
                            m,                          &
                            Laplace_Factored_COL(:,l),  &
                            Laplace_Factored_ROW(:,l),  &
                            WORK_VEC,                   &
                            iU                          )



    CALL NEUMANN_BC_CCS(    NUM_R_NODES,                &
                            Factored_NNZ,               &
                            l,                          &
                            m,                          &
                            WORK_ELEM_VAL,              &
                            Laplace_Factored_COL(:,l),  &
                            Laplace_Factored_ROW(:,l),  &
                            WORK_VEC                    )




    CALL CCS_Forward_Substitution(  NUM_R_NODES,                    &
                                    Factored_NNZ,                   &
                                    WORK_ELEM_VAL,                  &
                                    Laplace_Factored_COL(:,l),      &
                                    Laplace_Factored_ROW(:,l),      &
                                    WORK_VEC                        )


    CALL CCS_Back_Substitution(     NUM_R_NODES,                    &
                                    Factored_NNZ,                   &
                                    WORK_ELEM_VAL,                  &
                                    Laplace_Factored_COL(:,l),      &
                                    Laplace_Factored_ROW(:,l),      &
                                    WORK_VEC                        )


    FP_Update_Vector(:,lm_loc,iU) = WORK_VEC(:)-FP_Coeff_Vector_A(:,lm_loc,iU)
    FP_Coeff_Vector_A( :,lm_loc,iU) = WORK_VEC(:)


END DO  ! m Loop
END DO  ! l Loop


DEALLOCATE( Work_Elem_Val )

END SUBROUTINE XCFC_Solve_System_TypeA








!+301+###########################################################################!
!                                                                                !
!           Call Solve_FP_System                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Solve_System_TypeB(iU)


INTEGER, DIMENSION(3), INTENT(IN)                                   :: iU

INTEGER                                                             ::  INFO
INTEGER                                                             ::  iVB
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                      ::  WORK_VEC

!REAL(idp), DIMENSION(1:4)                                           ::  timer





IF ( Verbose_Flag ) THEN
    IF (iU(1) == iU_X1) THEN
        PRINT*,"--In XCFC Iteration, Begining X System Solve."
    ELSE IF ( iU(1) == iU_S1 ) THEN
        PRINT*,"--In XCFC Iteration, Begining Shift System Solve."
    ELSE

    END IF
END IF

IF ( iU(1) == iU_X1 ) THEN
    iVB = iVB_X
ELSE IF ( iU(1) == iU_S1 ) THEN
    iVB = iVB_S
ELSE
    WRITE(*,'(A)') "Incompatable iU value passed to XCFC_Solve_System_TypeB."
    WRITE(*,'(A,3I3.3)') "iU values received ",iU
    WRITE(*,'(A)') "iU must be 3,4,5 or 6,7,8."
END IF



IF ( .NOT. Beta_Factorized_Flag ) THEN
    CALL Factorize_Beta_Banded()
END IF



ALLOCATE( WORK_VEC( 1:Beta_Prob_Dim ) )
Work_Vec = FP_Source_Vector_B(:,iVB)


!PRINT*,Work_Vec


CALL DIRICHLET_BC_Beta_Banded(Beta_Prob_Dim, Work_Vec )
CALL Jacobi_PC_MVL_Banded_Vector( Work_Vec )


CALL ZGBTRS( 'N',                   &
             Beta_Prob_Dim,         &
             Beta_Diagonals,        &
             Beta_Diagonals,        &
             1,                     &
             Beta_MVL_Banded,       &
             3*Beta_Diagonals+1,    &
             Beta_IPIV,             &
             -Work_Vec,             &
             Beta_Prob_Dim,         &
             INFO                   )

IF (INFO .NE. 0) THEN
    print*,"ZGBTRS has failed with INFO = ",INFO
END IF



FP_Coeff_Vector_B(:,iVB) = Work_Vec(:)


DEALLOCATE( Work_Vec )




END SUBROUTINE XCFC_Solve_System_TypeB






END MODULE XCFC_System_Solvers_Module
