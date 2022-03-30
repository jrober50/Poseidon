   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poisson_Linear_Solve_Module                                           !##!
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
            ONLY :  Degree,             &
                    L_Limit,            &
                    Verbose_Flag

USE Variables_Mesh, &
            ONLY :  Num_R_Elements

USE Variables_Derived, &
            ONLY :  LM_LENGTH,                  &
                    Num_R_Nodes

USE Variables_Poisson, &
            ONLY :  Source_Vector,      &
                    Coefficient_Vector, &
                    STF_NNZ,            &
                    STF_Mat_Integrals,  &
                    STF_Elem_Val,       &
                    STF_Row_Ind,        &
                    STF_Col_Ptr,        &
                    Matrix_Cholesky_Factorized_Flag

USE Poisson_Cholesky_Module, &
            ONLY :  Cholesky_Factorization

USE Poisson_Boundary_Conditions_Module, &
            ONLY :  DIRICHLET_BC_CHOL,  &
                    NEUMANN_BC_CCS

USE Poisson_Forward_Substitution_Module, &
            ONLY :  CCS_Forward_Substitution

USE Poisson_Back_Substitution_Module, &
            ONLY :  CCS_Back_Substitution


IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!                                                     				!
!                                                                               !
!###############################################################################!
SUBROUTINE Poisson_Linear_Solve()


REAL(idp)                                               ::  SCALE_FACTOR
COMPLEX(idp), DIMENSION(0:NUM_R_NODES-1)                ::  WORK_VEC

REAL(idp), ALLOCATABLE, DIMENSION(:)                    ::  WORK_ELEM_VAL

INTEGER                                                 ::  NNZ
INTEGER                                                 ::  l, m, lm,k


IF ( Verbose_Flag ) THEN
    PRINT*,"- Begining Linear Solve."
END IF

NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1
ALLOCATE(WORK_ELEM_VAL(0:NNZ-1))






IF ( Matrix_Cholesky_Factorized_Flag .EQV. .FALSE. ) THEN


    !
    !   This only needs to be done everytime the radial mesh is defined/redefined.
    !   This performs Cholesky factorization on the stiffness matrix and overwrites
    !   the stiffness matrix variables STF_ELEM_VAL, STF_ROW_IND, and STF_COL_PTR to
    !   represent the factorization matrix, L.  This matrix can then be reused to
    !   solve the linear system using forward and backward substitution.
    !

    CALL Cholesky_Factorization()
    Matrix_Cholesky_Factorized_Flag = .TRUE.


END IF











DO l = 0,L_LIMIT
DO m = -l,l
    !#######################################################################!
    !                                                                       !
    !               CCS Cholesky Factorization Matrix Solver                !
    !                                                                       !
    !#######################################################################!

    lm = l*(l+1)+m+1

    WORK_VEC = Source_Vector(:,lm)


    CALL DIRICHLET_BC_CHOL( NUM_R_NODES,    &
                            STF_NNZ,        &
                            l, m,           &
                            STF_COL_PTR,    &
                            STF_ROW_IND,    &
                            WORK_VEC        )


    CALL NEUMANN_BC_CCS(    NUM_R_NODES,    &
                            STF_NNZ,        &
                            l, m,           &
                            WORK_ELEM_VAL,  &
                            STF_COL_PTR,    &
                            STF_ROW_IND,    &
                            WORK_VEC        )



    CALL CCS_Forward_Substitution(  NUM_R_NODES,        &
                                    STF_NNZ,            &
                                    STF_ELEM_VAL(:,l),  &
                                    STF_COL_PTR,        &
                                    STF_ROW_IND,        &
                                    WORK_VEC            )


    CALL CCS_Back_Substitution( NUM_R_NODES,        &
                                STF_NNZ,            &
                                STF_ELEM_VAL(:,l),  &
                                STF_COL_PTR,        &
                                STF_ROW_IND,        &
                                WORK_VEC            )




    Do k = 0,NUM_R_NODES - 1
        Coefficient_Vector(k,m,l) = WORK_VEC(k)
    END DO



END DO
END DO





DEALLOCATE(WORK_ELEM_VAL)


END SUBROUTINE Poisson_Linear_Solve









END MODULE Poisson_Linear_Solve_Module
