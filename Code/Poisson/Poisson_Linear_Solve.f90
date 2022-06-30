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

USE Poseidon_Message_Routines_Module, &
            ONLY :  Run_Message

USE Poseidon_Parameters, &
            ONLY :  DEGREE,                     &
                    L_LIMIT,                    &
                    Verbose_Flag

USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,                      &
                    iU_LF
USE Poseidon_IO_Parameters, &
            ONLY :  CFA_Var_Names

USE Variables_Derived, &
            ONLY :  Beta_Prob_Dim,              &
                    Num_R_Nodes,                &
                    LM_Length

USE Variables_MPI, &
            ONLY :  myID_Poseidon,              &
                    MasterID_Poseidon,          &
                    Poseidon_Comm_World,        &
                    nPROCS_Poseidon

USE Variables_Vectors,  &
            ONLY :  cVA_Coeff_Vector,           &
                    cVA_Load_Vector
                    
USE Variables_Matrices,  &
            ONLY :  Factored_NNZ,               &
                    Laplace_Factored_Val,       &
                    Laplace_Factored_Col,       &
                    Laplace_Factored_Row

USE Variables_FP,  &
            ONLY :  FP_Update_Vector

USE Matrix_Cholesky_Factorization_Module,   &
            ONLY :  CCS_Back_Substitution,          &
                    CCS_Forward_Substitution,       &
                    Cholesky_Factorization

USE Matrix_Boundary_Condition_Routines,  &
            ONLY :  DIRICHLET_BC_Beta_Banded,       &
                    Dirichlet_BC_CHOL,              &
                    Neumann_BC_CCS

USE Maps_Domain, &
            ONLY :  Map_To_lm

USE MPI_Communication_TypeA_Module,             &
            ONLY :  MPI_RTM_Source_TypeA,           &
                    MPI_BCast_Coeffs_TypeA

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_XCFC_Lapse_LinearSolve,   &
                    Timer_XCFC_ConFactor_LinearSolve

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Matrices_Flags,    &
                    iPF_Init_Matrices_Type_A_Cholesky

IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!                                                     				!
!                                                                               !
!###############################################################################!
SUBROUTINE Poisson_Linear_Solve()


COMPLEX(idp), DIMENSION(0:NUM_R_NODES-1)                ::  WORK_VEC

REAL(idp), ALLOCATABLE, DIMENSION(:)                    ::  WORK_ELEM_VAL

INTEGER                                                 ::  l, m, lm,k

INTEGER                                                 ::  iU

IF ( Verbose_Flag ) CALL Run_Message("Beginning Poisson Linear Solve.")







IF ( .NOT. lPF_Init_Matrices_Flags(iPF_Init_Matrices_Type_A_Cholesky) ) THEN


    !
    !   This only needs to be done everytime the radial mesh is defined/redefined.
    !   This performs Cholesky factorization on the stiffness matrix and overwrites
    !   the stiffness matrix variables STF_ELEM_VAL, STF_ROW_IND, and STF_COL_PTR to
    !   represent the factorization matrix, L.  This matrix can then be reused to
    !   solve the linear system using forward and backward substitution.
    !

    CALL Cholesky_Factorization()


END IF









IF ( myID_Poseidon == MasterID_Poseidon ) THEN

ALLOCATE( Work_Elem_Val(0:Factored_NNZ-1))

DO l = 0,L_LIMIT
DO m = -l,l
    !#######################################################################!
    !                                                                       !
    !               CCS Cholesky Factorization Matrix Solver                !
    !                                                                       !
    !#######################################################################!

    lm = Map_To_lm(l,m)
    WORK_VEC = -cVA_Load_Vector(:,lm_loc,iU)
    WORK_ELEM_VAL(:) = Laplace_Factored_VAL(:,l)

    CALL DIRICHLET_BC_CHOL( NUM_R_NODES,                &
                            Factored_NNZ,               &
                            l,                          &
                            m,                          &
                            Laplace_Factored_COL(:,l),  &
                            Laplace_Factored_ROW(:,l),  &
                            WORK_VEC,                   &
                            iU                           )


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




    cVA_Coeff_Vector( :,lm_loc,iU) = WORK_VEC(:)



END DO
END DO

END IF ! myID_Poseidon == MasterID_Poseidon



DEALLOCATE(WORK_ELEM_VAL)


END SUBROUTINE Poisson_Linear_Solve









END MODULE Poisson_Linear_Solve_Module
