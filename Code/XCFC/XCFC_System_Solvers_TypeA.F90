   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE XCFC_System_Solvers_TypeA_Module                                      !##!
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
            ONLY :  iU_CF,                      &
                    iU_LF

USE Variables_Derived, &
            ONLY :  Beta_Prob_Dim,              &
                    Num_R_Nodes,                &
                    LM_Length

USE Variables_MPI, &
            ONLY :  myID_Poseidon,              &
                    MasterID_Poseidon,          &
                    Poseidon_Comm_World,        &
                    nPROCS_Poseidon

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

USE Functions_Domain_Maps, &
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

IMPLICIT NONE


CONTAINS



!+301+###########################################################################!
!                                                                                !
!          XCFC_Solve_System_TypeA                                               !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Solve_System_TypeA(iU)

INTEGER, INTENT(IN)                                             ::  iU

COMPLEX(KIND = idp), DIMENSION(NUM_R_NODES)                     ::  WORK_VEC
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                  ::  WORK_ELEM_VAL

INTEGER                                                         ::  l, m, i
INTEGER                                                         ::  lm_loc, ierr

INTEGER                                                         ::  Lower_Limit
INTEGER                                                         ::  Upper_Limit



IF ( Verbose_Flag ) THEN
    IF( iU == iU_CF ) THEN
        WRITE(*,'(A)')"--In XCFC Iteration, Begining Conformal Factor System Solve.", myID_Poseidon
    ELSE IF ( iU == iU_LF ) THEN
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




IF ( iU == iU_CF ) THEN
    CALL TimerStart( Timer_XCFC_ConFactor_LinearSolve)
ELSEIF ( iU == iU_LF ) THEN
    CALL TimerStart( Timer_XCFC_Lapse_LinearSolve)
END IF



! Reduce Souce Vector To Master !
#ifdef POSEIDON_AMREX_FLAG
Lower_Limit = 1
Upper_Limit = Num_R_Nodes
CALL MPI_RTM_Source_TypeA(  iU,                     &
                            Lower_Limit,            &
                            Upper_Limit,            &
                            MasterID_Poseidon,      &
                            Poseidon_Comm_World,    &
                            ierr                    )

#endif







IF ( myID_Poseidon == MasterID_Poseidon ) THEN

    ALLOCATE( Work_Elem_Val(0:Factored_NNZ-1))

    DO l = 0,L_LIMIT
    DO m = -l,l

        lm_loc = Map_To_lm(l,m)
        WORK_VEC = -FP_Source_Vector_A(:,lm_loc,iU)
        WORK_ELEM_VAL(:) = Laplace_Factored_VAL(:,l,CFA_Var_Map(iU))

!        PRINT*,"Work_Vec, Type A l = ",l," m = ",m
!        PRINT*,Work_Vec
!        PRINT*,"+++++++++++++++++++++++++++++++++++"

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

!        PRINT*,FP_Coeff_Vector_A( :,lm_loc,iU)


    END DO  ! m Loop
    END DO  ! l Loop


    DEALLOCATE( Work_Elem_Val )
END IF







#ifdef POSEIDON_AMREX_FLAG

Lower_Limit = 1
Upper_Limit = Num_R_Nodes
CALL MPI_BCAST_Coeffs_TypeA(iU,                     &
                            Lower_Limit,            &
                            Upper_Limit,            &
                            MasterID_Poseidon,      &
                            Poseidon_Comm_World,    &
                            ierr                    )

!DO i = 0,nPROCS_Poseidon-1
!IF ( myID_Poseidon == i ) THEN
!    PRINT*,"myID_Poseidon ",i
!    DO lm_loc = 1,LM_Length
!        PRINT*,FP_Coeff_Vector_A(:,lm_loc, iU)
!    END DO
!END IF
!CALL MPI_Barrier(Poseidon_Comm_World,ierr)
!END DO

#endif




IF ( iU == iU_CF ) THEN
    CALL TimerStop( Timer_XCFC_ConFactor_LinearSolve)
ELSEIF ( iU == iU_LF ) THEN
    CALL TimerStop( Timer_XCFC_Lapse_LinearSolve)
END IF




!PRINT*,FP_Coeff_Vector_A( :,1,iU)




END SUBROUTINE XCFC_Solve_System_TypeA






END MODULE XCFC_System_Solvers_TypeA_Module
