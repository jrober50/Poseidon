   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Linear_System_Solvers_TypeA_Module                                    !##!
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



 !+101+############################################################!
!                                                                   !
!          Solve_Linear_System_TypeA                                !
!                                                                   !
 !#################################################################!
SUBROUTINE Solve_Linear_System_TypeA(iU)

INTEGER, INTENT(IN)                                             ::  iU

COMPLEX(KIND = idp), DIMENSION(NUM_R_NODES)                     ::  WORK_VEC
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                  ::  WORK_ELEM_VAL

INTEGER                                                         ::  l, m
INTEGER                                                         ::  lm_loc, ierr

INTEGER                                                         ::  Lower_Limit
INTEGER                                                         ::  Upper_Limit

CHARACTER(LEN = 300)                    ::  Message

IF ( Verbose_Flag ) THEN
    WRITE(Message,'(A,A,A)')'Beginning ',TRIM(CFA_Var_Names(iU)),' Linear Solve.'
    CALL Run_Message(TRIM(Message))
END IF



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
        WORK_VEC = -cVA_Load_Vector(:,lm_loc,iU)
        WORK_ELEM_VAL(:) = Laplace_Factored_VAL(:,l)

!        PRINT*,WORK_ELEM_VAL(:)


!        PRINT*,"Work_Vec, Type A l = ",l," m = ",m
!        PRINT*,Work_Vec
!        PRINT*,"+++++++++++++++++++++++++++++++++++"

!        STOP
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


!        PRINT*,"After"
!        PRINT*,Work_Vec(:)

        FP_Update_Vector(:,lm_loc,iU) = WORK_VEC(:)-cVA_Coeff_Vector(:,lm_loc,iU)
        cVA_Coeff_Vector( :,lm_loc,iU) = WORK_VEC(:)



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
!        PRINT*,cVA_Coeff_Vector(:,lm_loc, iU)
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




!PRINT*,cVA_Coeff_Vector( :,1,iU)




END SUBROUTINE Solve_Linear_System_TypeA






END MODULE Linear_System_Solvers_TypeA_Module