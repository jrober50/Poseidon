   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE XCFC_System_Solvers_TypeB_Module                                      !##!
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
                    Num_R_Nodes,                &
                    LM_Length

USE Variables_FP,  &
            ONLY :  FP_Coeff_Vector_A,          &
                    FP_Coeff_Vector_B,          &
                    FP_Source_Vector_A,         &
                    FP_Source_Vector_B,         &
                    Beta_Diagonals,             &
                    Beta_MVL_Banded,            &
                    Beta_IPIV,                  &
                    Beta_Factorized_Flag,       &
                    MCF_Flag,                   &
                    Factored_NNZ,               &
                    Laplace_Factored_Val,       &
                    Laplace_Factored_Col,       &
                    Laplace_Factored_Row,       &
                    FP_Update_Vector,           &
                    CFA_Var_Map

USE Poseidon_Cholesky_Module,   &
            ONLY :  CCS_Back_Substitution,      &
                    CCS_Forward_Substitution,   &
                    Cholesky_Factorization


USE FP_Factorize_Beta_Banded, &
            ONLY :  Factorize_Beta_Banded,      &
                    Jacobi_PC_MVL_Banded_Vector

USE FP_Functions_BC,  &
            ONLY :  DIRICHLET_BC_Beta_Banded,   &
                    Dirichlet_BC_CHOL,          &
                    Neumann_BC_CCS


USE Variables_MPI, &
            ONLY :  myID_Poseidon,              &
                    MasterID_Poseidon,          &
                    nPROCS_Poseidon,            &
                    Poseidon_Comm_World

USE Poseidon_MPI_Utilities_Module, &
            ONLY :  STOP_MPI,                   &
                    MPI_Master_Print,           &
                    MPI_All_Print

USE MPI_Communication_TypeB_Module,             &
            ONLY :  MPI_RTM_Source_TypeB,       &
                    MPI_BCast_Coeffs_TypeB


IMPLICIT NONE


CONTAINS





!+301+###########################################################################!
!                                                                                !
!           Call Solve_FP_System                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Solve_System_TypeB(iU, iVB)


INTEGER, DIMENSION(3), INTENT(IN)                                   :: iU
INTEGER,               INTENT(IN)                                   :: iVB

INTEGER                                                             ::  INFO
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                      ::  WORK_VEC

INTEGER                                                             ::  Lower_Limit
INTEGER                                                             ::  Upper_Limit
INTEGER                                                             ::  ierr, i





IF ( Verbose_Flag ) THEN
    IF ( iVB == iVB_X ) THEN
        PRINT*,"--In XCFC Iteration, Begining X System Solve."
    ELSE IF ( iVB == iVB_S ) THEN
        PRINT*,"--In XCFC Iteration, Begining Shift System Solve."
    ELSE
        WRITE(*,'(A)') "Incompatable iVB value passed to XCFC_Solve_System_TypeB."
        WRITE(*,'(A,3I3.3)') "iVB value received ",iVB
        WRITE(*,'(A)') "iVB must be 1 or 2."

    END IF
END IF





#ifdef POSEIDON_AMREX_FLAG
Lower_Limit = 1
Upper_Limit = Beta_Prob_Dim
CALL MPI_RTM_Source_TypeB(  iVB,                    &
                            Lower_Limit,            &
                            Upper_Limit,            &
                            MasterID_Poseidon,      &
                            Poseidon_Comm_World,    &
                            ierr                    )

#endif


IF ( myID_Poseidon == MasterID_Poseidon ) THEN


    IF ( .NOT. Beta_Factorized_Flag ) THEN
        CALL Factorize_Beta_Banded()
    END IF



    ALLOCATE( WORK_VEC( 1:Beta_Prob_Dim ) )
    Work_Vec = FP_Source_Vector_B(:,iVB)



    CALL DIRICHLET_BC_Beta_Banded(Beta_Prob_Dim, Work_Vec )
    CALL Jacobi_PC_MVL_Banded_Vector( Work_Vec )

!    PRINT*,Work_Vec

    CALL ZGBTRS( 'N',                   &
                 Beta_Prob_Dim,         &
                 Beta_Diagonals,        &
                 Beta_Diagonals,        &
                 1,                     &
                 Beta_MVL_Banded,       &
                 3*Beta_Diagonals+1,    &
                 Beta_IPIV,             &
                 Work_Vec,             &
                 Beta_Prob_Dim,         &
                 INFO                   )

    IF (INFO .NE. 0) THEN
        print*,"ZGBTRS has failed with INFO = ",INFO
    END IF

!    PRINT*,"Coeff_Vec"
!    PRINT*,Work_Vec

    FP_Coeff_Vector_B(:,iVB) = Work_Vec(:)


    DEALLOCATE( Work_Vec )
END IF

!PRINT*,"Stopping in XCFC_System_Solvers_TypeB"
!STOP



#ifdef POSEIDON_AMREX_FLAG
Lower_Limit = 1
Upper_Limit = Num_R_Nodes
CALL MPI_BCAST_Coeffs_TypeB(iVB,                     &
                            Lower_Limit,            &
                            Upper_Limit,            &
                            MasterID_Poseidon,      &
                            Poseidon_Comm_World,    &
                            ierr                    )


#endif


END SUBROUTINE XCFC_Solve_System_TypeB








END MODULE XCFC_System_Solvers_TypeB_Module
