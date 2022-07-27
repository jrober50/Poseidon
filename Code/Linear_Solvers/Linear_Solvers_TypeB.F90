   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Linear_System_Solvers_TypeB_Module                                    !##!
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
            ONLY :  Run_Message,            &
                    Warning_Message

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Parameters_Variable_Indices, &
            ONLY :  iVB_X,                      &
                    iVB_S,                      &
                    iU_X1,                      &
                    iU_X2,                      &
                    iU_X3,                      &
                    iU_S1,                      &
                    iU_S2,                      &
                    iU_S3

USE Poseidon_IO_Parameters, &
            ONLY :  CFA_VecVar_Names

USE Variables_Derived, &
            ONLY :  Beta_Prob_Dim,              &
                    Num_R_Nodes,                &
                    LM_Length

USE Variables_Vectors,  &
            ONLY :  cVB_Coeff_Vector,          &
                    cVB_Load_Vector

USE Variables_Matrices,  &
            ONLY :  Beta_Diagonals,             &
                    Beta_MVL_Banded,            &
                    Beta_IPIV


USE Matrix_Vector_Laplacian_Routines, &
            ONLY :  Factorize_Vector_Laplacian,     &
                    Jacobi_PC_MVL_Banded_Vector

USE Matrix_Boundary_Condition_Routines,  &
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

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Matrices_Flags,    &
                    iPF_Init_Matrices_Type_B_LU

IMPLICIT NONE


CONTAINS





 !+101+############################################################!
!                                                                   !
!          Solve_Linear_System_TypeB                                !
!                                                                   !
 !#################################################################!
SUBROUTINE Solve_Linear_System_TypeB(iU, iVB)


INTEGER, DIMENSION(3), INTENT(IN)                                   :: iU
INTEGER,               INTENT(IN)                                   :: iVB

INTEGER                                                             ::  INFO
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                      ::  WORK_VEC

INTEGER                                                             ::  Lower_Limit
INTEGER                                                             ::  Upper_Limit
INTEGER                                                             ::  ierr

CHARACTER(LEN = 300)                    ::  Message


IF ( Verbose_Flag ) THEN
    WRITE(Message,'(A,A,A)')'Beginning ',TRIM(CFA_VecVar_Names(iVB)),' Linear Solve.'
    CALL Run_Message(TRIM(Message))
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


    IF ( .NOT. lPF_Init_Matrices_Flags(iPF_Init_Matrices_Type_B_LU) ) THEN
        CALL Factorize_Vector_Laplacian()
    END IF



    ALLOCATE( WORK_VEC( 1:Beta_Prob_Dim ) )
    Work_Vec = cVB_Load_Vector(:,iVB)



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
        WRITE(Message,'(A,I1.1)')"ZGBTRS has failed with INFO = ",INFO
        CALL Warning_Message(TRIM(Message))
    END IF

!    PRINT*,"Coeff_Vec"
!    PRINT*,Work_Vec

    cVB_Coeff_Vector(:,iVB) = Work_Vec(:)


    DEALLOCATE( Work_Vec )
END IF

!PRINT*,"Stopping in XCFC_System_Solvers_TypeB"
!STOP



#ifdef POSEIDON_AMREX_FLAG
Lower_Limit = 1
Upper_Limit = Beta_Prob_Dim
CALL MPI_BCAST_Coeffs_TypeB(iVB,                     &
                            Lower_Limit,            &
                            Upper_Limit,            &
                            MasterID_Poseidon,      &
                            Poseidon_Comm_World,    &
                            ierr                    )


#endif


END SUBROUTINE Solve_Linear_System_TypeB








END MODULE Linear_System_Solvers_TypeB_Module
