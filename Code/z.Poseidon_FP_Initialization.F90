   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_FP_Initialization_Module                                                   !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the top level subroutines needed to inialize, run, and close       !##!
!##!        Poseidon.                                                               !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!    +101+   Poseidon_Initialize                                                 !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!


!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE Poseidon_Constants_Module, &
            ONLY :  idp, pi, fdp


USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM,                 &
                    DEGREE,                     &
                    L_LIMIT,                    &
                    Num_CFA_Eqs,                &
                    Solver_Name,                &
                    Solver_Type_Flag,           &
                    Num_R_Quad_Points,          &
                    Num_T_Quad_Points,          &
                    Num_P_Quad_Points,          &
                    Num_Block_Theta_Rows,       &
                    Num_Block_Phi_Columns,      &
                    Ratio_bndlperblck,          &
                    Ratio_T_bndlperblck,        &
                    Ratio_P_bndlperblck,        &
                    R_Coarsen_Factor,           &
                    T_Coarsen_Factor,           &
                    P_Coarsen_Factor,           &
                    Data_Dist_Mode



USE Poseidon_Variables_Module, &
            ONLY :  R_INNER, R_OUTER,           &
                    NUM_R_ELEMENTS,             &
                    NUM_T_ELEMENTS,             &
                    NUM_P_ELEMENTS,             &
                    NUM_LOC_R_ELEMENTS,         &
                    NUM_LOC_T_ELEMENTS,         &
                    NUM_LOC_P_ELEMENTS,         &
                    NUM_TP_QUAD_POINTS,         &
                    NUM_R_NODES,                &
                    NUM_R_NODESp1,              &
                    rlocs,                      &
                    tlocs,                      &
                    plocs,                      &
                    RADIAL_MESH_SET_FLAG,       &
                    THETA_MESH_SET_FLAG,        &
                    PHI_MESH_SET_FLAG,          &
                    INT_R_LOCATIONS,            &
                    INT_R_WEIGHTS,              &
                    INT_T_LOCATIONS,            &
                    INT_T_WEIGHTS,              &
                    INT_P_LOCATIONS,            &
                    INT_P_WEIGHTS,              &
                    INT_TP_WEIGHTS,             &
                    LM_LENGTH,                  &
                    M_VALUES,                   &
                    ierr,                       &
                    drlocs,                     &
                    dtlocs,                     &
                    dplocs,                     &
                    LM_Location


USE Poseidon_FP_Variables_Module, &
            ONLY :  CFA_EQ_Flags,               &
                    CFA_EQ_Map,                 &
                    Laplace_NNZ,                &
                    Num_Matrices

USE Poseidon_Quadrature_Module, &
            ONLY :  Initialize_LG_Quadrature


USE Poseidon_Mesh_Module, &
            ONLY :  Generate_Defined_Coarse_Mesh


USE Poseidon_Tables_Module, &
            ONLY :  Initialize_Ylm_Tables,              &
                    Initialize_Lagrange_Poly_Tables

USE Poseidon_Allocation_Module, &
            ONLY :  Allocate_Poseidon_CFA_Variables,    &
                    Deallocate_Poseidon_CFA_Variables

USE Poseidon_FP_Allocation_Module, &
            ONLY :  Allocate_Poseidon_FP_Variables,     &
                    Deallocate_Poseidon_FP_Variables

USE Poseidon_MPI_Module, &
            ONLY :  CREATE_POSEIDON_COMMUNICATORS

USE Poseidon_Parameter_Read_Module, &
            ONLY :  UNPACK_POSEIDON_PARAMETERS

USE Poseidon_FP_Mapping_Functions_Module, &
            ONLY :  FP_LM_Map

USE Poseidon_Initialization_Module, &
            ONLY : OUTPUT_SETUP_TABLE

USE Poseidon_IO_Parameters, &
            ONLY : Solver_Names

USE Poseidon_FP_Laplace_Matrix_Module, &
            ONLY : Initialize_Laplace_Matrices

USE mpi




IMPLICIT NONE




                    !*F&S*==========================================!
                    !                                               !
                    !           Functions & Subroutines             !
                    !                                               !
                    !===============================================!
CONTAINS










 !+101+####################################################################################!
!                                                                                           !
!       Poseidon_FP_Init_From_File                                                          !
!                                                                                           !
!===========================================================================================!
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Poseidon_FP_Init_From_File( mode,                                                                    &
                                        CFA_EQ_Flags_Input,                                                     &
                                        Inner_Radius, Outer_Radius,                                             &
                                        R_Elements_Input, T_Elements_Input, P_Elements_Input,                   &
                                        Local_R_Elements_Input, Local_T_Elements_Input, Local_P_Elements_Input, &
                                        Input_Delta_R_Vector, Input_Delta_T_Vector, Input_Delta_P_Vector        )



                                         !                          !
                                        !!      Input Variables     !!
                                         !

INTEGER, INTENT(IN)                                                             ::  mode,                   &
                                                                                    R_Elements_Input,       &
                                                                                    T_Elements_Input,       &
                                                                                    P_Elements_Input,       &
                                                                                    Local_R_Elements_Input, &
                                                                                    Local_T_Elements_Input, &
                                                                                    Local_P_Elements_Input

INTEGER, DIMENSION(1:5), INTENT(IN)                                             ::  CFA_EQ_Flags_Input

REAL(KIND = idp), INTENT(IN)                                                    ::  Inner_Radius,           &
                                                                                    Outer_Radius


REAL(KIND = idp), DIMENSION(1:R_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_R_Vector
REAL(KIND = idp), DIMENSION(1:T_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_T_Vector
REAL(KIND = idp), DIMENSION(1:P_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_P_Vector





                                         !                              !
                                        !!     Subroutine Variables     !!
                                         !                              !

REAL(KIND = idp), DIMENSION(1:R_Elements_Input)                                 ::  Delta_R_Vector
REAL(KIND = idp), DIMENSION(1:T_Elements_Input)                                 ::  Delta_T_Vector
REAL(KIND = idp), DIMENSION(1:P_Elements_Input)                                 ::  Delta_P_Vector


INTEGER                                                                         ::  td, pd
INTEGER                                                                         ::  nPROCS, tmpID
INTEGER                                                                         ::  l, i, j



REAL(KIND = idp)                                                                :: timea, timeb, timec



CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nPROCS,ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, tmpID, ierr)


 !                                      !
!!  Load Poseidon Parameters From File  !!
 !                                      !
CALL UNPACK_POSEIDON_PARAMETERS()


Solver_Name = Solver_Names(2)


 !                                          !
!!  Set Global Variables to Input Values    !!
 !                                          !
R_INNER = Inner_Radius
R_OUTER = Outer_Radius


NUM_R_ELEMENTS = R_Elements_Input/R_Coarsen_Factor
NUM_T_ELEMENTS = T_Elements_Input/T_Coarsen_Factor
NUM_P_ELEMENTS = P_Elements_Input/P_Coarsen_Factor

NUM_LOC_R_ELEMENTS = Local_R_Elements_Input
NUM_LOC_T_ELEMENTS = Local_T_Elements_Input
NUM_LOC_P_ELEMENTS = Local_P_Elements_Input

!
!   *_NUM_R_NODES - Number of total radial nodes in total, per block, per subshell
!
NUM_R_NODES             = DEGREE*NUM_R_ELEMENTS + 1
NUM_R_NODESp1           = NUM_R_NODES + 1


!
!   NUM_TP_QUAD_POINTS = number of angular quadrature points per radial point
!
NUM_TP_QUAD_POINTS = NUM_T_QUAD_POINTS*NUM_P_QUAD_POINTS



!
!   CFA_EQ_Flags = Turns on/off which equations are solved.
!
CFA_EQ_Flags = CFA_EQ_Flags_Input
NUM_CFA_Eqs = SUM(CFA_EQ_Flags)
Solver_Type_Flag = 2

CFA_EQ_Map = -1
j = 1
DO i = 1,5
    IF ( CFA_EQ_Flags(i) == 1 ) THEN
        CFA_EQ_Map(j) = i
        j = j+1
    END IF
END DO


!
!   Calculate the number of matrices to be created and stored. 
!
! The Psi and AlphaPsi equations can share the same Lapace Matrix.
! Each of the Shift Vector componenets will need its own matrix.
Num_Matrices = 0
IF ( (CFA_EQ_Flags(1) == 1) .OR. (CFA_EQ_Flags(2) == 1) ) THEN
    Num_Matrices = 1
END IF
DO i = 3,5
    IF ( CFA_EQ_Flags(i) == 1 ) THEN
        Num_Matrices = Num_Matrices + 1
    END IF
END DO




!
!   Associate the Correct Map Functions, and Set Spherical Harmonic Length
!
IF ( DOMAIN_DIM == 1 ) THEN

    LM_LENGTH = 1

ELSE IF ( DOMAIN_DIM == 2 ) THEN

    LM_LENGTH = L_LIMIT + 1

ELSE IF ( DOMAIN_DIM == 3 ) THEN

    LM_LENGTH = (L_LIMIT + 1)*(L_LIMIT + 1)

END IF





Laplace_NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1




                                 !                                      !
                                !!      Allocate Space for Poseidon     !!
                                 !                                      !
CALL Allocate_Poseidon_CFA_Variables()
CALL Allocate_Poseidon_FP_Variables()



                                 !                                          !
                                !!      Initialize Reusable Quadratures     !!
                                 !                                          !

CALL Initialize_LG_Quadrature(NUM_R_QUAD_POINTS, INT_R_LOCATIONS, INT_R_WEIGHTS)
CALL Initialize_LG_Quadrature(NUM_T_QUAD_POINTS, INT_T_LOCATIONS, INT_T_WEIGHTS)
CALL Initialize_LG_Quadrature(NUM_P_QUAD_POINTS, INT_P_LOCATIONS, INT_P_WEIGHTS)

IF ( NUM_T_QUAD_POINTS == 1 ) THEN
    ! 1D Cheat !
    INT_T_WEIGHTS = 4.0_idp/pi
END IF


DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS
        INT_TP_WEIGHTS( (td-1)*NUM_P_QUAD_POINTS+pd ) = INT_T_WEIGHTS(td)*INT_P_WEIGHTS(pd)
    END DO
END DO



                                 !                                          !
                                !!   Initialize Lagrange Polynomial Tables  !!
                                 !                                          !
CALL Initialize_Lagrange_Poly_Tables()








                                 !                                          !
                                !!      Set Initial Mesh (Optional)         !!
                                 !                                          !
CALL Generate_Defined_Coarse_Mesh(R_Elements_Input, NUM_R_ELEMENTS, R_Coarsen_Factor,     &
                                  R_INNER, Input_Delta_R_Vector, rlocs, drlocs)
RADIAL_MESH_SET_FLAG = .TRUE.


CALL Generate_Defined_Coarse_Mesh(T_Elements_Input, NUM_T_ELEMENTS, T_Coarsen_Factor,     &
                                  0.0_idp, Input_Delta_T_Vector, tlocs, dtlocs)
THETA_MESH_SET_FLAG = .TRUE.


CALL Generate_Defined_Coarse_Mesh(P_Elements_Input, NUM_P_ELEMENTS, P_Coarsen_Factor,     &
                                  0.0_idp, Input_Delta_P_Vector, plocs, dplocs)
PHI_MESH_SET_FLAG = .TRUE.





IF ( DOMAIN_DIM == 1 ) THEN
     M_VALUES = 0
ELSE IF ( DOMAIN_DIM == 2 ) THEN
     M_VALUES = 0
ELSE IF ( DOMAIN_DIM == 3 ) THEN
     M_VALUES = (/(l,l=0,L_LIMIT,1)/)
END IF

LM_Location => FP_LM_Map



Ratio_T_BNDLperBLCK = NUM_T_ELEMENTS/( NUM_BLOCK_THETA_ROWS*NUM_LOC_T_ELEMENTS )
Ratio_P_BNDLperBLCK = NUM_P_ELEMENTS/( NUM_BLOCK_PHI_COLUMNS*NUM_LOC_P_ELEMENTS )
Ratio_BNDLperBLCK = Ratio_T_BNDLperBLCK * Ratio_P_BNDLperBLCK

CALL CREATE_POSEIDON_COMMUNICATORS( DATA_DIST_MODE )

CALL Initialize_Ylm_Tables()


CALL Initialize_Laplace_Matrices()



!CALL OUTPUT_SETUP_TABLE( nPROCS, R_Elements_Input, T_Elements_Input, P_Elements_Input )




END SUBROUTINE Poseidon_FP_Init_From_File










END MODULE Poseidon_FP_Initialization_Module

