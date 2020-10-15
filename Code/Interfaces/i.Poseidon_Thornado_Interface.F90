   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Thornado_Interface                                                  !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the interfacing subroutine that will handle action between         !##!
!##!    Poseidon and CHIMERA.                                                       !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!



USE Poseidon_Constants_Module, &
            ONLY :  idp, pi


USE Units_Module, &
            ONLY :  Set_Units, Grav_Constant_G


USE DRIVER_Parameters,  &
            ONLY :  myID,                                               &
                    myID_Theta,                                         &
                    myID_Phi


USE Poseidon_Variables_Module, &
            ONLY :  NUM_R_ELEMENTS,                                     &
                    NUM_T_ELEMENTS,                                     &
                    NUM_P_ELEMENTS,                                     &
                    R_INNER, R_OUTER,                                   &
                    rlocs,                                              &
                    tlocs,                                              &
                    plocs


USE Poseidon_Main_Module, &
            ONLY :  Poseidon_Initialize,                                &
                    Poseidon_Run,                                       &
                    Poseidon_Close,                                     &
                    Poseidon_Set_Mesh,                                  &
                    Poseidon_CFA_Set_Uniform_Boundary_Conditions

USE Poseidon_Initialization_Module, &
            ONLY :  Poseidon_Initialize_1D,                             &
                    Poseidon_Initialize_From_File


USE Poseidon_Additional_Functions_Module, &
            ONLY :  Calc_1D_CFA_Values


USE Poseidon_Source_Module, &
            ONLY :   Poseidon_Input_Sources

USE Thornado_Guess_BC_Module,   &
            ONLY :   CALC_BC_VALUES,                                    &
                     Initialize_Flat_Space_Guess_Values


use mpi








IMPLICIT NONE

CONTAINS


 !+101+############################################################################!
!                                                                                   !
!                    Poseidon_Init                                                  !
!                                                                                   !
 !#################################################################################!
SUBROUTINE Poseidon_Init( mode,                                                                   &
                          FEM_Degree, SH_Limit,                                                   &
                          Inner_Radius, Outer_Radius,                                             &
                          R_Elements_Input, T_Elements_Input, P_Elements_Input,                   &
                          Local_R_Elements_Input, Local_T_Elements_Input, Local_P_Elements_Input, &
                          Num_R_Quad_Input, Num_T_Quad_Input, Num_P_Quad_Input,                   &
                          Input_Delta_R_Vector, Input_Delta_T_Vector, Input_Delta_P_Vector        )


INTEGER, INTENT(IN)                                                             ::  mode,                   &
                                                                                    FEM_Degree,             &
                                                                                    SH_Limit,               &
                                                                                    R_Elements_Input,       &
                                                                                    T_Elements_Input,       &
                                                                                    P_Elements_Input,       &
                                                                                    Local_R_Elements_Input, &
                                                                                    Local_T_Elements_Input, &
                                                                                    Local_P_Elements_Input, &
                                                                                    Num_R_Quad_Input,       &
                                                                                    Num_T_Quad_Input,       &
                                                                                    Num_P_Quad_Input

REAL(KIND = idp), INTENT(IN)                                                    ::  Inner_Radius,           &
                                                                                    Outer_Radius


REAL(KIND = idp), DIMENSION(1:R_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_R_Vector
REAL(KIND = idp), DIMENSION(1:T_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_T_Vector
REAL(KIND = idp), DIMENSION(1:P_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_P_Vector



CALL Set_Units( "G" )

CALL Poseidon_Initialize_1D(FEM_Degree,                     & ! FEM Degree
                            SH_Limit,                       & ! Spherical Harmonic Expansion L limit
                            Inner_Radius,                   & ! Inner_Radius
                            Outer_Radius,                   & ! Outer_Radius
                            R_Elements_Input,               & ! NUM_R_ELEMENTS
                            T_Elements_Input,               & ! NUM_T_ELEMENTS
                            P_Elements_Input,               & ! NUM_P_ELEMENTS
                            Local_R_Elements_Input,         & ! NUM_LOC_R_ELEMENTS
                            Local_T_Elements_Input,         & ! NUM_LOC_T_ELEMENTS
                            Local_P_Elements_Input,         & ! NUM_LOC_P_ELEMENTS
                            Num_R_Quad_Input,                     & ! NUM_R_QUAD_POINTS
                            Num_T_Quad_Input,                     & ! NUM_T_QUAD_POINTS
                            Num_P_Quad_Input,                     & ! NUM_P_QUAD_POINTS
                            Input_Delta_R_Vector,           & ! Delta_R_Vector
                            Input_Delta_T_Vector,           & ! Delta_T_Vector
                            Input_Delta_P_Vector            ) ! Delta_P_Vector

PRINT*,"POSIEDON INIT COMPLETE, Poseidon_Init"


END SUBROUTINE Poseidon_Init




 !+201+############################################################################!
!                                                                                   !
!                    Poseidon_CFA_1D                                                !
!                                                                                   !
 !#################################################################################!
SUBROUTINE Poseidon_CFA_1D( First_Run_Flag,                                     &
                            Num_Input_R_Quad,                                   &
                            Input_R_Quad,                                       &
                            Left_Limit, Right_Limit,                            &
                            Local_E, Local_S, Local_Si,                         &
                            CFA_Lapse, CFA_ConFactor, CFA_Shift                 )


INTEGER, INTENT(IN)                                             ::  First_Run_Flag

INTEGER, INTENT(IN)                                             ::  Num_Input_R_Quad
REAL(KIND = idp), DIMENSION(1:Num_Input_R_Quad), INTENT(IN)     ::  Input_R_Quad


REAL(KIND = idp), INTENT(IN)                                    ::  Left_Limit, &
                                                                    Right_Limit

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:Num_Input_R_Quad,                    &
                                         0:NUM_R_ELEMENTS-1,                    &
                                         0:NUM_T_ELEMENTS-1,                    &
                                         0:NUM_P_ELEMENTS-1  )  ::  Local_E,    &
                                                                    Local_S,    &
                                                                    Local_Si

REAL(KIND = idp), INTENT(OUT), DIMENSION(1:Num_Input_R_Quad,                    &
                                         0:NUM_R_ELEMENTS-1,                    &
                                         0:NUM_T_ELEMENTS-1,                    &
                                         0:NUM_P_ELEMENTS-1  )  ::  CFA_Lapse,  &
                                                                    CFA_ConFactor,  &
                                                                    CFA_Shift


REAL(KIND = idp)                                                ::  Lapse_BC,   &
                                                                    ConFact_BC, &
                                                                    Shift_BC

CHARACTER(LEN=1), DIMENSION(1:5)                            ::  INNER_BC_TYPES, OUTER_BC_TYPES
REAL(KIND = idp), DIMENSION(1:5)                            ::  INNER_BC_VALUES, OUTER_BC_VALUES


! Set Source Values

CALL Poseidon_Input_Sources( myID, myID_Theta, myID_Phi,                                &
                             Local_E, Local_S, Local_Si,                                &
                             NUM_R_ELEMENTS, NUM_T_ELEMENTS, NUM_P_ELEMENTS,            &
                             Num_Input_R_Quad,1,1,                                      &
                             Input_R_Quad, (/ pi/2.0_idp /), (/ 0.0_idp /),                  &
                             Left_Limit, Right_Limit                                    )

PRINT*,"AFTER Poseidon_Input_Sources"


!
! Set Boundary Conditions
!
CALL CALC_BC_VALUES(Num_Input_R_Quad, Left_Limit, Right_Limit,              &
                            Input_R_Quad, rlocs, NUM_R_ELEMENTS,            &
                            Local_E, Local_Si,                              &
                            Lapse_BC, ConFact_BC, Shift_BC                  )

PRINT*,"AFTER CALC_BC_VALUES"


INNER_BC_TYPES = (/"N", "N","N","N","N"/)
OUTER_BC_TYPES = (/"D", "D","D","D","D"/)

INNER_BC_VALUES = (/0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp /)
OUTER_BC_VALUES = (/ Lapse_BC, ConFact_BC , Shift_BC, 0.0_idp, 0.0_idp /)

CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("I", INNER_BC_TYPES, INNER_BC_VALUES)
CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("O", OUTER_BC_TYPES, OUTER_BC_VALUES)





!
! Set Initial Guess
!
IF ( First_Run_Flag == 1 ) THEN

    ! Create Initial Guess
    CALL Initialize_Flat_Space_Guess_Values()
END IF







!
! Run
!
CALL Poseidon_Run()



PRINT*,"AFTER Poseidon_Run"



! Calculate Results
CAll Calc_1D_CFA_Values( Num_R_Elements, Num_Input_R_Quad, Input_R_Quad,   &
                         Left_Limit, Right_Limit,                          &
                         CFA_Lapse, CFA_ConFactor, CFA_Shift               )






END SUBROUTINE Poseidon_CFA_1D













 !+301+############################################################################!
!                                                                                   !
!                    Poseidon_Shutdown                                              !
!                                                                                   !
 !#################################################################################!
SUBROUTINE Poseidon_Shutdown()

CALL Poseidon_Close()

END SUBROUTINE Poseidon_Shutdown










END MODULE Poseidon_Thornado_Interface
