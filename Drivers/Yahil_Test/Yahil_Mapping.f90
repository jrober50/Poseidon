  !################################################################################!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
PROGRAM Yahil_Mapping                                                               !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  !################################################################################!

USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Numbers_Module, &
            ONLY :  pi

USE Units_Module, &
            ONLY :  C_Square,        &
                    Set_Units,       &
                    Centimeter,      &
                    Gram

USE Initialization_Poseidon, &
            ONLY :  Initialize_Poseidon


USE Poseidon_XCFC_Interface_Module, &
            ONLY : Poseidon_Return_ExtrinsicCurvature


USE Variables_IO, &
            ONLY :  Write_Results_R_Samps,      &
                    Write_Results_T_Samps,      &
                    File_Suffix,                &
                    Report_Flags


USE Variables_MPI, &
            ONLY :  ierr

USE Variables_Functions, &
            ONLY :  Potential_Solution

USE FP_Functions_Results , &
            ONLY : Calc_1D_CFA_Values_FP

USE Variables_Yahil, &
            ONLY :  SelfSim_V_Switch

USE Poseidon_IO_Parameters, &
            ONLY :  Poseidon_Results_Dir


USE Functions_Mesh, &
            ONLY :  Create_3D_Mesh

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature_Locations

USE Functions_Mapping, &
            ONLY :  Map_From_X_Space

USE Poseidon_IO_Module, &
            ONLY :  Open_Run_Report_File,       &
                    Output_Final_Results,       &
                    Open_New_File,              &
                    OPEN_FILE_INQUISITION,      &
                    Output_Poseidon_Sources_3D

USE IO_Print_Results, &
            ONLY :  Print_Results

USE Poseidon_Main_Module, &
            ONLY :  Poseidon_Run,                                       &
                    Poseidon_Close

USE Initial_Guess_Module, &
            ONLY :  Poseidon_Input_Guess,           &
                    Poseidon_Init_FlatGuess

USE FP_IO_Module, &
            ONLY :  Output_FP_Timetable

USE Poseidon_Utilities_Module, &
            ONLY :  Poseidon_Calc_ADM_Mass,         &
                    Poseidon_Calc_ADM_Mass_Parts

USE Driver_SetSource_Module, &
            ONLY :  Driver_SetSource

USE Driver_SetBC_Module, &
            ONLY :  Driver_SetBC

USE Driver_SetGuess_Module, &
            ONLY :  Driver_SetGuess

!USE Timer_IO_Module, &
!            ONLY :  Output_Time_Report
!
!USE Timer_Routines_Module, &
!            ONLY :  TimerStart,     &
!                    TimerStop
!
!USE Timer_Variables_Module, &
!            ONLY :  Timer_Driver_SetSource,     &
!                    Timer_Driver_SetBC,         &
!                    Timer_Driver_SetGuess,      &
!                    Timer_Driver_Run,           &
!                    Timer_Driver_Extra

USE MPI


IMPLICIT NONE

!                                       !
!   Poseidon Initialization Variables   !
!                                       !

INTEGER                                                 ::  Dimension_Input

INTEGER                                                 ::  Mesh_Type
INTEGER                                                 ::  Solver_Type

CHARACTER(LEN = 1)                                      ::  Units_Input

INTEGER, DIMENSION(3)                                   ::  NE
INTEGER, DIMENSION(3)                                   ::  NQ
REAL(idp), DIMENSION(2)                                 ::  Domain_Edge
INTEGER                                                 ::  Num_DOF


REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  dx_c, dy_c, dz_c
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  x_e, y_e, z_e
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  x_c, y_c, z_c

LOGICAL                                                 ::  Verbose
LOGICAL                                                 ::  Print_Results_Flag

INTEGER,   DIMENSION(5)                                 ::  CFA_EQs




CHARACTER(LEN=10)                                       ::  Suffix_Input
CHARACTER(LEN=1)                                        ::  Suffix_Tail

REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Cur_R_Locs
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_R_Quad
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_T_Quad
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_P_Quad
REAL(idp)                                               ::  Left_Limit
REAL(idp)                                               ::  Right_Limit



INTEGER                                                 ::  myID

INTEGER                                                 ::  re
INTEGER                                                 ::  rq


INTEGER                                                 ::  Guess_Type

INTEGER, DIMENSION(:), ALLOCATABLE                      ::  RE_Table

INTEGER                                                 ::  RE_Index
INTEGER                                                 ::  RE_Index_Min
INTEGER                                                 ::  RE_Index_Max

INTEGER                                                 ::  Degree_Input
INTEGER                                                 ::  Degree_Min
INTEGER                                                 ::  Degree_Max

INTEGER                                                 ::  L_Limit_Input
INTEGER                                                 ::  L_Limit_Min
INTEGER                                                 ::  L_Limit_Max

INTEGER                                                 ::  M_Index
INTEGER                                                 ::  M_Index_Min
INTEGER                                                 ::  M_Index_Max

INTEGER                                                 ::  T_Index
INTEGER                                                 ::  T_Index_Min
INTEGER                                                 ::  T_Index_Max

REAL(idp)                                               ::  Kappa
REAL(idp)                                               ::  Gamma

INTEGER                                                 ::  Max_Iterations
REAL(idp)                                               ::  CC_Option

REAL(idp)                                               ::  Perturbation
REAL(idp)                                               ::  Offset

REAL(idp), DIMENSION(4)                                 ::  Yahil_Params

INTEGER, DIMENSION(1:8)                                 ::  Anderson_M_Values
CHARACTER(LEN=1), DIMENSION(1:10)                       ::  Letter_Table
REAL(idp), DIMENSION(1:6)                               ::  Time_Values
INTEGER, DIMENSION(1:2)                                 ::  L_Values

REAL(idp)                                               ::  ADM_Mass
REAL(idp)                                               ::  ADM_MassB
REAL(idp)                                               ::  ADM_Phys
REAL(idp)                                               ::  ADM_Curve

REAL(idp)                                               ::  Komar_Mass


REAL(idp), DIMENSION(:,:,:,:,:), ALLOCATABLE            ::  Output_Kij

CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

ALLOCATE( RE_Table(1:9) )



!############################################################!
!#                                                          #!
!#                      Test Parameters                     #!
!#                                                          #!
!############################################################!
Units_Input         = "G"
Solver_Type         = 3

RE_Table            = (/ 18, 128, 160, 240, 320, 400, 600, 256, 512 /)
Anderson_M_Values   = (/ 1, 2, 3, 4, 5, 10, 20, 50 /)
Time_Values         = (/ 51.0_idp, 15.0_idp, 5.0_idp, 1.50_idp, 0.5_idp, 0.05_idp /)
L_Values            = (/ 5, 10 /)

T_Index_Min         =  1
T_Index_Max         =  1

M_Index_Min         =  3
M_Index_Max         =  3

RE_Index_Min        =  2
RE_Index_Max        =  2

Degree_Min          =  1
Degree_Max          =  1

L_Limit_Min         =  0
L_Limit_Max         =  0

Guess_Type          =  1            !  1 = Flat, 2 = Educated, 3 = Perturbed Educated.
Perturbation        =  -0.01_idp    !  If Guess_Type == 3, rho is the perturbation parameter

Kappa               = 953946015514834.4
Gamma               = 1.30_idp

SelfSim_V_Switch    =  0

!Suffix_Tail         = "A"
!Convergence_Criteria = 1.0E-8_idp

Dimension_Input     = 3

Max_Iterations      = 10
CC_Option           = 1.0E-10_idp

Mesh_Type           = 4                         ! 1 = Uniform, 2 = Log, 3 = Split, 4 = Zoom
Domain_Edge(1)      = 0.0_idp                   ! Inner Radius (cm)
Domain_Edge(2)      = 1E9_idp                  ! Outer Radius (cm)


NE(1)               = 128 ! 1.5*128                       ! Number of Radial Elements
NE(2)               = 2                        ! Number of Theta Elements
NE(3)               = 2                        ! Number of Phi Elements

NQ(1)               = 5                        ! Number of Radial Quadrature Points
NQ(2)               = 5                        ! Number of Theta Quadrature Points
NQ(3)               = 1                         ! Number of Phi Quadrature Points


Verbose             = .TRUE.
!Verbose             = .FALSE.
Print_Results_Flag  = .TRUE.
!Print_Results_Flag  = .FALSE.
Suffix_Input        = "Params"

CFA_Eqs = (/ 1, 1, 1, 1, 1 /)


Letter_Table = (/ "A","B","C","D","E","F","G","H","I","J" /)


Write_Results_R_Samps = 256
Write_Results_T_Samps = 1

CALL Set_Units(Units_Input)



Domain_Edge = Domain_Edge*Centimeter


CALL Open_Run_Report_File()

DO M_Index = M_Index_Min, M_Index_Max
DO T_Index = T_Index_Min, T_Index_Max
DO RE_Index = RE_Index_Min, RE_Index_Max
DO Degree_Input = Degree_Min, Degree_Max
DO L_Limit_Input = L_Limit_Min, L_Limit_Max

    


    NE(1) = RE_Table(RE_Index)
!    NQ(3) = 2*L_Limit_Input + 1

    Suffix_Tail = Letter_Table(Mesh_Type)


    Num_DOF = NQ(1)*NQ(2)*NQ(3)

    ALLOCATE( Cur_R_Locs(1:NQ(1)) )
    ALLOCATE( Input_R_Quad(1:NQ(1)) )
    ALLOCATE( Input_T_Quad(1:NQ(2)) )
    ALLOCATE( Input_P_Quad(1:NQ(3)) )

    ALLOCATE( x_e(0:NE(1)), y_e(0:NE(2)), z_e(0:NE(3)) )
    ALLOCATE( x_c(1:NE(1)), y_c(1:NE(2)), z_c(1:NE(3)) )
    ALLOCATE( dx_c(1:NE(1)) )
    ALLOCATE( dy_c(1:NE(2)) )
    ALLOCATE( dz_c(1:NE(3)) )

    ALLOCATE( Output_Kij(NQ(1)*NQ(2)*NQ(3),NE(1),NE(2),NE(3),1:6) )


    Input_R_Quad = Initialize_LG_Quadrature_Locations(NQ(1))
    Input_T_Quad = Initialize_LG_Quadrature_Locations(NQ(2))
    Input_P_Quad = Initialize_LG_Quadrature_Locations(NQ(3))

    Left_Limit  = -0.50_idp
    Right_Limit = +0.50_idp

    Input_R_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_R_Quad)
    Input_T_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_T_Quad)
    Input_P_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_P_Quad)


    CALL Create_3D_Mesh( Mesh_Type,         &
                        Domain_Edge(1),    &
                        Domain_Edge(2),    &
                        NE(1),             &
                        NE(2),             &
                        NE(3),             &
                        x_e, x_c, dx_c,    &
                        y_e, y_c, dy_c,    &
                        z_e, z_c, dz_c,     &
                        Zoom = 1.032034864238313_idp )



    !############################################################!
    !#                                                          #!
    !#                   Initialize Poseidon                    #!
    !#                                                          #!
    !############################################################!

    CALL Initialize_Poseidon &
       (   Dimensions_Option           = Dimension_Input,               &
           FEM_Degree_Option           = Degree_Input,                  &
           L_Limit_Option              = L_Limit_Input,                 &
           Units_Option                = Units_Input,                   &
           Domain_Edge_Option          = Domain_Edge,                   &
           NE_Option                   = NE,                            &
           NQ_Option                   = NQ,                            &
           r_Option                    = x_e,                           &
           t_Option                    = y_e,                           &
           p_Option                    = z_e,                           &
           Suffix_Flag_Option          = Suffix_Input,                  &
           Suffix_Tail_Option          = Suffix_Tail,                   &
           Method_Flag_Option          = Solver_Type,                   &
           CFA_Eq_Flags_Option         = CFA_Eqs,                       &
           Max_Iterations_Option       = Max_Iterations,                &
           Convergence_Criteria_Option = CC_Option,                     &
           Anderson_M_Option           = Anderson_M_Values(M_Index),    &
           Verbose_Option              = Verbose,                       &
           WriteAll_Option             = .FALSE.,                       &
           Print_Setup_Option          = .TRUE.,                        &
           Write_Setup_Option          = .FALSE.,                       &
           Print_Results_Option        = Print_Results_Flag,            &
           Write_Results_Option        = .TRUE.,                        &
           Print_Timetable_Option      = .TRUE.,                       &
           Write_Timetable_Option      = .FALSE.,                       &
           Write_Sources_Option        = .FALSE.                        )





    !############################################################!
    !#                                                          #!
    !#               Create & Input Source Values               #!
    !#                                                          #!
    !############################################################!
!    CALL TimerStart( Timer_Driver_SetSource )

    Yahil_Params = [Time_Values(T_Index), Kappa, Gamma, 0.0_idp]
    CALL Driver_SetSource(  NE, NQ,             &
                            dx_c, x_e, y_e,     &
                            Input_R_Quad,       &
                            Input_T_Quad,       &
                            Input_P_Quad,       &
                            Left_Limit,         &
                            Right_Limit,        &
                            Solver_Type,        &
                            myID,               &
                            Yahil_Params        )

!    CALL TimerStop( Timer_Driver_SetSource )




    !############################################################!
    !#                                                          #!
    !#          Calculate and Set Boundary Conditions           #!
    !#                                                          #!
    !############################################################!
!    CALL TimerStart( Timer_Driver_SetBC )

    CALL Driver_SetBC( NE, x_e )

!    CALL TimerStop( Timer_Driver_SetBC )







    !############################################################!
    !#                                                          #!
    !#              Calculate and Set Initial Guess             #!
    !#                                                          #!
    !############################################################!
!    CALL TimerStart( Timer_Driver_SetGuess )

    CALL Driver_SetGuess(   NE, NQ,             &
                            dx_c, x_e,          &
                            Input_R_Quad,       &
                            Input_T_Quad,       &
                            Input_P_Quad,       &
                            Left_Limit,         &
                            Right_Limit,        &
                            Guess_Type          )

!    CALL TimerStop( Timer_Driver_SetGuess )


    !############################################################!
    !#                                                          #!
    !#                         Run Poseidon                     #!
    !#                                                          #!
    !############################################################!
!    CALL TimerStart( Timer_Driver_Run )

    Call Poseidon_Run()


!    CALL TimerStop( Timer_Driver_Run )

!    CALL TimerStart( Timer_Driver_Extra  )

!    CALL Poseidon_Calc_ADM_Mass( ADM_Mass )
!    CALL Poseidon_Calc_ADM_Mass_Parts( ADM_MassB, ADM_Phys, ADM_Curve)
!    CALL Poseidon_Calc_Komar_Mass( Komar_Mass )
    


    !############################################################!
    !#                                                          #!
    !#                       Output Results                     #!
    !#                                                          #!
    !############################################################!
    IF ((Print_Results_Flag .EQV. .TRUE.) .OR. (Verbose .EQV. .TRUE. )) THEN
        WRITE(*,'(A)')" Final Results "

!        WRITE(*,'(A,ES18.12,A)')"ADM Mass   : ", ADM_Mass / Gram, " grams"
!        WRITE(*,'(A,ES18.12,A)')"Komar Mass : ", Komar_Mass / Gram, " grams"

!        WRITE(*,'(A,ES18.12,A)')"ADM MassB: ", ADM_MassB / Gram, " grams"
!        WRITE(*,'(A,ES18.12,A)')"ADM Phys : ", ADM_Phys / Gram, " grams"
!        WRITE(*,'(A,ES18.12,A)')"ADM Crve : ", ADM_Curve / Gram, " grams"
!        WRITE(*,'(A,ES18.12,A)')"ADM Sum  : ", (ADM_Phys + ADM_Curve) / Gram, " grams"
        
        CALL Print_Results()
    END IF



!    CALL TimerStop( Timer_Driver_Extra )


    CALL Poseidon_Return_ExtrinsicCurvature( NE, NQ,                &
                                             Input_R_Quad,          &
                                             Input_T_Quad,          &
                                             Input_P_Quad,          &
                                             Left_Limit,            &
                                             Right_Limit,           &
                                             Output_Kij             )



    !############################################################!
    !#                                                          #!
    !#                      Close Poseidon                      #!
    !#                                                          #!
    !############################################################!

    CALL Poseidon_Close()

    DEALLOCATE( Cur_R_Locs )
    DEALLOCATE( Input_R_Quad )
    DEALLOCATE( Input_T_Quad )
    DEALLOCATE( Input_P_Quad )
    DEALLOCATE( x_e, y_e, z_e )
    DEALLOCATE( x_c, y_c, z_c )
    DEALLOCATE( dx_c, dy_c, dz_c )

    DEALLOCATE( Output_Kij )


END DO ! L_Limit
END DO ! Degree_Index
END DO ! RE_Index
END DO ! T_Index
END DO ! M_Index




CALL MPI_FINALIZE(ierr)

CONTAINS






END PROGRAM Yahil_Mapping



