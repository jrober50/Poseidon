  !################################################################################!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
PROGRAM Yahil_Movie                                                                 !##!
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

USE Poseidon_Units_Module, &
            ONLY :  C_Square,        &
                    Set_Units,       &
                    Centimeter,      &
                    Gram

USE Initialization_Poseidon, &
            ONLY :  Initialize_Poseidon

USE Variables_IO, &
            ONLY :  Write_Results_R_Samps,      &
                    Write_Results_T_Samps,      &
                    File_Suffix,                &
                    Report_Flags,               &
                    iRF_Time

USE Variables_MPI, &
            ONLY :  ierr

USE Variables_External, &
            ONLY :  SelfSim_V_Switch

USE Functions_Mesh, &
            ONLY :  Create_3D_Mesh

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature_Locations

USE Maps_X_Space, &
            ONLY :  Map_From_X_Space

USE IO_Write_Final_Results, &
            ONLY :  Write_Final_Results

USE IO_Print_Results, &
            ONLY :  Print_Results

USE Poseidon_Main_Module, &
            ONLY :  Poseidon_Run,                                       &
                    Poseidon_Close

USE Poseidon_Utilities_Module, &
            ONLY :  Poseidon_Calc_ADM_Mass,         &
                    Poseidon_Calc_ADM_Mass_Parts,   &
                    Poseidon_Calc_Komar_Mass

USE Driver_SetSource_Module, &
            ONLY :  Driver_SetSource

USE Driver_SetBC_Module, &
            ONLY :  Driver_SetBC

USE Driver_SetGuess_Module, &
            ONLY :  Driver_SetGuess

USE Timer_Routines_Module, &
            ONLY :  TimerStart,     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Driver_SetSource,     &
                    Timer_Driver_SetBC,         &
                    Timer_Driver_SetGuess,      &
                    Timer_Driver_Run,           &
                    Timer_Driver_Extra



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

INTEGER                                                 ::  Ylm_L_Limit
INTEGER                                                 ::  FEM_Degree

REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  dx_c, dy_c, dz_c
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  x_e, y_e, z_e
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  x_c, y_c, z_c

LOGICAL                                                 ::  Verbose
LOGICAL                                                 ::  Print_Results_Flag

INTEGER,   DIMENSION(5)                                 ::  CFA_EQs

CHARACTER(LEN=10)                                       ::  Suffix_Input
CHARACTER(LEN=1)                                        ::  Suffix_Tail

REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_R_Quad
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_T_Quad
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_P_Quad
REAL(idp)                                               ::  Left_Limit
REAL(idp)                                               ::  Right_Limit

INTEGER                                                 ::  myID

INTEGER                                                 ::  Guess_Type

INTEGER                                                 ::  T_Index
REAL(idp)                                               ::  Time_Start
REAL(idp)                                               ::  Time_Stop
INTEGER                                                 ::  Time_Steps
REAL(idp)                                               ::  dt
REAL(idp)                                               ::  Time

REAL(idp)                                               ::  Kappa
REAL(idp)                                               ::  Gamma

INTEGER                                                 ::  AMReX_Levels
INTEGER                                                 ::  Max_Iterations
REAL(idp)                                               ::  CC_Option

REAL(idp)                                               ::  Perturbation
REAL(idp)                                               ::  Offset

REAL(idp), DIMENSION(4)                                 ::  Yahil_Params

INTEGER                                                 ::  Anderson_M_Value
CHARACTER(LEN=1), DIMENSION(1:10)                       ::  Letter_Table



REAL(idp)                                               ::  ADM_Mass
REAL(idp)                                               ::  ADM_MassB
REAL(idp)                                               ::  ADM_Phys
REAL(idp)                                               ::  ADM_Curve

REAL(idp)                                               ::  Komar_Mass
REAL(idp), DIMENSION(:,:,:,:,:), ALLOCATABLE            ::  Output_Kij


CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)





!############################################################!
!#                                                          #!
!#                      Test Parameters                     #!
!#                                                          #!
!############################################################!
Units_Input         = "G"
Solver_Type         = 3

Time_Start          = 15.0 ! ms
Time_Stop           = 0.015  ! ms
Time_Steps          = 240

FEM_Degree          = 1
Ylm_L_Limit         = 0


NE(1)               = 256 ! 1.5*128            ! Number of Radial Elements
NE(2)               = 1                        ! Number of Theta Elements
NE(3)               = 1                        ! Number of Phi Elements

NQ(1)               = 5                        ! Number of Radial Quadrature Points
NQ(2)               = 1                        ! Number of Theta Quadrature Points
NQ(3)               = 2*Ylm_L_Limit + 1        ! Number of Phi Quadrature Points
Num_DOF             = NQ(1)*NQ(2)*NQ(3)


!Verbose             = .TRUE.
Verbose             = .FALSE.
Print_Results_Flag  = .TRUE.
!Print_Results_Flag  = .FALSE.


Anderson_M_Value    =  3
AMReX_Levels        =  0


Mesh_Type           = 2                         ! 1 = Uniform, 2 = Log, 3 = Split, 4 = Zoom
Domain_Edge(1)      = 0.0_idp                   ! Inner Radius (cm)
Domain_Edge(2)      = 1E9_idp                   ! Outer Radius (cm)


Guess_Type          =  1            !  1 = Flat, 2 = Educated, 3 = Perturbed Educated.
Perturbation        =  -0.01_idp    !  If Guess_Type == 3, rho is the perturbation parameter


Kappa               = 953946015514834.4
Gamma               = 1.30_idp
SelfSim_V_Switch    =  0


Dimension_Input     = 3

Max_Iterations      = 20
CC_Option           = 1.0E-15_idp


CFA_Eqs = (/ 1, 1, 1, 1, 1 /)

Suffix_Input    = "Frame"
!Suffix_Input    = "Params"
Suffix_Tail     = Letter_Table(1)
Letter_Table    = (/ "A","B","C","D","E","F","G","H","I","J" /)


Write_Results_R_Samps = 256
Write_Results_T_Samps = 1


CALL Set_Units(Units_Input)


Domain_Edge = Domain_Edge*Centimeter


ALLOCATE( Input_R_Quad(1:NQ(1)) )
ALLOCATE( Input_T_Quad(1:NQ(2)) )
ALLOCATE( Input_P_Quad(1:NQ(3)) )


ALLOCATE( x_e(0:NE(1)), y_e(0:NE(2)), z_e(0:NE(3)) )
ALLOCATE( x_c(1:NE(1)), y_c(1:NE(2)), z_c(1:NE(3)) )
ALLOCATE( dx_c(1:NE(1)) )
ALLOCATE( dy_c(1:NE(2)) )
ALLOCATE( dz_c(1:NE(3)) )


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
                    Zoom = 1.032034864238313_idp,   &
                    Levels_In = AMReX_Levels       )











!############################################################!
!#                                                          #!
!#                   Initialize Poseidon                    #!
!#                                                          #!
!############################################################!








!############################################################!
!#                                                          #!
!#                    Start Time Stepping                   #!
!#                                                          #!
!############################################################!
IF (Time_Steps == 0 ) THEN
    dt = 0.0_idp
ELSE
    dt = (Time_Stop - Time_Start)/(Time_Steps-1)
END IF


DO T_Index = 1,Time_Steps
    Time = Time_Start + (T_Index-1)*dt
    WRITE(*,'(A,I2.2,A,F4.2)') "Starting Loop ",T_Index," at time, ",Time




    CALL Initialize_Poseidon &
    (   Dimensions_Option           = Dimension_Input,               &
        FEM_Degree_Option           = FEM_Degree,                    &
        L_Limit_Option              = YLM_L_Limit,                   &
        Units_Option                = Units_Input,                   &
        Domain_Edge_Option          = Domain_Edge,                   &
        NE_Option                   = NE,                            &
        NQ_Option                   = NQ,                            &
        r_Option                    = x_e,                           &
        t_Option                    = y_e,                           &
        p_Option                    = z_e,                           &
        Suffix_Flag_Option          = Suffix_Input,                  &
        Frame_Option                = T_Index,                       &
        Method_Flag_Option          = Solver_Type,                   &
        CFA_Eq_Flags_Option         = CFA_Eqs,                       &
        Max_Iterations_Option       = Max_Iterations,                &
        Convergence_Criteria_Option = CC_Option,                     &
        Anderson_M_Option           = Anderson_M_Value,              &
        Verbose_Option              = Verbose,                       &
        WriteAll_Option             = .FALSE.,                       &
        Print_Setup_Option          = .TRUE.,                        &
        Write_Setup_Option          = .FALSE.,                       &
        Print_Results_Option        = Print_Results_Flag,            &
        Write_Results_Option        = .TRUE.,                        &
        Print_Timetable_Option      = .FALSE.,                       &
        Write_Timetable_Option      = .FALSE.,                       &
        Write_Sources_Option        = .FALSE.                        )


    !############################################################!
    !#                                                          #!
    !#               Create & Input Source Values               #!
    !#                                                          #!
    !############################################################!

    CALL TimerStart( Timer_Driver_SetSource )

    Yahil_Params = [Time, Kappa, Gamma, 0.0_idp]
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

    CALL TimerStop( Timer_Driver_SetSource )




    !############################################################!
    !#                                                          #!
    !#          Calculate and Set Boundary Conditions           #!
    !#                                                          #!
    !############################################################!
    CALL TimerStart( Timer_Driver_SetBC )

    CALL Driver_SetBC( NE, x_e )

    CALL TimerStop( Timer_Driver_SetBC )







    !############################################################!
    !#                                                          #!
    !#              Calculate and Set Initial Guess             #!
    !#                                                          #!
    !############################################################!
    CALL TimerStart( Timer_Driver_SetGuess )

    CALL Driver_SetGuess(   NE, NQ,             &
                            dx_c, x_e,          &
                            Input_R_Quad,       &
                            Input_T_Quad,       &
                            Input_P_Quad,       &
                            Left_Limit,         &
                            Right_Limit,        &
                            Guess_Type          )

    CALL TimerStop( Timer_Driver_SetGuess )


    !############################################################!
    !#                                                          #!
    !#                         Run Poseidon                     #!
    !#                                                          #!
    !############################################################!
    CALL TimerStart( Timer_Driver_Run )

    Call Poseidon_Run()

    CALL TimerStop( Timer_Driver_Run )

    


    !############################################################!
    !#                                                          #!
    !#                       Output Results                     #!
    !#                                                          #!
    !############################################################!
    CALL TimerStart( Timer_Driver_Extra  )
    IF ((Print_Results_Flag .EQV. .TRUE.) .OR. (Verbose .EQV. .TRUE. )) THEN
        WRITE(*,'(A)')" Final Results "
        
        CALL Print_Results()
        CALL Write_Final_Results()

    END IF

    CALL TimerStop( Timer_Driver_Extra )







    !############################################################!
    !#                                                          #!
    !#                      Close Poseidon                      #!
    !#                                                          #!
    !############################################################!
    CALL Poseidon_Close()




END DO ! T_Index




DEALLOCATE( Input_R_Quad )
DEALLOCATE( Input_T_Quad )
DEALLOCATE( Input_P_Quad )
DEALLOCATE( x_e, y_e, z_e )
DEALLOCATE( x_c, y_c, z_c )
DEALLOCATE( dx_c, dy_c, dz_c )

CALL MPI_FINALIZE(ierr)

CONTAINS






END PROGRAM Yahil_Movie



