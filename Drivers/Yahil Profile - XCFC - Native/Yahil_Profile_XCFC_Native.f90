  !################################################################################!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
PROGRAM Yahil_Profile_XCFC_Native                                                   !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  !################################################################################!

USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Units_Module, &
            ONLY :  Set_Units

USE Poseidon_Letters_Module, &
            ONLY :  Letter_Table_Upper

USE Poseidon_Interface_Initialization, &
            ONLY :  Initialize_Poseidon

USE Poseidon_Interface_Run, &
            ONLY :  Poseidon_Run

USE Poseidon_Interface_Close, &
            ONLY :  Poseidon_Close

USE Poseidon_Return_Routines_Module, &
            ONLY : Poseidon_Return_Extrinsic_Curvature

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

USE Driver_SetSource_Module, &
            ONLY :  Driver_SetSource

USE Driver_SetBC_Module, &
            ONLY :  Driver_SetBC

USE Driver_SetGuess_Module, &
            ONLY :  Driver_SetGuess

USE Driver_ConFactor_Loop_Module, &
            ONLY :  Driver_ConFactor_Loop

USE Timer_Routines_Module, &
            ONLY :  TimerStart,     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Driver_SetSource,     &
                    Timer_Driver_SetBC,         &
                    Timer_Driver_SetGuess,      &
                    Timer_Driver_Run,           &
                    Timer_Driver_Extra

USE ADM_Mass_Module, &
            ONLY :  Calc_ADM_Mass

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
LOGICAL                                                 ::  Print_Setup_Flag
LOGICAL                                                 ::  Print_Time_Flag

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

INTEGER                                                 ::  AMReX_Levels
INTEGER                                                 ::  Max_Iterations
REAL(idp)                                               ::  CC_Option


INTEGER                                                 ::  External_Iter_Max
REAL(idp)                                               ::  External_Tolerance

REAL(idp)                                               ::  Perturbation

REAL(idp), DIMENSION(4)                                 ::  Yahil_Params

INTEGER, DIMENSION(1:8)                                 ::  Anderson_M_Values
REAL(idp), DIMENSION(1:7)                               ::  Time_Values
INTEGER, DIMENSION(1:2)                                 ::  L_Values

REAL(idp)                                               ::  ADM_Mass
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


RE_Table            = (/ 10, 128, 256, 384, 512, 640, 768, 896, 1024, 4096 /)
Anderson_M_Values   = (/ 1, 2, 3, 4, 5, 10, 20, 50 /)
Time_Values         = (/ 51.0_idp, 15.0_idp, 5.0_idp, 1.50_idp, 0.51_idp, 0.275_idp, 0.15_idp /)
L_Values            = (/ 5, 10 /)

T_Index_Min         =  1
T_Index_Max         =  1

M_Index_Min         =  1
M_Index_Max         =  1

RE_Index_Min        =  5
RE_Index_Max        =  5

Degree_Min          =  1
Degree_Max          =  1

L_Limit_Min         =  0
L_Limit_Max         =  0

AMReX_Levels        =  0


External_Iter_Max   =  500
External_Tolerance  =  1.0E-10

Guess_Type          =  1            !  1 = Flat, 2 = Educated, 3 = Perturbed Educated.
Perturbation        =  -0.01_idp    !  If Guess_Type == 3, rho is the perturbation parameter

Kappa               = 953946015514834.4
Gamma               = 1.30_idp

SelfSim_V_Switch    =  0

!Suffix_Tail         = "A"
!Convergence_Criteria = 1.0E-8_idp

Dimension_Input     = 3

Max_Iterations      = 100
CC_Option           = 1.0E-10_idp

Mesh_Type           = 1                         ! 1 = Uniform, 2 = Log, 3 = Split, 4 = Zoom
Domain_Edge(1)      = 0.0_idp                   ! Inner Radius (cm)
Domain_Edge(2)      = 8.0E8_idp                   ! Outer Radius (cm)




NE(1)               = 128 ! 1.5*128            ! Number of Radial Elements
NE(2)               = 1                        ! Number of Theta Elements
NE(3)               = 1                        ! Number of Phi Elements

NQ(1)               = 5                        ! Number of Radial Quadrature Points
NQ(2)               = 1                        ! Number of Theta Quadrature Points
NQ(3)               = 1                        ! Number of Phi Quadrature Points


Verbose             = .TRUE.
!Verbose             = .FALSE.

Print_Results_Flag  = .TRUE.
!Print_Results_Flag  = .FALSE.

Print_Setup_Flag    = .TRUE.
!Print_Setup_Flag    = .FALSE.

Print_Time_Flag     = .TRUE.
!Print_Time_Flag     = .FALSE.

Suffix_Input        = "Params"

CFA_Eqs = (/ 1, 1, 1, 1, 1 /)





CALL Set_Units(Units_Input)





DO M_Index = M_Index_Min, M_Index_Max
DO T_Index = T_Index_Min, T_Index_Max
DO RE_Index = RE_Index_Min, RE_Index_Max
DO Degree_Input = Degree_Min, Degree_Max
DO L_Limit_Input = L_Limit_Min, L_Limit_Max

    
    IF ( Mesh_Type == 6 ) THEN
        NE(1) = INT(RE_Table(RE_Index)*(1.0_idp + REAL((AMReX_Levels-1),idp)/2.0_idp))
    ELSE
        NE(1) = RE_Table(RE_Index)
    END IF
    NQ(3) = 2*L_Limit_Input + 1

    Suffix_Tail = Letter_Table_Upper(T_Index)


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
                        Zoom = 1.032034864238313_idp,   &
                        Levels_In = AMReX_Levels       )

    


    !############################################################!
    !#                                                          #!
    !#                   Initialize Poseidon                    #!
    !#                                                          #!
    !############################################################!
CALL Initialize_Poseidon &
   (   Dimensions_Option            = Dimension_Input,              &
       FEM_Degree_Option            = Degree_Input,                 &
       L_Limit_Option               = L_Limit_Input,                &
       Source_NE                    = NE,                           &
       Domain_Edge_Option           = Domain_Edge,                  &
       Source_NQ                    = NQ,                           &
       Source_xL                    = [Left_Limit, Right_Limit],    &
       Source_RQ_xlocs              = Input_R_Quad,                 &
       Source_TQ_xlocs              = Input_T_Quad,                 &
       Source_PQ_xlocs              = Input_P_Quad,                 &
       Source_Units                 = Units_Input,                  &
       Source_Radial_Boundary_Units = "cm",                         &
       Integration_NQ_Option        = NQ,                           &
       Source_R_Option              = x_e,                          &
       Source_T_Option              = y_e,                          &
       Source_P_Option              = z_e,                          &
       Method_Flag_Option           = Solver_Type,                  &
       CFA_Eq_Flags_Option          = CFA_Eqs,                      &
       Max_Iterations_Option        = Max_Iterations,               &
       Convergence_Criteria_Option  = CC_Option,                    &
       Anderson_M_Option            = Anderson_M_Values(M_Index),   &
       Verbose_Option               = Verbose,                      &
       WriteAll_Option              = .FALSE.,                      &
       Print_Setup_Option           = Print_Setup_Flag,             &
       Write_Setup_Option           = .TRUE.,                       &
       Print_Results_Option         = Print_Results_Flag,           &
       Write_Results_Option         = .TRUE.,                       &
       Print_Timetable_Option       = Print_Time_Flag,              &
       Write_Timetable_Option       = .TRUE.,                       &
       Write_Sources_Option         = .TRUE.,                       &
       Print_Condition_Option       = .FALSE.,                      &
       Write_Condition_Option       = .TRUE.,                       &
       Suffix_Flag_Option           = Suffix_Input,                 &
       Suffix_Tail_Option           = Suffix_Tail                   )




    !############################################################!
    !#                                                          #!
    !#               Create & Input Source Values               #!
    !#                                                          #!
    !############################################################!

    CALL TimerStart( Timer_Driver_SetSource )

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

    CALL TimerStart( Timer_Driver_Extra  )



!    CALL Driver_ConFactor_Loop( NE, NQ,                    &
!                                dx_c, x_e, y_e,             &
!                                Input_R_Quad, Input_T_Quad, Input_P_Quad,    &
!                                Left_Limit, Right_Limit,   &
!                                Solver_Type,               &
!                                myID,                      &
!                                Yahil_Params,              &
!                                External_Tolerance,         &
!                                External_Iter_Max          )
!


    !############################################################!
    !#                                                          #!
    !#                       Output Results                     #!
    !#                                                          #!
    !############################################################!
    IF ((Print_Results_Flag .EQV. .TRUE.) .OR. (Verbose .EQV. .TRUE. )) THEN
        WRITE(*,'(A)')" Final Results "
        
        CALL Print_Results()
        CALL Write_Final_Results()

        CALL Calc_ADM_Mass(ADM_Mass)
        PRINT*,"ADM Mass",ADM_Mass
    END IF




    CALL TimerStop( Timer_Driver_Extra )


    CALL Poseidon_Return_Extrinsic_Curvature( NE, NQ,                &
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






END PROGRAM Yahil_Profile_XCFC_Native



!
