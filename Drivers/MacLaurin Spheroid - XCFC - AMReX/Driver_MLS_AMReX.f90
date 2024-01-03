  !################################################################################!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
PROGRAM Driver_Main                                                                 !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  !################################################################################!

USE amrex_base_module
USE amrex_init_module, &
                ONLY:  amrex_init
USE amrex_amrcore_module, &
            ONLY: amrex_amrcore_init

USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Numbers_Module, &
            ONLY :  pi

USE Poseidon_Letters_Module, &
            ONLY :  Letter_Table_Upper

USE Poseidon_Units_Module, &
            ONLY :  Set_Units,      &
                    Centimeter,     &
                    rho_Units

USE Poseidon_Interface_Initialization, &
            ONLY :  Initialize_Poseidon

USE Variables_IO, &
            ONLY :  Write_Results_R_Samps,      &
                    Write_Results_T_Samps,      &
                    Write_Results_P_Samps

USE Variables_MPI, &
            ONLY :  ierr

USE Functions_Mesh, &
            ONLY :  Create_3D_Mesh

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature_Locations

USE Maps_X_Space, &
            ONLY :  Map_From_X_Space


USE IO_Write_Final_Results, &
            ONLY :  Write_Final_Results,        &
                    Output_2D_Results

USE Poseidon_Interface_Run, &
            ONLY :  Poseidon_Run
            
USE Poseidon_Interface_Close, &
            ONLY :  Poseidon_Close

USE Driver_SetSource_Module, &
            ONLY:  Driver_SetSource
            
USE Driver_InitSource_Module, &
            ONLY :  Driver_InitSource

USE Driver_SetBC_Module, &
            ONLY:  Driver_SetBC

USE Driver_SetGuess_Module, &
            ONLY:  Driver_SetGuess
            
USE Driver_ConFactor_Loop_Module, &
            ONLY :  Driver_ConFactor_Loop

USE Poseidon_AMReX_Input_Parsing_Module, &
            ONLY : Init_AMReX_Parameters

USE Variables_MPI, &
            ONLY :  myID_Poseidon,      &
                    MasterID_Poseidon,  &
                    nPROCS_Poseidon,    &
                    Poseidon_Comm_World

USE Variables_Driver_AMReX, &
            ONLY :  xL,                 &
                    xR,                 &
                    nCells,             &
                    MaxGridSizeX,       &
                    nLevels

USE Driver_Variables, &
            ONLY :  Driver_NQ,          &
                    Driver_RQ_xLocs,    &
                    Driver_xL

USE Variables_AMReX_Core, &
            ONLY :  MF_Source,          &
                    AMReX_Num_Levels
                    
USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    Num_T_Elements,         &
                    Num_P_Elements,         &
                    rlocs,                  &
                    drlocs
                    
USE External_MLS_Profile_Module, &
            ONLY :  Set_MLS_Parameters

USE External_IO_Test_Results_Module, &
            ONLY :  Print_MacLaurin_Error
            
USE MPI


IMPLICIT NONE

!                                       !
!   Poseidon Initialization Variables   !
!                                       !
CHARACTER(LEN = 1)                                      ::  Units_Input

INTEGER,   DIMENSION(5)                                 ::  CFA_EQs

LOGICAL                                                 ::  Verbose
LOGICAL                                                 ::  Print_Results_Flag
LOGICAL                                                 ::  Print_Setup_Flag
LOGICAL                                                 ::  Print_Time_Flag
LOGICAL                                                 ::  Print_Cond_Flag


CHARACTER(LEN=10)                                       ::  Suffix_Input
CHARACTER(LEN=4)                                        ::  Suffix_Tail

INTEGER, DIMENSION(3)                                   ::  NQ
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_R_Quad
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_T_Quad
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_P_Quad
REAL(idp)                                               ::  Left_Limit
REAL(idp)                                               ::  Right_Limit


INTEGER                                                 ::  myID, nPROCS

INTEGER                                                 ::  M_Index
INTEGER                                                 ::  M_Index_Min
INTEGER                                                 ::  M_Index_Max

INTEGER                                                 ::  IRL
INTEGER                                                 ::  IFL

INTEGER, DIMENSION(1:8)                                 ::  Anderson_M_Values

INTEGER                                                 ::  i

REAL(idp)                                               ::  SemiMajor_Axis
REAL(idp)                                               ::  SemiMinor_Axis
REAL(idp)                                               ::  Rho

CHARACTER(LEN=7)                                        ::  SphereType

REAL(idp)                                               ::  Radius_Earth  = 6.371E8
REAL(idp)                                               ::  Radius_Sun    = 6.957e10
REAL(idp)                                               ::  Density_Earth = 5.5
REAL(idp)                                               ::  Density_Sun   = 1.41


!############################################################!
!#                                                          #!
!#                      Test Parameters                     #!
!#                                                          #!
!############################################################!
Units_Input         = "G"

Anderson_M_Values   = (/ 1, 2, 3, 4, 5, 10, 20, 50 /)

M_Index_Min         =  3
M_Index_Max         =  3

IFL                 =  0
IRL                 =  0

NQ(1)               =  5                        ! Number of Radial Quadrature Points
NQ(2)               =  10                        ! Number of Theta Quadrature Points
NQ(3)               =  3                        ! Number of Phi Quadrature Points

Left_Limit          = -0.50_idp
Right_Limit         = +0.50_idp

Verbose             = .TRUE.
!Verbose             = .FALSE.

Print_Results_Flag  = .TRUE.
!Print_Results_Flag  = .FALSE.

Print_Setup_Flag    = .TRUE.
!Print_Setup_Flag    = .FALSE.

!Print_Time_Flag     = .TRUE.
Print_Time_Flag     = .FALSE.

!Print_Cond_Flag     = .TRUE.
Print_Cond_Flag     = .FALSE.

Suffix_Input        = "Params"

SphereType          = "Oblate "
!SphereType          = "Prolate"

SemiMinor_Axis      = Radius_Sun/10
SemiMajor_Axis      = Radius_Sun*10

Rho                 = Density_Sun*5

CFA_Eqs = (/ 1, 1, 1, 1, 1 /)


Write_Results_R_Samps = 512
Write_Results_T_Samps = 256
Write_Results_P_Samps = 1



!############################################################!
!#                                                          #!
!#                 Initialize AMReX & MPI                   #!
!#                                                          #!
!############################################################!
CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nPROCS,ierr)

call amrex_init()
CALL amrex_amrcore_init()

CALL Init_AMReX_Parameters()
CALL Set_Units( Units_Input )

CALL Set_MLS_Parameters( SphereType,                    &
                         SemiMajor_Axis*Centimeter,     &
                         SemiMinor_Axis*Centimeter,     &
                         Rho*rho_Units )

Driver_xL(1) = Left_Limit
Driver_xL(2) = Right_Limit
Driver_NQ    = NQ
ALLOCATE(Driver_RQ_xLocs(Driver_NQ(1)))

!############################################################!
!#                                                          #!
!#                       Start Program                      #!
!#                                                          #!
!############################################################!
DO M_Index = M_Index_Min, M_Index_Max


    WRITE(Suffix_Tail,'(A)')Letter_Table_Upper(nLevels)

    ALLOCATE( Input_R_Quad(1:NQ(1)) )
    ALLOCATE( Input_T_Quad(1:NQ(2)) )
    ALLOCATE( Input_P_Quad(1:NQ(3)) )

    Input_R_Quad = Initialize_LG_Quadrature_Locations(NQ(1))
    Input_T_Quad = Initialize_LG_Quadrature_Locations(NQ(2))
    Input_P_Quad = Initialize_LG_Quadrature_Locations(NQ(3))

    Input_R_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_R_Quad)
    Input_T_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_T_Quad)
    Input_P_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_P_Quad)

    Driver_RQ_xlocs = Input_R_Quad


    
    !############################################################!
    !#                                                          #!
    !#                   Initialize Poseidon                    #!
    !#                                                          #!
    !############################################################!
    CALL Initialize_Poseidon &
       (    Source_NQ                           = NQ,                   &
            Source_xL                           = [Left_Limit, Right_Limit],    &
            Source_RQ_xlocs                     = Input_R_Quad,         &
            Source_TQ_xlocs                     = Input_T_Quad,         &
            Source_PQ_xlocs                     = Input_P_Quad,         &
            Source_Units                        = Units_Input,          &
            Source_Radial_Boundary_Units        = "cm",                 &
            Integration_NQ_Option               = NQ,                   &
            Eq_Flags_Option                     = CFA_Eqs,              &
            AMReX_FEM_Refinement_Option         = IFL,                  &
            AMReX_Integral_Refinement_Option    = IRL,                  &
            Verbose_Option                      = Verbose,              &
            WriteAll_Option                     = .FALSE.,              &
            Print_Setup_Option                  = Print_Setup_Flag,     &
            Write_Setup_Option                  = .TRUE.,              &
            Print_Results_Option                = Print_Results_Flag,   &
            Write_Results_Option                = .TRUE.,               &
            Print_Timetable_Option              = Print_Time_Flag,      &
            Write_Timetable_Option              = .TRUE.,               &
            Write_Sources_Option                = .TRUE.,              &
            Print_Condition_Option              = Print_Cond_Flag,      &
            Write_Condition_Option              = .FALSE.,              &
            Write_FP_Diagnostics_Option         = .TRUE.,               &
            Suffix_Flag_Option                  = Suffix_Input,         &
            Suffix_Tail_Option                  = Suffix_Tail           )



    !############################################################!
    !#                                                          #!
    !#               Create & Input Source Values               #!
    !#                                                          #!
    !############################################################!
    CALL Driver_InitSource( )
    CALL Driver_SetSource( )


    !############################################################!
    !#                                                          #!
    !#          Calculate and Set Boundary Conditions           #!
    !#                                                          #!
    !############################################################!
    CALL Driver_SetBC( )


    !############################################################!
    !#                                                          #!
    !#              Calculate and Set Initial Guess             #!
    !#                                                          #!
    !############################################################!
    CALL Driver_SetGuess( )
    

    !############################################################!
    !#                                                          #!
    !#                         Run Poseidon                     #!
    !#                                                          #!
    !############################################################!
!    Call Driver_ConFactor_Loopx( )
    CALL Poseidon_Run()

    CALL Print_MacLaurin_Error()
    CALL Output_2D_Results()
    !############################################################!
    !#                                                          #!
    !#                      Close Poseidon                      #!
    !#                                                          #!
    !############################################################!
    CALL Poseidon_Close()

    DEALLOCATE( Input_R_Quad )
    DEALLOCATE( Input_T_Quad )
    DEALLOCATE( Input_P_Quad )


END DO ! M_Index


!############################################################!
!#                                                          #!
!#                    Close AMReX & MPI                     #!
!#                                                          #!
!############################################################!
call amrex_finalize()
CALL MPI_Finalize(ierr)


END PROGRAM Driver_Main



