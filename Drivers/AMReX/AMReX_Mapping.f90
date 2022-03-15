  !################################################################################!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
PROGRAM AMReX_Mapping                                                               !##!
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

USE Poseidon_Units_Module, &
               ONLY :  Set_Units,      &
                       Centimeter

USE Initialization_AMReX, &
               ONLY :  Initialize_Poseidon_with_AMReX


USE Variables_IO, &
               ONLY :  Write_Results_R_Samps,      &
                       Write_Results_T_Samps

USE Variables_MPI, &
               ONLY :  ierr

USE Functions_Mesh, &
               ONLY :  Create_3D_Mesh

USE Functions_Quadrature, &
               ONLY :  Initialize_LG_Quadrature_Locations

USE Functions_Mapping, &
               ONLY :  Map_From_X_Space

USE Poseidon_IO_Module, &
               ONLY :  Open_Run_Report_File,       &
                       Output_Poseidon_Sources_3D

USE IO_Print_Results, &
               ONLY :  Print_Results

USE Poseidon_Main_Module, &
               ONLY :  Poseidon_Run,                                       &
                       Poseidon_Close

USE Driver_SetSource_Module, &
                ONLY:  Driver_SetSource

USE Driver_SetBC_Module, &
                ONLY:  Driver_SetBC

USE Driver_SetGuess_Module, &
                ONLY:  Driver_SetGuess

USE FP_IO_Module, &
               ONLY :  Output_FP_Timetable

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

USE Variables_AMReX_Core, &
            ONLY :  MF_Source,          &
                    AMReX_Num_Levels

USE Variables_FP, &
            ONLY :  FP_Coeff_Vector_A,      &
                    FP_Coeff_Vector_B

USE Variables_Mesh, &
ONLY :  Num_R_Elements,         &
        Num_T_Elements,         &
        Num_P_Elements,         &
        rlocs,                  &
        drlocs

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

INTEGER,   DIMENSION(5)                                 ::  CFA_EQs


CHARACTER(LEN=10)                                       ::  Suffix_Input
CHARACTER(LEN=1)                                        ::  Suffix_Tail

REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Cur_R_Locs
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_R_Quad
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_T_Quad
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_P_Quad
REAL(idp)                                               ::  Left_Limit
REAL(idp)                                               ::  Right_Limit


INTEGER                                                 ::  myID, nPROCS

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
REAL(idp), DIMENSION(3)                                 ::  Yahil_Params


INTEGER                                                 ::  Max_Iterations
REAL(idp)                                               ::  CC_Option

INTEGER                                                 ::  IRL
INTEGER                                                 ::  IFL

INTEGER, DIMENSION(1:8)                                 ::  Anderson_M_Values
CHARACTER(LEN=1), DIMENSION(1:10)                       ::  Letter_Table
REAL(idp), DIMENSION(1:6)                               ::  Time_Values
INTEGER, DIMENSION(1:2)                                 ::  L_Values

INTEGER                                                 ::  i
INTEGER                                                 ::  Level


CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nPROCS,ierr)


ALLOCATE( RE_Table(1:9) )


!############################################################!
!#                                                          #!
!#                      Test Parameters                     #!
!#                                                          #!
!############################################################!
Units_Input         = "G"
Solver_Type         = 3

RE_Table            = (/ 32, 128, 160, 240, 320, 400, 600, 256, 512 /)
Anderson_M_Values   = (/ 1, 2, 3, 4, 5, 10, 20, 50 /)
Time_Values         = (/ 51.0_idp, 15.0_idp, 5.0_idp, 1.50_idp, 0.5_idp, 0.05_idp /)
L_Values            = (/ 5, 10 /)

T_Index_Min         =  5
T_Index_Max         =  5

M_Index_Min         =  3
M_Index_Max         =  3

RE_Index_Min        =  2
RE_Index_Max        =  2

Degree_Min          =  1
Degree_Max          =  1

L_Limit_Min         =  0
L_Limit_Max         =  0

IFL                 =  0
IRL                 =  0

Guess_Type          =  1            !  1 = Flat, 2 = Educated, 3 = Perturbed Educated.

Kappa               = 953946015514834.4
Gamma               = 1.30_idp

!Suffix_Tail         = "A"
!Convergence_Criteria = 1.0E-8_idp

Dimension_Input     = 3

Max_Iterations      = 10
CC_Option           = 1.0E-10_idp

Mesh_Type           = 1                         ! 1 = Uniform, 2 = Log, 3 = Split, 4 = Zoom
Domain_Edge(1)      = 0.0_idp                   ! Inner Radius (cm)
Domain_Edge(2)      = 1E9_idp                  ! Outer Radius (cm)


NE(1)               = 128                      ! Number of Radial Elements
NE(2)               = 1                        ! Number of Theta Elements
NE(3)               = 1                        ! Number of Phi Elements

NQ(1)               = 5                        ! Number of Radial Quadrature Points
NQ(2)               = 1                        ! Number of Theta Quadrature Points
NQ(3)               = 1                        ! Number of Phi Quadrature Points


!Verbose             = .TRUE.
Verbose             = .FALSE.
Suffix_Input        = "Params"



CFA_Eqs = (/ 1, 1, 1, 1, 1 /)


Letter_Table = (/ "A","B","C","D","E","F","G","H","I","J" /)


Write_Results_R_Samps = 256
Write_Results_T_Samps = 1

CALL Set_Units(Units_Input)

call amrex_init()
CALL amrex_amrcore_init()

CALL Init_AMReX_Parameters()

!Domain_Edge = Domain_Edge*Centimeter

Domain_Edge(1) = xL(1)*Centimeter
Domain_Edge(2) = xR(1)*Centimeter

DO M_Index = M_Index_Min, M_Index_Max
DO T_Index = T_Index_Min, T_Index_Max
DO RE_Index = RE_Index_Min, RE_Index_Max
DO Degree_Input = Degree_Min, Degree_Max
DO L_Limit_Input = L_Limit_Min, L_Limit_Max

!    NE(1) = RE_Table(RE_Index)
    NE = nCells
    NQ(3) = 2*L_Limit_Input + 1

    Suffix_Tail = Letter_Table(nLevels)



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



    Input_R_Quad = Initialize_LG_Quadrature_Locations(NQ(1))
    Input_T_Quad = Initialize_LG_Quadrature_Locations(NQ(2))
    Input_P_Quad = Initialize_LG_Quadrature_Locations(NQ(3))

    Left_Limit  = -0.50_idp
    Right_Limit = +0.50_idp

    Input_R_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_R_Quad)
    Input_T_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_T_Quad)
    Input_P_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_P_Quad)




    CALL Open_Run_Report_File()

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
    CALL Initialize_Poseidon_with_AMReX &
       (    FEM_Degree_Option                   = Degree_Input,         &
            L_Limit_Option                      = L_Limit_Input,        &
            Units_Option                        = Units_Input,          &
            Domain_Edge_Option                  = Domain_Edge,          &
            Coarse_NE_Option                    = NE,                   &
            NQ_Option                           = NQ,                   &
            Max_Iterations_Option               = Max_Iterations,       &
            Convergence_Criteria_Option         = CC_Option,            &
            Anderson_M_Option                   = Anderson_M_Values(M_Index),   &
            CFA_Eq_Flags_Option                 = CFA_Eqs,              &
            AMReX_Max_Levels_Option             = nlevels-1,            &
            AMReX_Max_Grid_Size_Option          = MaxGridSizeX,         &
            AMReX_FEM_Refinement_Option         = IFL,                  &
            AMReX_Integral_Refinement_Option    = IRL,                  &
            Poisson_Mode_Option                 = .FALSE.,              &
            Verbose_Option                      = Verbose,              &
            WriteAll_Option                     = .FALSE.,              &
            Print_Setup_Option                  = .TRUE.,               &
            Write_Setup_Option                  = .FALSE.,              &
            Print_Results_Option                = .TRUE.,               &
            Write_Results_Option                = .TRUE.,               &
            Print_Timetable_Option              = .TRUE.,               &
            Write_Timetable_Option              = .TRUE.,               &
            Write_Sources_Option                = .FALSE.,              &
            Suffix_Flag_Option                   = Suffix_Input,        &
            Suffix_Tail_Option                   = Suffix_Tail          )



    !############################################################!
    !#                                                          #!
    !#               Create & Input Source Values               #!
    !#                                                          #!
    !############################################################!
    Yahil_Params = [Time_Values(T_Index), Kappa, Gamma]
    CALL Driver_SetSource(  NQ, Yahil_Params, nLevels )


    !############################################################!
    !#                                                          #!
    !#          Calculate and Set Boundary Conditions           #!
    !#                                                          #!
    !############################################################!
!    PRINT*,"Before Driver_SetBC"
    CALL Driver_SetBC( x_e(NE(1)) )



    !############################################################!
    !#                                                          #!
    !#              Calculate and Set Initial Guess             #!
    !#                                                          #!
    !############################################################!
!    PRINT*,"Before Driver_SetGuess"

    ! These values are established during source input.
    ! As the original NE accounts for only the coarsest level,
    ! we need to do this to get the total number of leaf elements
    NE(1) = Num_R_Elements
    NE(2) = Num_T_Elements
    NE(3) = Num_P_Elements


    CALL Driver_SetGuess(   NE, NQ,         &
                            drlocs, rlocs,      &
                            Input_R_Quad,   &
                            Input_T_Quad,   &
                            Input_P_Quad,   &
                            Left_Limit,     &
                            Right_Limit,    &
                            Guess_Type      )
    
    !############################################################!
    !#                                                          #!
    !#                         Run Poseidon                     #!
    !#                                                          #!
    !############################################################!
!    PRINT*,"Before Poseidon_Run"

    Call Poseidon_Run()

!    PRINT*,FP_Coeff_Vector_A

    !############################################################!
    !#                                                          #!
    !#                       Output Results                     #!
    !#                                                          #!
    !############################################################!
    IF (Verbose .EQV. .TRUE. ) THEN
        CALL MPI_Barrier(Poseidon_Comm_World, ierr )
        DO i = 0,nPROCS_Poseidon-1
            IF (myID_Poseidon == i ) THEN
                PRINT*,"Final Results ",myID_Poseidon
                CALL Print_Results()
            END IF
            CALL MPI_Barrier(Poseidon_Comm_World, ierr)
        END DO
    END IF


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




END DO ! L_Limit
END DO ! Degree_Index
END DO ! RE_Index
END DO ! T_Index
END DO ! M_Index



call amrex_finalize()
CALL MPI_Finalize(ierr)


END PROGRAM AMReX_Mapping


