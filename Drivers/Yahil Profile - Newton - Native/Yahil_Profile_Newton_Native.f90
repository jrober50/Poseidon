  !################################################################################!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
PROGRAM Yahil_Profile_Poisson_Native                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  !################################################################################!


USE Poseidon_Kinds_Module, &
               ONLY :  idp

USE Poseidon_Letters_Module, &
            ONLY :  Letter_Table_Upper

USE Poseidon_Interface_Initialization, &
            ONLY :  Initialize_Poseidon

USE Poseidon_Interface_Run, &
            ONLY :  Poseidon_Run

USE Poseidon_Interface_Close, &
            ONLY :  Poseidon_Close

USE Variables_MPI, &
               ONLY :  ierr

USE Functions_Mesh, &
               ONLY :  Create_3D_Mesh

USE Functions_Quadrature, &
               ONLY :  Initialize_LG_Quadrature_Locations

USE Maps_X_Space, &
               ONLY :  Map_From_X_Space

USE Driver_SetSource_Module, &
                ONLY:  Driver_SetSource

USE Driver_SetBC_Module, &
                ONLY:  Driver_SetBC


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
CHARACTER(LEN=4)                                        ::  Suffix_Tail

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

INTEGER, DIMENSION(1:8)                                 ::  Anderson_M_Values
CHARACTER(LEN=1), DIMENSION(1:10)                       ::  Letter_Table
REAL(idp), DIMENSION(1:6)                               ::  Time_Values
INTEGER, DIMENSION(1:2)                                 ::  L_Values

INTEGER                                                 ::  i

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
Time_Values         = (/ 151.0_idp, 51.0_idp, 15.0_idp, 5.0_idp, 1.50_idp, 0.15_idp /)
L_Values            = (/ 5, 10 /)

T_Index_Min         =  1
T_Index_Max         =  1

M_Index_Min         =  3
M_Index_Max         =  3

RE_Index_Min        =  9
RE_Index_Max        =  9

Degree_Min          =  1
Degree_Max          =  1

L_Limit_Min         =  0
L_Limit_Max         =  0

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


NE(1)               = 128                       ! Number of Radial Elements
NE(2)               = 1                        ! Number of Theta Elements
NE(3)               = 1                         ! Number of Phi Elements

NQ(1)               = 10                        ! Number of Radial Quadrature Points
NQ(2)               = 1                         ! Number of Theta Quadrature Points
NQ(3)               = 1                         ! Number of Phi Quadrature Points


Verbose             = .TRUE.
!Verbose             = .FALSE.
Suffix_Input        = "Params"

CFA_Eqs = (/ 1, 1, 1, 1, 1 /)


DO M_Index = M_Index_Min, M_Index_Max
DO T_Index = T_Index_Min, T_Index_Max
DO RE_Index = RE_Index_Min, RE_Index_Max
DO Degree_Input = Degree_Min, Degree_Max
DO L_Limit_Input = L_Limit_Min, L_Limit_Max

    NE(1) = RE_Table(RE_Index)
    NQ(3) = 2*L_Limit_Input + 1

!    Suffix_Tail = Letter_Table(Solver_Type)
    WRITE(Suffix_Tail,'(A)')TRIM(Letter_Table(Solver_Type))


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
        (   Newtonian_Mode_Option        = .TRUE.,                       &
            Dimensions_Option            = Dimension_Input,              &
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
            Eq_Flags_Option              = CFA_Eqs,                      &
            Max_Iterations_Option        = Max_Iterations,               &
            Convergence_Criteria_Option  = CC_Option,                    &
            Anderson_M_Option            = Anderson_M_Values(M_Index),   &
            Verbose_Option               = Verbose,                      &
            WriteAll_Option              = .FALSE.,                      &
            Print_Setup_Option           = .TRUE.,                       &
            Write_Setup_Option           = .TRUE.,                       &
            Print_Results_Option         = .TRUE.,                       &
            Write_Results_Option         = .FALSE.,                       &
            Print_Timetable_Option       = .TRUE.,                       &
            Write_Timetable_Option       = .TRUE.,                       &
            Write_Sources_Option         = .TRUE.,                       &
            Print_Condition_Option       = .FALSE.,                       &
            Write_Condition_Option       = .TRUE.,                       &
            Suffix_Flag_Option           = Suffix_Input,                 &
            Suffix_Tail_Option           = Suffix_Tail                   )



    !############################################################!
    !#                                                          #!
    !#               Create & Input Source Values               #!
    !#                                                          #!
    !############################################################!
    Yahil_Params = [Time_Values(T_Index), Kappa, Gamma]
    CALL Driver_SetSource(  NE, NQ,                 &
                            dx_c, x_e,              &
                            Input_R_Quad,           &
                            Input_T_Quad,           &
                            Input_P_Quad,           &
                            Left_Limit,             &
                            Right_Limit,            &
                            Yahil_Params            )




    !############################################################!
    !#                                                          #!
    !#          Calculate and Set Boundary Conditions           #!
    !#                                                          #!
    !############################################################!
    CALL Driver_SetBC( )






    
    !############################################################!
    !#                                                          #!
    !#                         Run Poseidon                     #!
    !#                                                          #!
    !############################################################!
    Call Poseidon_Run()






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




CALL MPI_Finalize(ierr)





END PROGRAM Yahil_Profile_Poisson_Native



