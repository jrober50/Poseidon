   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
PROGRAM HCT_Mapping                                                                 !##!
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
            ONLY :  C_Square,       &
                    Set_Units,      &
                    Centimeter,     &
                    Gram

USE Initialization_Poseidon, &
            ONLY :  Initialize_Poseidon

USE Variables_IO, &
            ONLY :  Write_Results_R_Samps,      &
                    Write_Results_T_Samps,      &
                    File_Suffix

USE Variables_MPI, &
            ONLY :  ierr

USE Functions_Mesh, &
            ONLY :  Create_3D_Mesh

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature_Locations

USE Maps_X_Space, &
            ONLY :  Map_From_X_Space

USE IO_Print_Results, &
            ONLY :  Print_Results

USE IO_Write_Final_Results, &
            ONLY :  Write_Final_Results

USE IO_Convergence_Output, &
            ONLY :  Output_Convergence_Reports

USE Poseidon_Main_Module, &
            ONLY :  Poseidon_Run,                                       &
                    Poseidon_Close

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

USE IO_Write_Final_Results, &
            ONLY :  Write_Final_Results


USE MPI


IMPLICIT NONE

!                                       !
!   Poseidon Initialization Variables   !
!                                       !
INTEGER                                                 ::  Test_Number

INTEGER                                                 ::  Dimension_Input

INTEGER                                                 ::  Mesh_Type
INTEGER                                                 ::  Solver_Type

INTEGER                                                 ::  Max_Iterations

CHARACTER(LEN = 1)                                      ::  Units_Input

INTEGER, DIMENSION(3)                                   ::  NE
INTEGER, DIMENSION(3)                                   ::  NQ
REAL(idp), DIMENSION(2)                                 ::  Domain_Edge
REAL(idp)                                               ::  Star_Surface
INTEGER                                                 ::  Num_DOF


REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  dx_c, dy_c, dz_c
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  x_e, y_e, z_e
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  x_c, y_c, z_c

LOGICAL                                                 ::  Verbose
LOGICAL                                                 ::  Print_Results_Flag
LOGICAL                                                 ::  Flat_Guess

INTEGER,   DIMENSION(5)                                 ::  CFA_EQs



REAL(idp)                                               ::  Shift_Vector_BC
CHARACTER(LEN=1), DIMENSION(1:5)                        ::  INNER_BC_TYPES, OUTER_BC_TYPES
REAL(idp), DIMENSION(1:5)                               ::  INNER_BC_VALUES, OUTER_BC_VALUES

CHARACTER(LEN=10)                                       ::  Suffix_Input
CHARACTER(LEN=1)                                        ::  Suffix_Tail

REAL(idp)                                               ::  Alpha
REAL(idp)                                               ::  Star_Radius


REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_R_Quad
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_T_Quad
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_P_Quad
REAL(idp)                                               ::  Left_Limit
REAL(idp)                                               ::  Right_Limit

REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Local_E
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Local_S
REAL(idp), DIMENSION(:,:,:,:,:), ALLOCATABLE            ::  Local_Si

REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Psi_Guess
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  AlphaPsi_Guess
REAL(idp), DIMENSION(:,:,:,:,:), ALLOCATABLE            ::  Beta_Guess

INTEGER                                                 ::  myID

INTEGER                                                 ::  Iter_Max

INTEGER                                                 ::  re, rq, i

INTEGER                                                 ::  Guess_Type
INTEGER                                                 ::  AMReX_Levels
REAL(idp)                                               ::  CC_Option

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

INTEGER                                                 ::  Tol_Index
INTEGER                                                 ::  Tolerance_Index_Min
INTEGER                                                 ::  Tolerance_Index_Max

INTEGER                                                 ::  External_Iter
INTEGER                                                 ::  External_Iter_Max


REAL(idp)                                               ::  Perturbation

INTEGER, DIMENSION(1:8)                                 ::  Anderson_M_Values
REAL(idp), DIMENSION(1:7)                               ::  Tolerance_Values
CHARACTER(LEN=1), DIMENSION(1:7)                        ::  Tolerance_Letters


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

RE_Table          = (/ 2, 4, 6, 8, 16, 32, 512, 1024 /)
Anderson_M_Values   = (/ 1, 2, 3, 4, 5, 10, 20, 50 /)
Tolerance_Values  = (/ 1.0E-6, 1.0E-8, 1.0E-10, 1.0E-12, 1.0E-14, 1.0E-15, 1.0E-30 /)

Tolerance_Index_Min =  6
Tolerance_Index_Max =  6

External_Iter_Max   =  500

M_Index_Min         =  3
M_Index_Max         =  3

RE_Index_Min        =  8
RE_Index_Max        =  8

Degree_Min          =  1
Degree_Max          =  1

L_Limit_Min         =  0
L_Limit_Max         =  0

Guess_Type          =  1            !  1 = Flat, 2 = Educated, 3 = Perturbed Educated.
Perturbation        =  -0.01_idp    !  If Guess_Type == 3, rho is the perturbation parameter

Alpha               =  sqrt(5.0_idp)
Star_Radius         =  1.0E+5_idp               ! (cm)

!Suffix_Tail         = "A"
!Convergence_Criteria = 1.0E-8_idp

Dimension_Input     = 3

Max_Iterations      = 100
CC_Option           = 1.0E-12_idp


Mesh_Type           = 1
Domain_Edge(1)      = 0.0_idp                   ! Inner Radius (cm)
Domain_Edge(2)      = 8.0_idp * Star_Radius     ! Outer Radius (cm)


NE(1)               = 128                       ! Number of Radial Elements
NE(2)               = 1                         ! Number of Theta Elements
NE(3)               = 1                         ! Number of Phi Elements

NQ(1)               = 10                        ! Number of Radial Quadrature Points
NQ(2)               = 1                         ! Number of Theta Quadrature Points
NQ(3)               = 1                         ! Number of Phi Quadrature Points


Verbose             = .TRUE.
!Verbose             = .FALSE.
!Print_Results_Flag  = .TRUE.
Print_Results_Flag  = .FALSE.
Suffix_Input        = "Params"

!CFA_Eqs = (/ 1, 0, 0, 0, 0 /)
CFA_Eqs = (/ 1, 1, 1, 1, 1 /)

Tolerance_Letters = (/ "A","B","C","D","E","F","G" /)


Write_Results_R_Samps = 256
Write_Results_T_Samps = 1

CALL Set_Units(Units_Input)



Star_Radius = Star_Radius*Centimeter
Domain_Edge = Domain_Edge*Centimeter


DO Tol_Index = Tolerance_Index_Min, Tolerance_Index_Max
!
!if ( Tol_Index == 6 ) then
!    cycle
!end if

DO M_Index = M_Index_Min, M_Index_Max
DO RE_Index = RE_Index_Min, RE_Index_Max
DO Degree_Input = Degree_Min, Degree_Max
DO L_Limit_Input = L_Limit_Min, L_Limit_Max

    NE(1) = RE_Table(RE_Index)
    NQ(3) = 2*L_Limit_Input + 1
    Suffix_Tail = Tolerance_Letters(Tol_Index)

    

    Num_DOF = NQ(1)*NQ(2)*NQ(3)

    ALLOCATE( Input_R_Quad(1:NQ(1)) )
    ALLOCATE( Input_T_Quad(1:NQ(2)) )
    ALLOCATE( Input_P_Quad(1:NQ(3)) )

    ALLOCATE( x_e(0:NE(1)), y_e(0:NE(2)), z_e(0:NE(3)) )
    ALLOCATE( x_c(1:NE(1)), y_c(1:NE(2)), z_c(1:NE(3)) )
    ALLOCATE( dx_c(1:NE(1)) )
    ALLOCATE( dy_c(1:NE(2)) )
    ALLOCATE( dz_c(1:NE(3)) )

    ALLOCATE( Local_E(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )       )
    ALLOCATE( Local_S(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )       )
    ALLOCATE( Local_Si(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1, 1:3)  )

    ALLOCATE( Psi_Guess(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 ) )
    ALLOCATE( AlphaPsi_Guess(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 ) )
    ALLOCATE( Beta_Guess(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1,1:3 ) )

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
           Print_Setup_Option           = .TRUE.,                       &
           Write_Setup_Option           = .TRUE.,                       &
           Print_Results_Option         = Print_Results_Flag,           &
           Write_Results_Option         = .TRUE.,                       &
           Print_Timetable_Option       = .FALSE.,                      &
           Write_Timetable_Option       = .TRUE.,                       &
           Write_Sources_Option         = .TRUE.,                       &
           Print_Condition_Option       = .TRUE.,                       &
           Write_Condition_Option       = .TRUE.,                       &
           Suffix_Flag_Option           = Suffix_Input,                 &
           Suffix_Tail_Option           = Suffix_Tail                   )






        !############################################################!
        !#                                                          #!
        !#              Calculate and Set Initial Guess             #!
        !#                                                          #!
        !############################################################!
        CALL TimerStart( Timer_Driver_SetGuess )

        CALL Driver_SetGuess( )

        CALL TimerStop( Timer_Driver_SetGuess )


    !############################################################!
    !#                                                          #!
    !#               Create & Input Source Values               #!
    !#                                                          #!
    !############################################################!
    CALL TimerStart( Timer_Driver_SetSource )

    CALL Driver_SetSource(  NE, NQ,             &
                            dx_c, x_e, y_e,     &
                            Input_R_Quad,       &
                            Input_T_Quad,       &
                            Input_P_Quad,       &
                            Left_Limit,         &
                            Right_Limit,        &
                            myID,               &
                            Alpha,              &
                            Star_Radius         )

    CALL TimerStop( Timer_Driver_SetSource )







    !############################################################!
    !#                                                          #!
    !#          Calculate and Set Boundary Conditions           #!
    !#                                                          #!
    !############################################################!
    CALL TimerStart( Timer_Driver_SetBC )

    CALL Driver_SetBC( Domain_Edge(2), Alpha, Star_Radius )

    CALL TimerStop( Timer_Driver_SetBC )










    !############################################################!
    !#                                                          #!
    !#                         Run Poseidon                     #!
    !#                                                          #!
    !############################################################!
    Call Poseidon_Run()


    

    !############################################################!
    !#                                                          #!
    !#                          Loop Solve                      #!
    !#                                                          #!
    !############################################################!
    CALL Driver_ConFactor_Loop( NE, NQ,                    &
                                dx_c, x_e, y_e,            &
                                Input_R_Quad,               &
                                Input_T_Quad,               &
                                Input_P_Quad,               &
                                Left_Limit, Right_Limit,    &
                                myID,                       &
                                Alpha,                      &
                                Star_Radius,                &
                                Tolerance_Values(5),        &
                                External_Iter_Max           )


    !############################################################!
    !#                                                          #!
    !#                       Output Results                     #!
    !#                                                          #!
    !############################################################!

    IF ( Verbose ) THEN
        CALL Print_Results()
    END IF


    !Write_Results_Flag = 1
    IF ( .TRUE. ) THEN
!        CALL Output_Convergence_Reports()
        CALL Write_Final_Results( Output_Locations_Flag = 2)
    END IF






    !############################################################!
    !#                                                          #!
    !#                      Close Poseidon                      #!
    !#                                                          #!
    !############################################################!

    CALL Poseidon_Close()

    DEALLOCATE( Input_R_Quad )
    DEALLOCATE( Input_T_Quad )
    DEALLOCATE( Input_P_Quad )
    DEALLOCATE( x_e, y_e, z_e )
    DEALLOCATE( x_c, y_c, z_c )
    DEALLOCATE( dx_c, dy_c, dz_c )

    DEALLOCATE( Local_E, Local_S, Local_Si )

    DEALLOCATE( Psi_Guess )
    DEALLOCATE( AlphaPsi_Guess )
    DEALLOCATE( Beta_Guess )




END DO ! L_Limit
END DO ! Degree_Index
END DO ! RE_Index
END DO ! Tol_Index
END DO ! M_Index


CONTAINS

!############################################################!
!#                                                          #!
!#                   HCT_Solution Function                  #!
!#                                                          #!
!############################################################!
REAL FUNCTION HCT_Solution( r, Alpha, Beta, C, Star_Radius )

REAL(idp), INTENT(IN)                      ::   r
REAL(idp), INTENT(IN)                      ::   Alpha
REAL(idp), INTENT(IN)                      ::   Beta
REAL(idp), INTENT(IN)                      ::   C
REAL(idp), INTENT(IN)                      ::   Star_Radius


IF ( r .LE. Star_Radius ) THEN

    HCT_Solution = C*SQRT( (Alpha*Star_Radius)/(r*r + (Alpha*Star_Radius)**2 ) )

ELSE

    HCT_Solution = Beta/r + 1.0_idp

END IF

 

END FUNCTION HCT_Solution




!############################################################!
!#                                                          #!
!#                   HCT_Solution Function                  #!
!#                                                          #!
!############################################################!
REAL FUNCTION HCT_Perturbed_Solution( r, Alpha, Beta, C, Star_Radius, rho, Domain_Edge )

REAL(idp), INTENT(IN)                       ::  r
REAL(idp), INTENT(IN)                       ::  Alpha
REAL(idp), INTENT(IN)                       ::  Beta
REAL(idp), INTENT(IN)                       ::  C
REAL(idp), INTENT(IN)                       ::  Star_Radius
REAL(idp), INTENT(IN)                       ::  rho
REAL(idp), INTENT(IN)                       ::  Domain_Edge

REAL(idp)                                   ::  Perturbation
REAL(idp)                                   ::  Sol

Perturbation = (1.0_idp - r/Domain_Edge) * rho

IF ( r .LE. Star_Radius ) THEN

    Sol = C*SQRT( (Alpha*Star_Radius)/(r*r + (Alpha*Star_Radius)**2 ) )

ELSE

    Sol = Beta/r + 1.0_idp

END IF

HCT_Perturbed_Solution = Sol * ( 1.0_idp + Perturbation )
 

END FUNCTION HCT_Perturbed_Solution








END PROGRAM HCT_Mapping


