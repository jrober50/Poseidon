  !################################################################################!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
PROGRAM MacLaurin_Mapping                                                           !##!
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
                       Shift_Units,    &
                       Centimeter,     &
                       Second,         &
                       Gram,            &
                       Grav_Constant_G

USE Initialization_Poseidon, &
               ONLY :  Initialize_Poseidon


USE Variables_IO, &
               ONLY :  Write_Flags,                 &
                       Write_Results_R_Samps,       &
                       Write_Results_T_Samps,       &
                       File_Suffix

USE Variables_Derived, &
                ONLY : LM_Length

USE Variables_MPI, &
               ONLY :  ierr

USE Variables_Functions, &
                ONLY :  Potential_Solution

USE FP_Functions_Results , &
                ONLY : Calc_1D_CFA_Values_FP

USE Poseidon_IO_Parameters, &
               ONLY :  Poseidon_Results_Dir


USE Functions_Mesh, &
               ONLY :  Create_3D_Mesh

USE Functions_Quadrature, &
               ONLY :  Initialize_LG_Quadrature_Locations

USE Functions_Mapping, &
               ONLY :  Map_From_X_Space

USE Poseidon_IO_Module, &
               ONLY :  Open_Run_Report_File,        &
                       Open_New_File,               &
                       Output_Poseidon_Sources_3D,  &
                       OPEN_FILE_INQUISITION

USE IO_Write_Final_Results, &
            ONLY :  Write_Final_Results

USE IO_Print_Results, &
           ONLY :  Print_Results

USE Poseidon_Main_Module, &
               ONLY :  Poseidon_CFA_Set_Uniform_Boundary_Conditions,       &
                       Poseidon_Run,                                       &
                       Poseidon_Close

USE Initial_Guess_Module, &
                ONLY :  Poseidon_Input_Guess

USE Source_Input_Module, &
                ONLY :  Poseidon_Input_Sources,          &
                        Poseidon_Input_Sources_Serial,  &
                        Poseidon_Input_MCL_Sources

USE MacLaurin_Module, &
                ONLY :  Initialize_MacLaurin_Sources,   &
                        MacLaurin_Potential_Sub

USE SelfSimilar_Module, &
                ONLY :  Initialize_Yahil_Sources


USE MPI


IMPLICIT NONE

!                                       !
!   Poseidon Initialization Variables   !
!                                       !
INTEGER                                                 ::  Test_Number

INTEGER                                                 ::  Dimension_Input

INTEGER                                                 ::  Mesh_Type
INTEGER                                                 ::  Method_Type

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
LOGICAL                                                 ::  Flat_Guess

INTEGER,   DIMENSION(5)                                 ::  CFA_EQs


REAL(idp)                                               ::  Psi_BC
REAL(idp)                                               ::  AlphaPsi_BC
REAL(idp)                                               ::  Shift_Vector_BC
CHARACTER(LEN=1), DIMENSION(1:5)                        ::  INNER_BC_TYPES, OUTER_BC_TYPES
REAL(idp), DIMENSION(1:5)                               ::  INNER_BC_VALUES, OUTER_BC_VALUES

CHARACTER(LEN=10)                                       ::  Suffix_Input
CHARACTER(LEN=1)                                        ::  Suffix_Tail

REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Cur_R_Locs
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

INTEGER                                                 ::  re, rq, tq, pq, q, i

INTEGER                                                 ::  Guess_Type

INTEGER, DIMENSION(:), ALLOCATABLE                      ::  RE_Table
INTEGER, DIMENSION(:), ALLOCATABLE                      ::  TE_Table

INTEGER                                                 ::  Ecc_Index
INTEGER                                                 ::  Ecc_Index_Min
INTEGER                                                 ::  Ecc_Index_Max


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

INTEGER                                                 ::  Max_Iterations
REAL(idp)                                               ::  CC_Option

REAL(idp)                                               ::  Perturbation
REAL(idp)                                               ::  Offset

INTEGER                                                 ::  HCT_Fileid
CHARACTER(LEN = 100)                                    ::  HCT_Filename

INTEGER, DIMENSION(1:8)                                 ::  Anderson_M_Values
CHARACTER(LEN=1), DIMENSION(1:7)                        ::  Letter_Table
REAL(idp), DIMENSION(1:7)                               ::  Time_Values
INTEGER, DIMENSION(1:2)                                 ::  L_Values

REAL(idp), DIMENSION(1:7)                               ::  Ecc_Values
REAL(idp), DIMENSION(1:7)                               ::  Ecc_Ratios

INTEGER                                                 ::  nc_input

REAL(idp)                                               ::  Potential
CHARACTER(LEN=7)                                        ::  Spheroid_Type
REAL(idp)                                               ::  SemiMajor_Axis
REAL(idp)                                               ::  SemiMinor_Axis
REAL(idp)                                               ::  Rho_O

INTEGER                                                 ::  ui, lm


REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              :: Tmp_Lapse
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              :: Tmp_ConFac
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              :: Tmp_Shift

REAL(idp)                                               ::  Kappa
REAL(idp)                                               ::  Gamma

INTEGER                                                 ::  AMReX_Levels

CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

ALLOCATE( RE_Table(1:9) )
ALLOCATE( TE_Table(1:6) )

111 FORMAT (A,I4.4,/,A,I1.1,/,A,I2.2,/,A,I2.2,/,A,F6.3,A,/)
112 FORMAT (A,I4.4,/A,I4.4,/,A,I1.1,/,A,I2.2,/,A,I2.2,/,A,F6.3,A,/)




!############################################################!
!#                                                          #!
!#                      Test Parameters                     #!
!#                                                          #!
!############################################################!
Units_Input         = "G"
Method_Type         = 3     ! 1 = NR, 2 = FP, 3 = XCFC

RE_Table            = (/ 10, 156, 1024, 160, 240, 320, 400, 600, 4096 /)
TE_Table            = (/ 8, 30, 90, 120, 150, 180/)
Anderson_M_Values   = (/ 1, 2, 3, 4, 5, 10, 20, 50 /)
Ecc_Values          = (/ 0.00_idp, 0.250_idp, 0.50_idp, 0.750_idp, 0.950_idp, 0.990_idp, 1.4142E-4_idp /)
Ecc_Ratios          = (/ 0.99498743710_idp, 0.994987_idp, 0.968246_idp, 0.866025_idp,      &
                         0.661438_idp, 0.31225_idp, 0.141067_idp /)
L_Values            = (/ 5, 10 /)

AMReX_Levels        =  2

Ecc_Index_Min       =  1
Ecc_Index_Max       =  1

RE_Index_Min        =  1
RE_Index_Max        =  RE_Index_Min
NE(2)               =  TE_Table(1)

Degree_Min          =  1
Degree_Max          =  1

L_Limit_Min         =  2
L_Limit_Max         =  2

M_Index_Min         =  3
M_Index_Max         =  3

Guess_Type          =  1            !  1 = Flat, 2 = Educated, 3 = Perturbed Educated.
Perturbation        =  -0.01_idp    !  If Guess_Type == 3, rho is the perturbation parameter

Rho_O               =   1.0E+20_idp  ! cm/g^3
SemiMajor_Axis      =   250.0_idp   ! cm
SemiMinor_Axis      =   2500.0_idp
Domain_Edge(2)      =   1.0E+4_idp  ! cm


Spheroid_Type       =  'Oblate '
!Spheroid_Type       =  'Prolate'

Dimension_Input     = 3
Max_Iterations      = 10
CC_Option           = 1.0E-10_idp

Mesh_Type           = 6                         ! 1 = Uniform, 2 = Log, 3 = Split, 4 = Zoom, 5 = ML
Domain_Edge(1)      = 0.0_idp                   ! Inner Radius (cm)
                  ! Outer Radius (cm)
nc_input            = RE_Table(RE_Index_Min)/2.0

!SemiMajor_Axis      = Domain_Edge(2)/2.0_idp


NE(1)               = 128                       ! Number of Radial Elements
!NE(2)               = 120                       ! Number of Theta Elements
NE(3)               = 1                         ! Number of Phi Elements

NQ(1)               = 5                       ! Number of Radial Quadrature Points
NQ(2)               = 5                         ! Number of Theta Quadrature Points
NQ(3)               = 1                         ! Number of Phi Quadrature Points


Verbose             = .TRUE.
!Verbose             = .TRUE.
Suffix_Input        = "Params"

CFA_Eqs = (/ 1, 1, 1, 1, 1 /)


Letter_Table = (/ "A","B","C","D","E","F","G" /)


Write_Results_R_Samps = 512
Write_Results_T_Samps = 256

CALL Set_Units(Units_Input)



Domain_Edge    = Domain_Edge*Centimeter
SemiMajor_Axis = SemiMajor_Axis * Centimeter
SemiMinor_Axis = SemiMinor_Axis * Centimeter
Rho_O          = Rho_O * (Gram / Centimeter**3)



DO M_Index = M_Index_Min, M_Index_Max
DO Ecc_Index = Ecc_Index_Min, Ecc_Index_Max
DO RE_Index = RE_Index_Min, RE_Index_Max
DO Degree_Input = Degree_Min, Degree_Max
DO L_Limit_Input = L_Limit_Min, L_Limit_Max

    NE(1) = RE_Table(RE_Index)
    NQ(3) = 2*L_Limit_Input + 1
    


    Suffix_Tail = Letter_Table(AMReX_Levels)



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
                        nc = nc_input,           &
                        Zoom = 1.032034864238313_idp,   &
                        Core_Radius = SemiMinor_Axis,   &
                        SemiMajor = SemiMajor_Axis,     &
                        SemiMinor = SemiMinor_Axis,     &
                        Levels_In = AMReX_Levels       )


    !############################################################!
    !#                                                          #!
    !#                   Initialize Poseidon                    #!
    !#                                                          #!
    !############################################################!
    CALL Initialize_Poseidon &
       (   Dimensions_Option           = Dimension_Input,  &
           FEM_Degree_Option           = Degree_Input,     &
           L_Limit_Option              = L_Limit_Input,    &
           Units_Option                = Units_Input,      &
           Domain_Edge_Option          = Domain_Edge,      &
           NE_Option                   = NE,               &
           NQ_Option                   = NQ,               &
           r_Option                    = x_e,              &
           t_Option                    = y_e,              &
           p_Option                    = z_e,              &
    !        dr_Option                   = dx_c,             &
    !        dt_Option                   = dy_c,             &
    !        dp_Option                   = dz_c              &
           Suffix_Flag_Option          = Suffix_Input,    &
           Suffix_Tail_Option          = Suffix_Tail,     &
           Method_Flag_Option          = Method_Type,     &
           CFA_Eq_Flags_Option         = CFA_Eqs,         &
           Max_Iterations_Option       = Max_Iterations,  &
           Convergence_Criteria_Option = CC_Option,  &
           Anderson_M_Option           = Anderson_M_Values(M_Index),   &
           Verbose_Option              = Verbose,           &
           WriteAll_Option             = .FALSE.,           &
           Print_Setup_Option          = .TRUE.,            &
           Write_Setup_Option          = .FALSE.,           &
           Print_Results_Option        = .TRUE.,            &
           Write_Results_Option        = .TRUE.,           &
           Print_Timetable_Option      = .FALSE.,            &
           Write_Timetable_Option      = .FALSE.,                       &
           Write_Sources_Option        = .TRUE.                        )






    !############################################################!
    !#                                                          #!
    !#               Create & Input Source Values               #!
    !#                                                          #!
    !############################################################!

    CALL Initialize_MacLaurin_Sources(  Rho_O, Spheroid_Type,                       &
                                        SemiMajor_Axis, SemiMinor_Axis,             &
                                        NQ,                                         &
                                        Input_R_Quad, Input_T_Quad, Input_P_Quad,   &
                                        Left_Limit, Right_Limit,                    &
                                        NE,                                         &
                                        x_e, y_e, z_e,                              &
                                        Local_E, Local_S, Local_Si                  )




    CALL Poseidon_Input_Sources_Serial( Local_E, Local_S, Local_Si,                 &
                                        NE(1), NE(2), NE(3),                        &
                                        NQ(1), NQ(2), NQ(3),                        &
                                        Input_R_Quad, Input_T_Quad, Input_P_Quad,   &
                                        Left_Limit, Right_Limit                     )


    IF ( .TRUE. ) THEN
    CALL Output_Poseidon_Sources_3D( Local_E, Local_S, Local_Si,                &
                                     NE(1), NE(2), NE(3),                       &
                                     NQ(1), NQ(2), NQ(3),                       &
                                     Input_R_Quad, Input_T_Quad, Input_P_Quad,  &
                                     Left_Limit, Right_Limit                    )

    END IF


    !############################################################!
    !#                                                          #!
    !#          Calculate and Set Boundary Conditions           #!
    !#                                                          #!
    !############################################################!
    
    IF ( SemiMajor_Axis == SemiMinor_Axis ) THEN
        
        Potential = -((4.0_idp/3.0_idp)*pi*SemiMajor_Axis**3*Rho_O)*Grav_Constant_G/Domain_Edge(2)        &
                  * gram / centimeter
    ELSE
        CALL MacLaurin_Potential_Sub(Domain_Edge(2), 0.5_idp*pi, 0.0_idp,Potential)
    END IF

    Psi_BC      = 1.0_idp - 0.5_idp*Potential/C_Square
    AlphaPsi_BC = 1.0_idp + 0.5_idp*Potential/C_Square



!    Psi_BC = 1.0000000015550770
!    AlphaPSi_BC = 0.99999999844492315

    INNER_BC_TYPES = (/"N", "N","N","N","N"/)
    OUTER_BC_TYPES = (/"D", "D","D","D","D"/)


    INNER_BC_VALUES = (/0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp /)
    OUTER_BC_VALUES = (/Psi_BC,  AlphaPsi_BC, 0.0_idp, 0.0_idp, 0.0_idp /)
!    OUTER_BC_VALUES = (/Psi_BC,  AlphaPsi_BC, Shift_Vector_BC, 0.0_idp, 0.0_idp /)

 
    CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("I", INNER_BC_TYPES, INNER_BC_VALUES)
    CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("O", OUTER_BC_TYPES, OUTER_BC_VALUES)









    !############################################################!
    !#                                                          #!
    !#              Calculate and Set Initial Guess             #!
    !#                                                          #!
    !############################################################!

    IF ( Guess_Type == 1 ) THEN

        Psi_Guess = 1.0_idp
        AlphaPsi_Guess = 1.0_idp
        Beta_Guess = 0.0_idp
   
    ELSE IF ( Guess_Type == 2 ) THEN
        !     Lapse function Coefficients in 1D correspond to the value of the function
        !     at the location of the FEM nodes.

        Beta_Guess = 0.0_idp

        DO re = 1,NE(1)
        Cur_R_Locs(:) = dx_c(re)*(Input_R_Quad(:) + Left_Limit)+  x_e(re)
        DO rq = 1,NQ(1)
            CALL MacLaurin_Potential_Sub(cur_R_Locs(rq), 0.5_idp*pi, 0.0_idp,Potential)
            Psi_Guess(rq, re-1, 0, 0) = 1.0_idp    &
                                    - 0.5_idp * Potential / C_Square

            AlphaPsi_Guess(rq, re-1, 0, 0) = 1.0_idp    &
                                    + 0.5_idp*Potential / C_Square

        END DO ! rq
        END DO ! re


    ELSE IF ( Guess_Type == 3 ) THEN
        !     Lapse function Coefficients in 1D correspond to the value of the function
        !     at the location of the FEM nodes.  These have been perturbed.


         Beta_Guess = 0.0_idp


         DO re = 1,NE(1)
         Cur_R_Locs(:) = dx_c(re)*(Input_R_Quad(:) + Left_Limit)+  x_e(re)
         DO rq = 1,NQ(1)

             Psi_Guess(rq, re-1, 0, 0) = 1.0_idp    &
                                     - 0.5_idp*Potential_Solution(Cur_R_Locs(rq), 0.0_idp, 0.0_idp)/C_Square

             AlphaPsi_Guess(rq, re-1, 0, 0) = 1.0_idp    &
                                     + 0.5_idp*Potential_Solution(Cur_R_Locs(rq), 0.0_idp, 0.0_idp)/C_Square


            Offset = (Cur_R_locs(rq)/Domain_Edge(2)) * (1.0_idp + Perturbation)
            Psi_Guess(rq, re-1, 0, 0) = Psi_Guess(rq, re-1, 0, 0)*Offset
            AlphaPsi_Guess(rq, re-1, 0, 0) = AlphaPsi_Guess(rq, re-1, 0, 0)*Offset

         END DO ! rq
         END DO ! re


    END IF


    CALL Poseidon_Input_Guess(  Psi_Guess,                                  &
                                AlphaPsi_Guess,                             &
                                Beta_Guess,                                 &
                                NE(1), NE(2), NE(3),                        &
                                NQ(1), NQ(2), NQ(3),                        &
                                Input_R_Quad, Input_T_Quad, Input_P_Quad,   &
                                Left_Limit, Right_Limit                     )


    !############################################################!
    !#                                                          #!
    !#                         Run Poseidon                     #!
    !#                                                          #!
    !############################################################!
    Call Poseidon_Run()



    !############################################################!
    !#                                                          #!
    !#                       Output Results                     #!
    !#                                                          #!
    !############################################################!

    IF ( Verbose ) THEN
        PRINT*,"Final Results "
        CALL Print_Results()
        CALL Write_Final_Results()
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

    DEALLOCATE( Local_E, Local_S, Local_Si )

    DEALLOCATE( Psi_Guess )
    DEALLOCATE( AlphaPsi_Guess )
    DEALLOCATE( Beta_Guess )



END DO ! L_Limit
END DO ! Degree_Index
END DO ! RE_Index
END DO ! T_Index
END DO ! M_Index



CONTAINS






END PROGRAM MacLaurin_Mapping



