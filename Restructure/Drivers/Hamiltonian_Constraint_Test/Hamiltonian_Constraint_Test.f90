   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
PROGRAM Hamiltonian_Constraint_Test                                                 !##!
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
                ONLY :  C_Square,       &
                        Set_Units,      &
                        Centimeter,     &
                        Gram

USE Initialization_Poseidon, &
                ONLY :  Initialize_Poseidon

USE Poseidon_Parameters, &
                ONLY :  Convergence_Criteria

USE Variables_IO, &
                ONLY :  Write_Results_Flag,         &
                        Write_Results_R_Samps,      &
                        Write_Results_T_Samps,      &
                        File_Suffix

USE Variables_MPI, &
                ONLY :  ierr

USE Variables_FP, &
                ONLY :  FP_Coeff_Vector,            &
                        FP_Coeff_Vector_Beta,       &
                        FP_Source_Vector_Beta

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
                        Open_New_File

USE IO_Print_Results, &
                ONLY :  Print_Results

USE IO_Convergence_Output, &
                ONLY :  Output_Convergence_Data

USE Poseidon_Main_Module, &
                ONLY :  Poseidon_CFA_Set_Uniform_Boundary_Conditions,       &
                        Poseidon_Run,                                       &
                        Poseidon_Close

USE FP_Method_Module, &
                ONLY :  Solve_FP_System,        &
                        Solve_FP_System_Beta

USE FP_Initial_Guess_Module, &
                ONLY :  Init_FP_Guess_Flat,     &
                        Init_FP_Guess_Informed, &
                        Input_FP_Guess


USE Source_Input_Module, &
                ONLY :  Poseidon_Input_Sources





USE MPI


IMPLICIT NONE

!                                       !
!   Poseidon Initialization Variables   !
!                                       !
INTEGER                                                 ::  Test_Number

INTEGER                                                 ::  FEM_Degree_Input
INTEGER                                                 ::  L_Limit_Input

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
LOGICAL                                                 ::  Flat_Guess

INTEGER,   DIMENSION(5)                                 ::  CFA_EQs



REAL(idp)                                               ::  Shift_Vector_BC
CHARACTER(LEN=1), DIMENSION(1:5)                        ::  INNER_BC_TYPES, OUTER_BC_TYPES
REAL(idp), DIMENSION(1:5)                               ::  INNER_BC_VALUES, OUTER_BC_VALUES

CHARACTER(LEN=10)                                       ::  Suffix_Input

REAL(idp)                                               ::  Alpha
REAL(idp)                                               ::  Beta
REAL(idp)                                               ::  C
REAL(idp)                                               ::  rho_o
REAL(idp)                                               ::  uaR
REAL(idp)                                               ::  fofalpha

REAL(idp)                                               ::  Star_Radius

REAL(idp)                                               ::  Psi_BC

REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Cur_R_Locs
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_R_Quad
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_T_Quad
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_P_Quad
REAL(idp)                                               ::  Left_Limit
REAL(idp)                                               ::  Right_Limit

REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Local_E
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Local_S
REAL(idp), DIMENSION(:,:,:,:,:), ALLOCATABLE            ::  Local_Si

REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Test_Solution

INTEGER                                                 ::  myID

INTEGER                                                 ::  re, rq


INTEGER                                                 ::  HCT_Fileid
CHARACTER(LEN = 100)                                    ::  HCT_Filename

CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)




!############################################################!
!#                                                          #!
!#                      Test Parameters                     #!
!#                                                          #!
!############################################################!
Units_Input         = "G"
Solver_Type         = 2


FEM_Degree_Input    = 1
L_Limit_Input       = 0


Alpha               =  sqrt(5.0_idp)
Star_Radius         =  1.0E+9_idp               ! (cm)

Convergence_Criteria = 1.0E-19_idp

Dimension_Input     = 3

Max_Iterations      = 300

Mesh_Type           = 1
Domain_Edge(1)      = 0.0_idp                   ! Inner Radius (cm)
Domain_Edge(2)      = 8.0_idp * Star_Radius     ! Outer Radius (cm)


NE(1)               = 128                       ! Number of Radial Elements
NE(2)               = 1                         ! Number of Theta Elements
NE(3)               = 1                         ! Number of Phi Elements

NQ(1)               = 10                        ! Number of Radial Quadrature Points
NQ(2)               = 1                         ! Number of Theta Quadrature Points
NQ(3)               = 1                         ! Number of Phi Quadrature Points

!Verbose             = .FALSE.
Verbose             = .TRUE.
Flat_Guess          = .FALSE.
Suffix_Input        = "Params"

CFA_Eqs = [1, 0, 0, 0, 0]


Write_Results_R_Samps = 256
Write_Results_T_Samps = 1


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

ALLOCATE( Test_Solution(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 ) )

Input_R_Quad = Initialize_LG_Quadrature_Locations(NQ(1))
Input_T_Quad = Initialize_LG_Quadrature_Locations(NQ(2))
Input_P_Quad = Initialize_LG_Quadrature_Locations(NQ(3))

Left_Limit  = -0.50_idp
Right_Limit = +0.50_idp

Input_R_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_R_Quad)
Input_T_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_T_Quad)
Input_P_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_P_Quad)




CALL Open_Run_Report_File()


CALL Set_Units(Units_Input)



Star_Radius = Star_Radius*Centimeter
Domain_Edge = Domain_Edge*Centimeter



CALL Create_3D_Mesh( Mesh_Type,         &
                     Domain_Edge(1),    &
                     Domain_Edge(2),    &
                     NE(1),             &
                     NE(2),             &
                     NE(3),             &
                     x_e, x_c, dx_c,    &
                     y_e, y_c, dy_c,    &
                     z_e, z_c, dz_c     )



!############################################################!
!#                                                          #!
!#                   Initialize Poseidon                    #!
!#                                                          #!
!############################################################!
CALL Initialize_Poseidon &
    (   Dimensions_Option       = Dimension_Input,  &
        FEM_Degree_Option       = FEM_Degree_Input, &
        L_Limit_Option          = L_Limit_Input,    &
        Units_Option            = Units_Input,      &
        Domain_Edge_Option      = Domain_Edge,      &
        NE_Option               = NE,               &
        NQ_Option               = NQ,               &
        r_Option                = x_e,              &
        t_Option                = y_e,              &
        p_Option                = z_e,              &
!        dr_Option               = dx_c,             &
!        dt_Option               = dy_c,             &
!        dp_Option               = dz_c              &
        Suffix_Flag_Option       = Suffix_Input,    &
        Solver_Type_Option       = Solver_Type,     &
        CFA_Eq_Flags_Option      = CFA_Eqs,         &
        Max_Iterations_Option    = Max_Iterations,  &
        Verbose_Option           = Verbose          )






!############################################################!
!#                                                          #!
!#               Create & Input Source Values               #!
!#                                                          #!
!############################################################!


WRITE(HCT_Filename,'(A,A,A,A)') Poseidon_Results_Dir,"HCT_Params_",TRIM(File_Suffix),".out"
Call Open_New_File(HCT_Filename,HCT_Fileid)
WRITE(HCT_Fileid,'(2ES22.15)') Alpha, Star_Radius/Centimeter
CLOSE( HCT_Fileid)



fofalpha            =  Alpha**5/(1.0_idp+Alpha*Alpha)**3

rho_o               =  (3.0_idp/(2.0_idp*pi*Star_Radius*Star_Radius) )*fofalpha*fofalpha
uaR                 =  sqrt(Alpha/((1.0_idp+Alpha*Alpha)*Star_Radius))
C                   =  1.0_idp/sqrt(sqrt( (2.0_idp/3.0_idp)*pi*rho_o  ) )
Beta                =  (C*uaR-1.0_idp)*Star_Radius

Print*,"fofalpha",fofalpha
print*,"rho_o",rho_o
Print*,"uaR",uaR
print*,"C",C
print*,"Beta",beta



DO re = 1,NE(1)
    Cur_R_Locs(:) = dx_c(re)*(Input_R_Quad(:) + Left_Limit)+  x_e(re)

    DO rq = 1,NQ(1)
        IF ( cur_r_locs(rq) .LE. Star_Radius ) THEN
            Local_E(rq, re-1, 0, 0) = rho_o
!            Print*,Local_E(rq, re-1, 0, 0)
        ELSE
            Local_E(rq, re-1, 0, 0) = 0.0_idp
        END IF
!        PRINT*,Cur_r_locs(rq), STar_radius,Local_E(rq, re-1, 0, 0)
    END DO ! rq
END DO ! re


Local_S  = 0.0_idp
Local_Si = 0.0_idp

CALL Poseidon_Input_Sources(    myID, myID, myID,                           &
                                Local_E, Local_S, Local_Si,                 &
                                NE(1), NE(2), NE(3),                        &
                                NQ(1), NQ(2), NQ(3),                        &
                                Input_R_Quad, Input_T_Quad, Input_P_Quad,   &
                                Left_Limit, Right_Limit                     )










!############################################################!
!#                                                          #!
!#          Calculate and Set Boundary Conditions           #!
!#                                                          #!
!############################################################!


Psi_BC = HCT_Solution( Domain_Edge(2), Alpha, Beta, C, Star_Radius )


INNER_BC_TYPES = (/"N", "N","N","N","N"/)
OUTER_BC_TYPES = (/"D", "D","D","D","D"/)


INNER_BC_VALUES = (/0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp /)
OUTER_BC_VALUES = (/Psi_BC,  0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp /)


CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("I", INNER_BC_TYPES, INNER_BC_VALUES)
CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("O", OUTER_BC_TYPES, OUTER_BC_VALUES)









!############################################################!
!#                                                          #!
!#              Calculate and Set Initial Guess             #!
!#                                                          #!
!############################################################!

IF ( Flat_Guess ) THEN

    CALL Init_FP_Guess_flat()

ELSE
!     Lapse function Coefficients in 1D correspond to the value of the function
!     at the location of the FEM nodes.

    DO re = 1,NE(1)
        Cur_R_Locs(:) = dx_c(re)*(Input_R_Quad(:) + Left_Limit)+  x_e(re)
        DO rq = 1,NQ(1)

            Test_Solution(rq, re-1, 0, 0) = HCT_Solution( Cur_R_Locs(rq), Alpha, Beta, C, Star_Radius )
!            PRINT*,Cur_R_locs(rq),Test_Solution(rq, re-1,0,0)

        END DO ! rq
    END DO ! re



    CALL Input_FP_Guess( Test_Solution,                             &
                         NE(1), NE(2), NE(3),                       &
                         NQ(1), NQ(2), NQ(3),                       &
                         Input_R_Quad, Input_T_Quad, Input_P_Quad,  &
                         Left_Limit, Right_Limit                    )

END IF





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

IF (Verbose .EQV. .TRUE. ) THEN
    CALL Print_Results()
END IF


!Write_Results_Flag = 1
IF ( Write_Results_Flag == 1 ) THEN
    CALL Output_Convergence_Data()
    CALL Output_Final_Results()
END IF



!############################################################!
!#                                                          #!
!#                      Close Poseidon                      #!
!#                                                          #!
!############################################################!

CALL Poseidon_Close()

DEALLOCATE( Cur_R_Locs )
DEALLOCATE( x_e, y_e, z_e )
DEALLOCATE( x_c, y_c, z_c )
DEALLOCATE( dx_c, dy_c, dz_c )

DEALLOCATE( Local_E, Local_S, Local_Si )

DEALLOCATE( Test_Solution )





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













END PROGRAM Hamiltonian_Constraint_Test

