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
            ONLY:  amrex_amrcore_init

USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Numbers_Module, &
            ONLY :  pi

USE Poseidon_Units_Module, &
            ONLY :  Set_Units

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

USE Maps_X_Space, &
            ONLY :  Map_From_X_Space

USE Poseidon_Main_Module, &
            ONLY :  Poseidon_Run,                                       &
                    Poseidon_Close

USE Driver_SetSource_Module, &
            ONLY :  Driver_SetSource

USE Driver_SetBC_Module, &
            ONLY :  Driver_SetBC

USE Driver_SetGuess_Module, &
            ONLY :  Driver_SetGuess

USE Driver_ReturnTest,  &
            ONLY :  Return_Test

USE Poseidon_AMReX_Input_Parsing_Module, &
            ONLY :  Init_AMReX_Parameters

USE Allocation_Yahil_Profile, &
            ONLY :  Deallocate_Yahil_Profile

USE Variables_MPI, &
            ONLY :  myID_Poseidon,      &
                    MasterID_Poseidon,  &
                    nPROCS_Poseidon,    &
                    Poseidon_Comm_World

USE Variables_Driver_AMReX, &
            ONLY :  nLevels


USE MPI


IMPLICIT NONE

!                                       !
!   Poseidon Initialization Variables   !
!                                       !
CHARACTER(LEN = 1)                                      ::  Units_Input

INTEGER,   DIMENSION(5)                                 ::  CFA_EQs

LOGICAL                                                 ::  Verbose
CHARACTER(LEN=10)                                       ::  Suffix_Input
CHARACTER(LEN=1)                                        ::  Suffix_Tail

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

INTEGER                                                 ::  T_Index
INTEGER                                                 ::  T_Index_Min
INTEGER                                                 ::  T_Index_Max

REAL(idp)                                               ::  Kappa
REAL(idp)                                               ::  Gamma
REAL(idp), DIMENSION(3)                                 ::  Yahil_Params


INTEGER                                                 ::  IRL
INTEGER                                                 ::  IFL

CHARACTER(LEN=1), DIMENSION(1:10)                       ::  Letter_Table
REAL(idp), DIMENSION(1:6)                               ::  Time_Values
INTEGER, DIMENSION(1:2)                                 ::  L_Values

INTEGER                                                 ::  i

Letter_Table = (/ "A","B","C","D","E","F","G","H","I","J" /)





!############################################################!
!#                                                          #!
!#                      Test Parameters                     #!
!#                                                          #!
!############################################################!
Units_Input         = "G"

Time_Values         = (/ 51.0_idp, 15.0_idp, 5.0_idp, 1.50_idp, 0.5_idp, 0.05_idp /)
L_Values            = (/ 5, 10 /)

T_Index_Min         =  4
T_Index_Max         =  4

M_Index_Min         =  3
M_Index_Max         =  3

IFL                 =  0
IRL                 =  0


Kappa               = 953946015514834.4
Gamma               = 1.30_idp


NQ(1)               = 5                        ! Number of Radial Quadrature Points
NQ(2)               = 1                        ! Number of Theta Quadrature Points
NQ(3)               = 1                        ! Number of Phi Quadrature Points


Verbose             = .TRUE.
!Verbose             = .FALSE.
Suffix_Input        = "Params"



CFA_Eqs = (/ 1, 1, 1, 1, 1 /)


Write_Results_R_Samps = 256
Write_Results_T_Samps = 1




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





!############################################################!
!#                                                          #!
!#                       Start Program                      #!
!#                                                          #!
!############################################################!
CALL Init_AMReX_Parameters()




DO M_Index = M_Index_Min, M_Index_Max
DO T_Index = T_Index_Min, T_Index_Max


    Suffix_Tail = Letter_Table(nLevels)




    ALLOCATE( Input_R_Quad(1:NQ(1)) )
    ALLOCATE( Input_T_Quad(1:NQ(2)) )
    ALLOCATE( Input_P_Quad(1:NQ(3)) )

    Input_R_Quad = Initialize_LG_Quadrature_Locations(NQ(1))
    Input_T_Quad = Initialize_LG_Quadrature_Locations(NQ(2))
    Input_P_Quad = Initialize_LG_Quadrature_Locations(NQ(3))

    Left_Limit  = -0.50_idp
    Right_Limit = +0.50_idp

    Input_R_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_R_Quad)
    Input_T_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_T_Quad)
    Input_P_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_P_Quad)



    
    !############################################################!
    !#                                                          #!
    !#                   Initialize Poseidon                    #!
    !#                                                          #!
    !############################################################!
    CALL Initialize_Poseidon_with_AMReX &
       (    Source_NQ                           = NQ,                   &
            Source_xL                           = [Left_Limit, Right_Limit],    &
            Source_RQ_xlocs                     = Input_R_Quad,         &
            Source_TQ_xlocs                     = Input_T_Quad,         &
            Source_PQ_xlocs                     = Input_P_Quad,         &
            Source_Units                        = Units_Input,          &
            Source_Radial_Boundary_Units        = "cm",                 &
            Integration_NQ_Option               = NQ,                   &
            CFA_Eq_Flags_Option                 = CFA_Eqs,              &
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
            Print_Condition_Option              = .TRUE.,               &
            Write_Condition_Option              = .TRUE.,               &
            Suffix_Flag_Option                  = Suffix_Input,         &
            Suffix_Tail_Option                  = Suffix_Tail           )



    !############################################################!
    !#                                                          #!
    !#               Create & Input Source Values               #!
    !#                                                          #!
    !############################################################!
    Yahil_Params = [Time_Values(T_Index), Kappa, Gamma]
    CALL Driver_SetSource(  Yahil_Params, nLevels )


    !############################################################!
    !#                                                          #!
    !#          Calculate and Set Boundary Conditions           #!
    !#                                                          #!
    !############################################################!
!    PRINT*,"Before Driver_SetBC"
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
    CALL Poseidon_Run()



    !############################################################!
    !#                                                          #!
    !#                       Output Results                     #!
    !#                                                          #!
    !############################################################!
    !CALL Return_Test(nLevels, NQ, MF_Source)


    !############################################################!
    !#                                                          #!
    !#                      Close Poseidon                      #!
    !#                                                          #!
    !############################################################!
    CALL Poseidon_Close()
    CALL Deallocate_Yahil_Profile()

    DEALLOCATE( Input_R_Quad )
    DEALLOCATE( Input_T_Quad )
    DEALLOCATE( Input_P_Quad )




END DO ! T_Index
END DO ! M_Index




!############################################################!
!#                                                          #!
!#                     Close AMReX & MPI                    #!
!#                                                          #!
!############################################################!
call amrex_finalize()
CALL MPI_Finalize(ierr)


END PROGRAM AMReX_Mapping


