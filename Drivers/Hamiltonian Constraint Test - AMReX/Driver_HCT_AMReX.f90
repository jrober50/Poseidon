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

USE Initialization_AMReX, &
            ONLY :  Initialize_Poseidon_with_AMReX

USE Variables_IO, &
           ONLY :  Write_Results_R_Samps,      &
                   Write_Results_T_Samps

USE Variables_MPI, &
            ONLY :  ierr

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature_Locations

USE Maps_X_Space, &
            ONLY :  Map_From_X_Space

USE Poseidon_Main_Module, &
            ONLY :  Poseidon_Close

USE Driver_Run_Module, &
            ONLY :  Driver_Run

USE Driver_SetSource_Module, &
            ONLY:   Driver_SetSource

USE Driver_SetBC_Module, &
            ONLY:   Driver_SetBC

USE Driver_SetGuess_Module, &
            ONLY:   Driver_SetGuess

USE Poseidon_AMReX_Input_Parsing_Module, &
            ONLY :  Init_AMReX_Parameters

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

INTEGER                                                 ::  M_Index
INTEGER                                                 ::  M_Index_Min
INTEGER                                                 ::  M_Index_Max


REAL(idp)                                               ::  Alpha
REAL(idp)                                               ::  Star_Radius


INTEGER                                                 ::  IRL
INTEGER                                                 ::  IFL
INTEGER, DIMENSION(1:8)                                 ::  Anderson_M_Values
CHARACTER(LEN=1), DIMENSION(1:10)                       ::  Letter_Table

Letter_Table = (/ "A","B","C","D","E","F","G","H","I","J" /)


!############################################################!
!#                                                          #!
!#                      Test Parameters                     #!
!#                                                          #!
!############################################################!
Units_Input         = "G"

Anderson_M_Values = (/ 1, 2, 3, 4, 5, 10, 20, 50 /)

M_Index_Min         =  3
M_Index_Max         =  3

IFL                 =  0
IRL                 =  0

Alpha               =  sqrt(5.0_idp)
Star_Radius         =  1.0E+5_idp               ! (cm)

NQ(1)               = 5                        ! Number of Radial Quadrature Points
NQ(2)               = 1                        ! Number of Theta Quadrature Points
NQ(3)               = 1                        ! Number of Phi Quadrature Points

CFA_Eqs = (/ 1, 1, 1, 1, 1 /)

Verbose             = .TRUE.
!Verbose             = .FALSE.
Suffix_Input        = "Params"

Write_Results_R_Samps = 256
Write_Results_T_Samps = 1




!############################################################!
!#                                                          #!
!#                 Initialize AMReX & MPI                   #!
!#                                                          #!
!############################################################!
CALL MPI_INIT(ierr)

CALL amrex_init()
CALL amrex_amrcore_init()
CALL Init_AMReX_Parameters()



!############################################################!
!#                                                          #!
!#                       Start Program                      #!
!#                                                          #!
!############################################################!
DO M_Index = M_Index_Min, M_Index_Max



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
    CALL Driver_SetSource( Alpha, Star_Radius )


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
    Call Driver_Run( )


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


