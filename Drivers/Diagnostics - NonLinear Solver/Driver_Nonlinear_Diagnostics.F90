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

USE Poseidon_Interface_Initialization, &
            ONLY :  Initialize_Poseidon

USE Poseidon_Interface_Run, &
            ONLY :  Poseidon_Run

USE Poseidon_Interface_Close, &
            ONLY :  Poseidon_Close

USE Variables_IO, &
            ONLY :  Write_Results_R_Samps,      &
                    Write_Results_T_Samps

USE Functions_Mesh, &
            ONLY :  Create_3D_Mesh

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature_Locations

USE Maps_X_Space, &
            ONLY :  Map_From_X_Space

USE Driver_SetSource_Module, &
            ONLY :  Driver_SetSource

USE Driver_SetBC_Module, &
            ONLY :  Driver_SetBC

USE Driver_SetGuess_Module, &
            ONLY :  Driver_SetGuess

USE Poseidon_AMReX_Input_Parsing_Module, &
            ONLY :  Init_AMReX_Parameters

USE Allocation_Yahil_Profile, &
            ONLY :  Deallocate_Yahil_Profile

USE Variables_MPI, &
            ONLY :  ierr

USE Variables_Driver_AMReX, &
            ONLY :  nLevels,            &
                    MF_Driver_Source

USE Variables_AMReX_Core, &
            ONLY :  MF_Source

USE ADM_Mass_Module, &
            ONLY :  Calc_ADM_Mass
            
USE Poseidon_Memory_Routines, &
            ONLY :  Poseidon_Mark_Memory
            
            
#ifdef POSEIDON_MEMORY_FLAG
USE Memory_Variables_Module, &
            ONLY :  Memory_Start,               &
                    Memory_After_MPI_Init,      &
                    Memory_After_AMReX_Init,    &
                    Memory_After_AMReX_Finalize,&
                    Memory_End,                 &
                    Memory_Loop_Start,          &
                    Memory_Loop_Before_Init,    &
                    Memory_Loop_After_Init,     &
                    Memory_Loop_Before_Run,     &
                    Memory_Loop_After_Run,      &
                    Memory_Loop_Before_Close,   &
                    Memory_Loop_End,            &
                    Memory_HWM
                    
#endif
                    
USE Memory_IO_Module, &
            ONLY :  Output_Poseidon_Memory_Loop_Report, &
                    Output_Poseidon_Memory_Total_Report


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

INTEGER                                                 ::  Loop
INTEGER                                                 ::  Loop_Min
INTEGER                                                 ::  Loop_Max

REAL(idp)                                               ::  Kappa
REAL(idp)                                               ::  Gamma
REAL(idp), DIMENSION(3)                                 ::  Yahil_Params

REAL(idp)                                               ::  ADM_Mass

INTEGER                                                 ::  IRL
INTEGER                                                 ::  IFL

CHARACTER(LEN=1), DIMENSION(1:10)                       ::  Letter_Table
REAL(idp), DIMENSION(1:6)                               ::  Time_Values
INTEGER, DIMENSION(1:2)                                 ::  L_Values

Letter_Table = (/ "A","B","C","D","E","F","G","H","I","J" /)





!############################################################!
!#                                                          #!
!#                      Test Parameters                     #!
!#                                                          #!
!############################################################!
Units_Input         = "G"

Time_Values         = (/ 51.0_idp, 15.0_idp, 5.0_idp, 1.50_idp, 0.5_idp, 0.05_idp /)
L_Values            = (/ 5, 10 /)

T_Index_Min         =  1
T_Index_Max         =  1

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

Print_Results_Flag  = .TRUE.
!Print_Results_Flag  = .FALSE.

!Print_Setup_Flag    = .TRUE.
Print_Setup_Flag    = .FALSE.

!Print_Time_Flag     = .TRUE.
Print_Time_Flag     = .FALSE.

!Print_Cond_Flag     = .TRUE.
Print_Cond_Flag     = .FALSE.


Suffix_Input        = "Params"



CFA_Eqs = (/ 1, 1, 1, 1, 1 /)


Write_Results_R_Samps = 256
Write_Results_T_Samps = 1

#ifdef POSEIDON_MEMORY_FLAG
        CALL Poseidon_Mark_Memory(Memory_Start,Memory_HWM)
        PRINT*,"Beginning                    : ",Memory_Start
#endif


!############################################################!
!#                                                          #!
!#                 Initialize AMReX & MPI                   #!
!#                                                          #!
!############################################################!
CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nPROCS,ierr)


#ifdef POSEIDON_MEMORY_FLAG
        CALL Poseidon_Mark_Memory(Memory_After_MPI_Init,Memory_HWM)
        PRINT*,"After MPI Init               : ",Memory_After_MPI_Init
#endif


call amrex_init()
CALL amrex_amrcore_init()

#ifdef POSEIDON_MEMORY_FLAG
        CALL Poseidon_Mark_Memory(Memory_After_AMReX_Init,Memory_HWM)
        PRINT*,"After AMReX Init             : ",Memory_After_AMReX_Init
#endif



!############################################################!
!#                                                          #!
!#                       Start Program                      #!
!#                                                          #!
!############################################################!
CALL Init_AMReX_Parameters()


DO M_Index = M_Index_Min, M_Index_Max
DO T_Index = T_Index_Min, T_Index_Max

#ifdef POSEIDON_MEMORY_FLAG
        CALL Poseidon_Mark_Memory(Memory_Loop_Start,Memory_HWM)
        PRINT*,"Beginning of Loop            : ",Memory_Loop_Start
#endif

    Suffix_Tail = Letter_Table(T_Index)




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
#ifdef POSEIDON_MEMORY_FLAG
    CALL Poseidon_Mark_Memory(Memory_Loop_Before_Init,Memory_HWM)
    PRINT*,"Before Init                  : ",Memory_Loop_Before_Init
#endif
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
            Fixed_Point_Diagnostics_Option      = .TRUE.,               &
            Verbose_Option                      = Verbose,              &
            WriteAll_Option                     = .FALSE.,              &
            Print_Setup_Option                  = Print_Setup_Flag,     &
            Write_Setup_Option                  = .FALSE.,              &
            Print_Results_Option                = Print_Results_Flag,   &
            Write_Results_Option                = .TRUE.,               &
            Print_Timetable_Option              = Print_Time_Flag,      &
            Write_Timetable_Option              = .TRUE.,               &
            Write_Sources_Option                = .FALSE.,              &
            Print_Condition_Option              = Print_Cond_Flag,      &
            Write_Condition_Option              = .FALSE.,              &
            Write_FP_Diagnostics_Option         = .TRUE.,               &
            Suffix_Flag_Option                  = Suffix_Input,         &
            Suffix_Tail_Option                  = Suffix_Tail           )

#ifdef POSEIDON_MEMORY_FLAG
    CALL Poseidon_Mark_Memory(Memory_Loop_After_Init,Memory_HWM)
    PRINT*,"After Init                   : ",Memory_Loop_After_Init
#endif

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
#ifdef POSEIDON_MEMORY_FLAG
    CALL Poseidon_Mark_Memory(Memory_Loop_Before_Run,Memory_HWM)
    PRINT*,"Before Poseidon_Run          : ",Memory_Loop_Before_Run
#endif
    
    CALL Poseidon_Run()
#ifdef POSEIDON_MEMORY_FLAG
    CALL Poseidon_Mark_Memory(Memory_Loop_After_Run,Memory_HWM)
    PRINT*,"After Poseidon_Run           : ",Memory_Loop_After_Run
#endif



    !############################################################!
    !#                                                          #!
    !#                      Close Poseidon                      #!
    !#                                                          #!
    !############################################################!
#ifdef POSEIDON_MEMORY_FLAG
    CALL Poseidon_Mark_Memory(Memory_Loop_Before_Close,Memory_HWM)
    PRINT*,"Before Close                 : ",Memory_Loop_Before_Close
#endif
    
    CALL Poseidon_Close()
    CALL Deallocate_Yahil_Profile()
    
    DEALLOCATE( MF_Driver_Source )
    DEALLOCATE( Input_R_Quad )
    DEALLOCATE( Input_T_Quad )
    DEALLOCATE( Input_P_Quad )


#ifdef POSEIDON_MEMORY_FLAG
    CALL Poseidon_Mark_Memory(Memory_Loop_End,Memory_HWM)
    PRINT*,"After Deallocation, Loop End : ",Memory_Loop_End
#endif

    CALL Output_Poseidon_Memory_Loop_Report( Loop )

END DO ! T_Index
END DO ! M_Index



!############################################################!
!#                                                          #!
!#                     Close AMReX & MPI                    #!
!#                                                          #!
!############################################################!
call amrex_finalize()
#ifdef POSEIDON_MEMORY_FLAG
    CALL Poseidon_Mark_Memory(Memory_After_AMReX_Finalize,Memory_HWM)
    PRINT*," After AMReX_Finalize        : ",Memory_After_AMReX_Finalize
#endif

CALL MPI_Finalize(ierr)
    


#ifdef POSEIDON_MEMORY_FLAG
    CALL Poseidon_Mark_Memory(Memory_End,Memory_HWM)
    PRINT*," Fin                         : ",Memory_End,Memory_HWM
#endif

CALL Output_Poseidon_Memory_Total_Report()

END PROGRAM AMReX_Mapping


