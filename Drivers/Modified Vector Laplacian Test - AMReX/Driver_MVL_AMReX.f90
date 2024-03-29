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
            ONLY:   amrex_init
USE amrex_amrcore_module, &
            ONLY:   amrex_amrcore_init

USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Numbers_Module, &
            ONLY :  pi

USE Poseidon_Units_Module, &
            ONLY :  Set_Units,      &
                    Centimeter

USE Parameters_Variable_Indices, &
            ONLY :  iVB_X,                      &
                    iVB_S,                      &
                    iU_CF,                      &
                    iU_LF,                      &
                    iU_S1,                      &
                    iU_S2,                      &
                    iU_S3,                      &
                    iU_X1,                      &
                    iU_X2,                      &
                    iU_X3

USE Poseidon_Interface_Initialization, &
            ONLY :  Initialize_Poseidon

USE Poseidon_Interface_Run, &
            ONLY :  Poseidon_Run

USE Poseidon_Interface_Close, &
            ONLY :  Poseidon_Close

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

USE IO_Print_Results, &
            ONLY :  Print_Results

USE Driver_InitSource_Module, &
            ONLY :  Driver_InitSource
            
USE Driver_SetSource_Module, &
            ONLY :  Driver_SetSource

USE Driver_SetBC_Module, &
            ONLY :  Driver_SetBC

USE Driver_SetGuess_Module, &
            ONLY :  Driver_SetGuess

USE IO_Print_Results, &
            ONLY :  Print_Single_Var_Results
            
USE Driver_Variables, &
            ONLY :  Driver_NQ,          &
                    Driver_RQ_xLocs,    &
                    Driver_xL

USE IO_Write_Final_Results, &
            ONLY :  Write_Final_Results

USE Poseidon_AMReX_Input_Parsing_Module, &
            ONLY :  Init_AMReX_Parameters

USE Variables_Driver_AMReX, &
            ONLY :  nLevels

USE External_IO_Test_Results_Module, &
            ONLY :  Print_MVL_Error


USE MPI


IMPLICIT NONE

!                                       !
!   Poseidon Initialization Variables   !
!                                       !
INTEGER                                                 ::  Dimension_Input

INTEGER                                                 ::  Mesh_Type
INTEGER                                                 ::  Solver_Type

CHARACTER(LEN = 1)                                      ::  Units_Input

INTEGER, DIMENSION(3)                                   ::  NQ


LOGICAL                                                 ::  Verbose
LOGICAL                                                 ::  Print_Results_Flag
LOGICAL                                                 ::  Print_Setup_Flag
LOGICAL                                                 ::  Print_Time_Flag
LOGICAL                                                 ::  Print_Cond_Flag

INTEGER,   DIMENSION(5)                                 ::  CFA_EQs


CHARACTER(LEN=10)                                       ::  Suffix_Input
CHARACTER(LEN=4)                                        ::  Suffix_Tail

REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_R_Quad
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_T_Quad
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_P_Quad
REAL(idp)                                               ::  Left_Limit
REAL(idp)                                               ::  Right_Limit

INTEGER                                                 ::  i
INTEGER                                                 ::  myID, nPROCS

INTEGER                                                 ::  Guess_Type

INTEGER                                                 ::  M_Index
INTEGER                                                 ::  M_Index_Min
INTEGER                                                 ::  M_Index_Max


INTEGER                                                 ::  IRL
INTEGER                                                 ::  IFL

INTEGER, DIMENSION(1:8)                                 ::  Anderson_M_Values
CHARACTER(LEN=1), DIMENSION(1:10)                       ::  Letter_Table





!############################################################!
!#                                                          #!
!#                      Test Parameters                     #!
!#                                                          #!
!############################################################!
Units_Input         = "U"
Solver_Type         = 3

Anderson_M_Values   = (/ 1, 2, 3, 4, 5, 10, 20, 50 /)

M_Index_Min         =  3
M_Index_Max         =  3

IFL                 =  0
IRL                 =  0

Guess_Type          =  1            !  1 = Flat, 2 = Educated, 3 = Perturbed Educated.

Dimension_Input     =  3

NQ(1)               =  5            ! Number of Radial Quadrature Points
NQ(2)               =  1            ! Number of Theta Quadrature Points
NQ(3)               =  1            ! Number of Phi Quadrature Points

Left_Limit  = -0.50_idp
Right_Limit = +0.50_idp


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



CFA_Eqs =  (/ 0, 0, 0, 0, 0 /)


Letter_Table = (/ "A","B","C","D","E","F","G","H","I","J" /)


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

CALL amrex_init()
CALL amrex_amrcore_init()


CALL Init_AMReX_Parameters()
CALL Set_Units( Units_Input )


Driver_xL(1) = Left_Limit
Driver_xL(2) = Right_Limit
Driver_NQ = NQ
ALLOCATE(Driver_RQ_xLocs(Driver_NQ(1)))



!############################################################!
!#                                                          #!
!#                       Start Program                      #!
!#                                                          #!
!############################################################!
DO M_Index = M_Index_Min, M_Index_Max


    WRITE(Suffix_Tail,'(A)') Letter_Table(nLevels)


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
    CALL Driver_InitSource( )
    
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
            Write_Setup_Option                  = .FALSE.,              &
            Print_Results_Option                = Print_Results_Flag,   &
            Write_Results_Option                = .FALSE.,               &
            Print_Timetable_Option              = Print_Time_Flag,      &
            Write_Timetable_Option              = .FALSE.,               &
            Write_Sources_Option                = .FALSE.,              &
            Print_Condition_Option              = Print_Cond_Flag,      &
            Write_Condition_Option              = .FALSE.,              &
            Suffix_Flag_Option                  = Suffix_Input,         &
            Suffix_Tail_Option                  = Suffix_Tail           )





    !############################################################!
    !#                                                          #!
    !#               Create & Input Source Values               #!
    !#                                                          #!
    !############################################################!
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
    CALL Poseidon_Run()



    !############################################################!
    !#                                                          #!
    !#                       Output Results                     #!
    !#                                                          #!
    !############################################################!
    IF ( Print_Results_Flag ) THEN
        CALL mpi_barrier(MPI_COMM_WORLD, ierr)
        IF ( myID == 1 ) THEN
            WRITE(*,'(A)')" Final Results "
        END IF
        DO i = 0,nProcs-1
            IF ( myID == i ) THEN
                WRITE(*,'(A,I1.1)')" Proc = ",myID
!                CALL Print_Single_Var_Results( iU_X1, iVB_X )
                CALL Print_MVL_Error()
                WRITE(*,'(A)')"-------------------------------"
            END IF
            CALL mpi_barrier(MPI_COMM_WORLD, ierr)
        END DO
        CALL mpi_barrier(MPI_COMM_WORLD, ierr)
    END IF
    CALL Write_Final_Results(u_Overide = (/ 1, 1, 1, 1, 1 /))


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


