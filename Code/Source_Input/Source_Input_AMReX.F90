  !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Source_Input_AMReX_Module                                               !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  !################################################################################!


USE Poseidon_Kinds_Module, &
               ONLY :  idp

#ifdef POSEIDON_AMREX_FLAG
use amrex_base_module

USE amrex_amrcore_module, &
            ONLY :  amrex_get_numlevels

USE amrex_box_module, &
            ONLY :  amrex_box
            
USE amrex_multifab_module, &
            ONLY :  amrex_multifab,         &
                    amrex_multifab_build,   &
                    amrex_multifab_destroy, &
                    amrex_mfiter,           &
                    amrex_mfiter_build,     &
                    amrex_mfiter_destroy

USE Variables_AMReX_Core, &
            ONLY :  MF_Source,          &
                    AMReX_Num_Levels,   &
                    AMReX_Old_Levels,   &
                    MF_Source_nComps

USE Initialization_XCFC_with_AMReX_Module, &
            ONLY :  Initialization_XCFC_with_AMReX,     &
                    Reinitialization_XCFC_with_AMReX
                    
USE Initialization_Poisson, &
            ONLY :  Initialization_Poisson_with_AMReX,     &
                    Reinitialization_Poisson_with_AMReX

USE Poseidon_AMReX_Utilities_Module, &
            ONLY :  Destroy_MF,         &
                    Build_MF_from_MF
#endif


USE Poseidon_Parameters, &
            ONLY :  Domain_Dim,         &
                    Verbose_Flag
                    
USE Poseidon_IO_Parameters , &
            ONLY :  Mode_Names

USE Poseidon_Message_Routines_Module, &
            ONLY :  Run_Message


USE Parameters_Variable_Indices, &
            ONLY :  iS_E,               &
                    iS_S,               &
                    iS_S1,              &
                    iS_S2,              &
                    iS_S3

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Num_T_Quad_Points,      &
                    Num_P_Quad_Points,      &
                    Local_Quad_DOF,         &
                    Int_R_Locations,        &
                    Int_T_Locations,        &
                    Int_P_Locations,        &
                    xLeftLimit,             &
                    xRightLimit

USE Functions_Translation_Matrix_Module, &
            ONLY :  Create_Translation_Matrix

USE Timer_Routines_Module, &
            ONLY :  TimerStart,     &
                    TimerSTop

USE Timer_Variables_Module, &
            ONLY :  Timer_GR_SourceInput

USE Variables_Interface, &
            ONLY :  Caller_Quad_DOF,            &
                    Translation_Matrix

USE Flags_Source_Input_Module, &
            ONLY :  lPF_SI_Flags,       &
                    iPF_SI_MF_Ready

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Flags,     &
                    iPF_Init_Method_Vars
                    
USE Flags_Core_Module, &
            ONLY :  iPF_Core_Flags,             &
                    iPF_Core_Method_Mode,       &
                    iPF_Core_Method_XCFC,       &
                    iPF_Core_Method_Newtonian

use mpi



IMPLICIT NONE

CONTAINS


#ifdef POSEIDON_AMREX_FLAG
!+101+##########################################################################!
!                                                                               !
!          Poseidon_Input_Sources_AMREX                                         !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources_AMREX(   MF_Src_Input,          &
                                                Num_Levels,            &
                                                Input_NQ,              &
                                                Input_R_Quad,          &
                                                Input_T_Quad,          &
                                                Input_P_Quad,          &
                                                Input_xL,               &
                                                Remesh_Flag_Option      )

TYPE(amrex_multifab),                   INTENT(IN)  ::  MF_Src_Input(0:Num_Levels-1)
INTEGER,                                INTENT(IN)  ::  Num_Levels

INTEGER,    DIMENSION(3),               INTENT(IN)  ::  Input_NQ
REAL(idp),  DIMENSION(Input_NQ(1)),     INTENT(IN)  ::  Input_R_Quad
REAL(idp),  DIMENSION(Input_NQ(2)),     INTENT(IN)  ::  Input_T_Quad
REAL(idp),  DIMENSION(Input_NQ(3)),     INTENT(IN)  ::  Input_P_Quad
REAL(idp),  DIMENSION(2),               INTENT(IN)  ::  Input_xL
LOGICAL,                    OPTIONAL,   INTENT(IN)  ::  Remesh_Flag_Option


INTEGER                                             ::  level
LOGICAL                                             ::  All_Flag
LOGICAL,    DIMENSION(0:AMReX_Num_Levels-1)         ::  Remesh_Flag

CHARACTER(LEN=300)                                  ::  MessageBox

IF ( Verbose_Flag ) THEN
    WRITE(MessageBox,'(A,A,A)')'Receiving ',trim(Mode_Names(iPF_Core_Flags(iPF_Core_Method_Mode))),' Sources. Container : AMReX Multifab.'
    CALL Run_Message(TRIM(MessageBox))
END IF

        
CALL TimerStart(Timer_GR_SourceInput)

IF ( PRESENT(Remesh_Flag_Option) ) THEN
    Remesh_Flag = Remesh_Flag_Option
ELSE
    Remesh_Flag = .FALSE.
END IF




! Define Interpolation Matrix
Caller_Quad_DOF = Input_NQ(1)*Input_NQ(2)*Input_NQ(3)


ALLOCATE(Translation_Matrix(1:Caller_Quad_DOF, 1:Local_Quad_DOF))


Translation_Matrix = Create_Translation_Matrix( Input_NQ,               &
                                                Input_xL,               &
                                                Input_R_Quad,           &
                                                Input_T_Quad,           &
                                                Input_P_Quad,           &
                                                Caller_Quad_DOF,        &
                                                [Num_R_Quad_Points,     &
                                                    Num_T_Quad_Points,  &
                                                    Num_P_Quad_Points ],&
                                                [xLeftLimit,            &
                                                    xRightLimit ],      &
                                                Int_R_Locations,        &
                                                Int_R_Locations,        &
                                                Int_R_Locations,        &
                                                Local_Quad_DOF          )


!
!   If the mesh is being defined or being redefined, MF_Source will need to be
!   built/rebuilt.
!
IF ( .NOT. lPF_SI_Flags(iPF_SI_MF_Ready) ) THEN
    DO level = 0,AMReX_Num_Levels-1

        CALL amrex_multifab_build(  MF_Source(level),           &
                                    MF_Src_Input(Level)%BA,     &
                                    MF_Src_Input(Level)%DM,     &
                                    MF_Source_nComps, 1         )

    END DO
    lPF_SI_Flags(iPF_SI_MF_Ready) = .TRUE.
END IF



! Transfer data from MF_Src_Input to MF_Source.
All_Flag = .TRUE.
CALL Copy_Source_Data( MF_Src_Input, All_Flag )


CALL TimerStop(Timer_GR_SourceInput)


IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_XCFC ) THEN
    CALL Initialization_XCFC_with_AMReX()
ELSE IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_Newtonian ) THEN
    CALL Initialization_Poisson_with_AMReX()
END IF



END SUBROUTINE Poseidon_Input_Sources_AMREX


!+102+##########################################################################!
!                                                                               !
!          Poseidon_Input_Sources_AMREX_Caller                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources_AMREX_Caller( MF_Src_Input )

TYPE(amrex_multifab),                   INTENT(IN)  ::  MF_Src_Input(0:)

INTEGER                                             ::  level
LOGICAL                                             ::  All_Flag
LOGICAL,  DIMENSION(0:AMReX_Num_Levels-1)           ::  Remesh_Flag
CHARACTER(LEN=300)                                  ::  MessageBox

IF ( Verbose_Flag ) THEN
    WRITE(MessageBox,'(A,A,A)')'Receiving ',trim(Mode_Names(iPF_Core_Flags(iPF_Core_Method_Mode))),' Sources. Container : AMReX Multifab.'
    CALL Run_Message(TRIM(MessageBox))
END IF


CALL TimerStart(Timer_GR_SourceInput)




! Check if MF_Source exists.
! If not, created it to match MF_Src_Input and copy data.
IF ( .NOT. lPF_SI_Flags(iPF_SI_MF_Ready) ) THEN

    AMReX_Num_Levels = amrex_get_numlevels()
    DO level = 0,AMReX_Num_Levels-1

        CALL amrex_multifab_build(  MF_Source(level),           &
                                    MF_Src_Input(Level)%BA,     &
                                    MF_Src_Input(Level)%DM,     &
                                    MF_Source_nComps, 0         )

        
    END DO
    
        ! Transfer data from MF_Src_Input to MF_Source.

    All_Flag = .TRUE.
    CALL Copy_Source_Data( MF_Src_Input, All_Flag )
    
    lPF_SI_Flags(iPF_SI_MF_Ready) = .TRUE.

    CALL TimerStop(Timer_GR_SourceInput)


    IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_XCFC ) THEN
        CALL Initialization_XCFC_with_AMReX()
    ELSE IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_Newtonian ) THEN
        CALL Initialization_Poisson_with_AMReX()
    END IF



ELSE    ! Check if MF_Source has the same domain decomop as MF_Src_Input.
        ! If the do not match, then Poseidon will assumed the mesh has changed,
        ! and will preform a remeshing.
    

    Remesh_Flag = .FALSE.
    DO level = 0,AMReX_Num_Levels-1
        Remesh_Flag(level) = Multifab_Issame( MF_Source(level),   &
                                              MF_Src_Input(level) )
    END DO

    IF ( .NOT. ALL(Remesh_Flag) ) THEN

        ! Destroy the Old
        DO level = 0,AMReX_Num_Levels-1
            CALL amrex_multifab_destroy( MF_Source(level) )
        END DO
        
        AMReX_Num_Levels = amrex_get_numlevels()

        ! Rebuild to match the new.
        DO level = 0,AMReX_Num_Levels-1
            CALL amrex_multifab_build(  MF_Source(level),           &
                                        MF_Src_Input(Level)%BA,     &
                                        MF_Src_Input(Level)%DM,     &
                                        MF_Source_nComps, 0         )
        END DO

        ! Transfer data from MF_Src_Input to MF_Source.
        All_Flag = .TRUE.
        CALL Copy_Source_Data( MF_Src_Input, All_Flag )
        
        CALL TimerStop(Timer_GR_SourceInput)

        IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_XCFC ) THEN
            CALL Reinitialization_XCFC_with_AMReX()
        ELSE IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_Newtonian ) THEN
            CALL Reinitialization_Poisson_with_AMReX()
        END IF

    
    
    
    ELSE  ! MF_Source exists and matches the MF_Src_Input's decomposition.
          ! All that needs to be done is copy over the new data.

        ! Transfer data from MF_Src_Input to MF_Source.
        All_Flag = .TRUE.
        CALL Copy_Source_Data( MF_Src_Input, All_Flag )
        
        CALL TimerStop(Timer_GR_SourceInput)
    
    END IF ! All(Remesh_Flag)

END IF ! Existence of MF_Source

END SUBROUTINE Poseidon_Input_Sources_AMREX_Caller









!+201+##########################################################################!
!                                                                               !
!          Poseidon_Input_Sources_Part1_AMReX                                   !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources_Part1_AMReX(  MF_Src_Input,           &
                                                Num_Levels,             &
                                                Input_NQ,               &
                                                Input_R_Quad,           &
                                                Input_T_Quad,           &
                                                Input_P_Quad,           &
                                                Input_xL,               &
                                                Remesh_Flag_Option      )

TYPE(amrex_multifab),                   INTENT(IN)  ::  MF_Src_Input(0:Num_Levels-1)
INTEGER,                                INTENT(IN)  ::  Num_Levels

INTEGER,    DIMENSION(3),               INTENT(IN)  ::  Input_NQ
REAL(idp),  DIMENSION(Input_NQ(1)),     INTENT(IN)  ::  Input_R_Quad
REAL(idp),  DIMENSION(Input_NQ(2)),     INTENT(IN)  ::  Input_T_Quad
REAL(idp),  DIMENSION(Input_NQ(3)),     INTENT(IN)  ::  Input_P_Quad
REAL(idp),  DIMENSION(2),               INTENT(IN)  ::  Input_xL

LOGICAL,                    OPTIONAL,   INTENT(IN)  ::  Remesh_Flag_Option


INTEGER                                             ::  level
LOGICAL                                             ::  All_Flag
LOGICAL,    DIMENSION(0:AMReX_Num_Levels-1)         ::  Remesh_Flag
CHARACTER(LEN=300)                                  ::  MessageBox

IF ( Verbose_Flag ) THEN
    WRITE(MessageBox,'(A,A,A)')'Receiving ',trim(Mode_Names(iPF_Core_Flags(iPF_Core_Method_Mode))),' Part 1 Sources. Container : AMReX Multifab.'
    CALL Run_Message(TRIM(MessageBox))
END IF

CALL TimerStart(Timer_GR_SourceInput)




! Define Interpolation Matrix
Caller_Quad_DOF = Input_NQ(1)*Input_NQ(2)*Input_NQ(3)


ALLOCATE(Translation_Matrix(1:Caller_Quad_DOF, 1:Local_Quad_DOF))

Translation_Matrix = Create_Translation_Matrix( Input_NQ,               &
                                                Input_xL,               &
                                                Input_R_Quad,           &
                                                Input_T_Quad,           &
                                                Input_P_Quad,           &
                                                Caller_Quad_DOF,        &
                                                [Num_R_Quad_Points,     &
                                                    Num_T_Quad_Points,  &
                                                    Num_P_Quad_Points ],&
                                                [xLeftLimit,            &
                                                    xRightLimit ],      &
                                                Int_R_Locations,        &
                                                Int_R_Locations,        &
                                                Int_R_Locations,        &
                                                Local_Quad_DOF          )


IF ( .NOT. lPF_SI_Flags(iPF_SI_MF_Ready) ) THEN
    DO level = 0,AMReX_Num_Levels-1
        CALL amrex_multifab_build(  MF_Source(level),           &
                                    MF_Src_Input(Level)%BA,     &
                                    MF_Src_Input(Level)%DM,     &
                                    MF_Source_nComps, 1         )

        lPF_SI_Flags(iPF_SI_MF_Ready) = .TRUE.
    END DO
END IF



! Transfer data from MF_Src_Input to MF_Source.
All_Flag = .False.
CALL Copy_Source_Data( MF_Src_Input, All_Flag )




CALL TimerStop(Timer_GR_SourceInput)

IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_XCFC ) THEN
    CALL Initialization_XCFC_with_AMReX()
ELSE IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_Newtonian ) THEN
    CALL Initialization_Poisson_with_AMReX()
END IF

END SUBROUTINE Poseidon_Input_Sources_Part1_AMReX





!+101+##########################################################################!
!                                                                               !
!          Poseidon_Input_Sources_Part1_AMReX_Caller                            !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources_Part1_AMReX_Caller( MF_Src_Input   )

TYPE(amrex_multifab),               INTENT(IN)  ::  MF_Src_Input(0:)

INTEGER                                         ::  level
LOGICAL                                         ::  All_Flag
LOGICAL, ALLOCATABLE, DIMENSION(:)              ::  Remesh_Flag
INTEGER                                         ::  AMReX_New_Levels
CHARACTER(LEN=300)                              ::  MessageBox

IF ( Verbose_Flag ) THEN
    WRITE(MessageBox,'(A,A,A)')'Receiving ',trim(Mode_Names(iPF_Core_Flags(iPF_Core_Method_Mode))),' Part 1 Sources. Container : AMReX Multifab.'
    CALL Run_Message(TRIM(MessageBox))
END IF
CALL TimerStart(Timer_GR_SourceInput)



! Check if MF_Source exists.
! If not, create it to match MF_Src_Input and copy data.
IF ( .NOT. lPF_SI_Flags(iPF_SI_MF_Ready) ) THEN

    AMReX_Num_Levels = amrex_get_numlevels()
    DO level = 0,AMReX_Num_Levels-1
        CALL amrex_multifab_build(  MF_Source(level),           &
                                    MF_Src_Input(Level)%BA,     &
                                    MF_Src_Input(Level)%DM,     &
                                    MF_Source_nComps, 0         )

        
    END DO
    ! Transfer data from MF_Src_Input to MF_Source.
    All_Flag = .False.
    CALL Copy_Source_Data( MF_Src_Input, All_Flag )
    
    
    lPF_SI_Flags(iPF_SI_MF_Ready) = .TRUE.

    CALL TimerStop(Timer_GR_SourceInput)

    IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_XCFC ) THEN
        CALL Initialization_XCFC_with_AMReX()
    ELSE IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_Newtonian ) THEN
        CALL Initialization_Poisson_with_AMReX()
    END IF

ELSE    ! Check if MF_Source has the same domain decomop as MF_Src_Input.
        ! If the do not match, then Poseidon will assumed the mesh has changed,
        ! and will preform a remeshing.
    
    AMReX_New_Levels = amrex_get_numlevels()
    ALLOCATE( Remesh_Flag(0:AMReX_New_Levels-1) )
    Remesh_Flag = .FALSE.
    DO level = 0,AMReX_New_Levels-1
        Remesh_Flag(level) = Multifab_Issame( MF_Source(level),   &
                                              MF_Src_Input(level) )
    END DO


    IF ( .NOT. ALL(Remesh_Flag) ) THEN
        ! Destroy the Old
        DO level = 0,AMReX_Num_Levels-1
            CALL amrex_multifab_destroy( MF_Source(level) )
        END DO
        
        AMReX_Num_Levels = amrex_get_numlevels()
        ! Rebuild to match the new.
        DO level = 0,AMReX_Num_Levels-1
            CALL amrex_multifab_build(  MF_Source(level),           &
                                        MF_Src_Input(Level)%BA,     &
                                        MF_Src_Input(Level)%DM,     &
                                        MF_Source_nComps, 0         )
        END DO
        
        
        ! Transfer data from MF_Src_Input to MF_Source.
        All_Flag = .False.
        CALL Copy_Source_Data( MF_Src_Input, All_Flag )
        
        CALL TimerStop(Timer_GR_SourceInput)

        CALL Reinitialization_XCFC_with_AMReX()
    
    
    
    ELSE  ! MF_Source exists and matches the MF_Src_Input's decomposition.
          ! All that needs to be done is copy over the new data.
    
        ! Transfer data from MF_Src_Input to MF_Source.
        All_Flag = .False.
        CALL Copy_Source_Data( MF_Src_Input, All_Flag )
        
        CALL TimerStop(Timer_GR_SourceInput)
    
    END IF ! All(Remesh_Flag)

END IF ! Existence of MF_Source





END SUBROUTINE Poseidon_Input_Sources_Part1_AMReX_Caller






 !+301+####################################################!
!                                                           !
!          Copy_Source_Data                                 !
!                                                           !
 !#########################################################!
SUBROUTINE Copy_Source_Data( MF_Src_Input,          &
                             ALL_Flag               )

TYPE(amrex_multifab),           INTENT(IN)  ::  MF_Src_Input(0:AMReX_Num_Levels-1)
LOGICAL, OPTIONAL,              INTENT(IN)  ::  All_Flag

REAL(idp), CONTIGUOUS, POINTER                      ::  My_PTR(:,:,:,:)
REAL(idp), CONTIGUOUS, POINTER                      ::  Their_PTR(:,:,:,:)

TYPE(amrex_mfiter)                                  ::  mfi
TYPE(amrex_box)                                     ::  Box

INTEGER                                             ::  RE, TE, PE
INTEGER,    DIMENSION(3)                            ::  iEL, iEU
INTEGER                                             ::  level


INTEGER                                             ::  si
INTEGER                                             ::  Index
INTEGER                                             ::  Here
INTEGER                                             ::  There
INTEGER                                             ::  Local_Here

DO level = 0,AMReX_Num_Levels-1
    CALL amrex_mfiter_build(mfi, MF_Source(level), tiling = .true. )

    DO WHILE(mfi%next())
        Their_PTR => MF_Src_Input(level)%dataPtr(mfi)
        My_PTR    => MF_Source(level)%dataPtr(mfi)

        My_PTR = 0.0_idp

        Box = mfi%tilebox()

        iEL = Box%lo
        iEU = Box%hi


        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)



        si = iS_E
        DO Local_Here = 1,Local_Quad_DOF

            Here  = (si-1)*Caller_Quad_DOF+1
            There = si*Caller_Quad_DOF

            Index = (si-1)*Local_Quad_DOF+Local_Here

            My_PTR(re,te,pe,Index) = DOT_PRODUCT( Translation_Matrix(:,Local_Here), &
                                                  Their_PTR(re,te,pe,Here:There)    )
                                                  
        END DO  ! Local_Here
        
        IF ( iPF_Core_Flags(iPF_Core_Method_Mode) .NE. iPF_Core_Method_Newtonian ) THEN
            IF ( All_Flag ) THEN
                si = iS_S
                DO Local_Here = 1,Local_Quad_DOF

                    Here  = (si-1)*Caller_Quad_DOF+1
                    There = si*Caller_Quad_DOF

                    Index = (si-1)*Local_Quad_DOF+Local_Here

                    My_PTR(re,te,pe,Index) = DOT_PRODUCT( Translation_Matrix(:,Local_Here), &
                                                          Their_PTR(re,te,pe,Here:There)    )
                END DO  ! Local_Here
            END IF ! All_Flag


        
            DO si = iS_S1,iS_S1+Domain_Dim-1
            DO Local_Here = 1,Local_Quad_DOF

                Here  = (si-1)*Caller_Quad_DOF+1
                There = si*Caller_Quad_DOF

                Index = (si-1)*Local_Quad_DOF+Local_Here
        
                My_PTR(re,te,pe,Index) = DOT_PRODUCT( Translation_Matrix(:,Local_Here), &
                                                      Their_PTR(re,te,pe,Here:There)    )

            END DO  ! Local_Here
            END DO  ! si
        END IF ! Newtonian Mode
        
        END DO  ! pe
        END DO  ! te
        END DO  ! re
    END DO      ! mfi
END DO          ! level


END SUBROUTINE Copy_Source_Data





 !+401+####################################################!
!                                                           !
!            Multifab_Issame                                !
!                                                           !
 !#########################################################!
LOGICAL FUNCTION Multifab_Issame(MFA, MFB)

TYPE(amrex_multifab),       INTENT(IN)  ::  MFA
TYPE(amrex_multifab),       INTENT(IN)  ::  MFB
LOGICAL                                 ::  Flag

Flag = .FALSE.
IF (MFA%owner .AND. MFB%owner) THEN
    IF (amrex_boxarray_issame(MFA%BA,MFB%BA) .AND.      &
        amrex_distromap_issame(MFA%DM,MFB%DM) ) THEN
        Flag = .TRUE.
    END IF
END IF

Multifab_Issame = Flag

END FUNCTION Multifab_Issame







#else


SUBROUTINE Poseidon_Input_Sources_AMReX( )
STOP "Warning: You are attempting to use an AMReX function while the pre-compiler flag, POSEIDON_AMREX_FLAG, is false."
END SUBROUTINE Poseidon_Input_Sources_AMReX

SUBROUTINE Poseidon_Input_Sources_Part1_AMReX( )
STOP "Warning: You are attempting to use an AMReX function while the pre-compiler flag, POSEIDON_AMREX_FLAG, is false."
END SUBROUTINE Poseidon_Input_Sources_Part1_AMReX



SUBROUTINE Poseidon_Input_Sources_AMReX_Caller( A_Difference )
INTEGER     :: A_Difference
STOP "Warning: You are attempting to use an AMReX function while the pre-compiler flag, POSEIDON_AMREX_FLAG, is false."
END SUBROUTINE Poseidon_Input_Sources_AMReX_Caller

SUBROUTINE Poseidon_Input_Sources_Part1_AMReX_Caller( A_Difference )
INTEGER     :: A_Difference
STOP "Warning: You are attempting to use an AMReX function while the pre-compiler flag, POSEIDON_AMREX_FLAG, is false."
END SUBROUTINE Poseidon_Input_Sources_Part1_AMReX_Caller

#endif











END MODULE Source_Input_AMReX_Module
