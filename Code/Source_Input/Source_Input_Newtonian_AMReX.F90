   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE  Source_Input_Newtonian_AMReX_Module                                  !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!    +101+   Poseidon_Input_Sources_Newtonian_Native                      !##!
!##!    +201+   Poseidon_Input_Sources_Newtonian_Native_Caller               !##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!                                                                         !##!
!###############################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !##########################################################################!


!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!

USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Poseidon_Message_Routines_Module, &
            ONLY :  Run_Message

USE Variables_Interface, &
            ONLY :  Caller_Set,                     &
                    Caller_NQ,                      &
                    Caller_Quad_DOF,                &
                    Caller_xL,                      &
                    Caller_RQ_xlocs,                &
                    Caller_TQ_xlocs,                &
                    Caller_PQ_xlocs,                &
                    Translation_Matrix

USE Variables_Source, &
            ONLY :  Source_Rho

USE Variables_Mesh, &
            ONLY :  NUM_R_ELEMENTS,     &
                    NUM_T_ELEMENTS,     &
                    NUM_P_ELEMENTS

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Num_T_Quad_Points,      &
                    Num_P_Quad_Points,      &
                    Num_TP_Quad_Points,     &
                    Local_Quad_DOF,         &
                    xLeftLimit,             &
                    xRightLimit,            &
                    Int_R_Locations

USE Functions_Translation_Matrix_Module, &
            ONLY :  Create_Translation_Matrix

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                         &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Poisson_SourceInput,          &
                    Timer_Poisson_SourceInput_PartA,    &
                    Timer_Poisson_SourceInput_PartB


#ifdef POSEIDON_AMREX_FLAG
use amrex_base_module

USE amrex_box_module,       ONLY: &
  amrex_box
USE amrex_multifab_module,  ONLY: &
  amrex_multifab, &
  amrex_multifab_build, &
  amrex_multifab_destroy, &
  amrex_mfiter, &
  amrex_mfiter_build, &
  amrex_mfiter_destroy


USE Variables_AMReX_Core, &
ONLY :  MF_Source,          &
        AMReX_Num_Levels,   &
        MF_Source_nComps
#endif


IMPLICIT NONE




CONTAINS
#ifdef POSEIDON_AMREX_FLAG
!+101+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources_Newtonian_AMReX(  MF_Src_Input,           &
                                                    MF_Src_Input_nComps,    &
                                                    Num_Levels,             &
                                                    Input_NQ,               &
                                                    Input_R_Quad,           &
                                                    Input_T_Quad,           &
                                                    Input_P_Quad,           &
                                                    Input_xL                )

TYPE(amrex_multifab),                   INTENT(IN)  ::  MF_Src_Input(0:Num_Levels-1)
INTEGER,                                INTENT(IN)  ::  MF_Src_Input_nComps
INTEGER,                                INTENT(IN)  ::  Num_Levels

INTEGER,    DIMENSION(3),               INTENT(IN)  ::  Input_NQ
REAL(idp),  DIMENSION(Input_NQ(1)),     INTENT(IN)  ::  Input_R_Quad
REAL(idp),  DIMENSION(Input_NQ(2)),     INTENT(IN)  ::  Input_T_Quad
REAL(idp),  DIMENSION(Input_NQ(3)),     INTENT(IN)  ::  Input_P_Quad
REAL(idp),  DIMENSION(2),               INTENT(IN)  ::  Input_xL

INTEGER                                             ::  RE, TE, PE
INTEGER,    DIMENSION(3)                            ::  iEL, iEU
INTEGER                                             ::  level


INTEGER                                             ::  Index
INTEGER                                             ::  Here
INTEGER                                             ::  There
INTEGER                                             ::  Local_Here

REAL(idp), CONTIGUOUS, POINTER                      ::  My_PTR(:,:,:,:)
REAL(idp), CONTIGUOUS, POINTER                      ::  Their_PTR(:,:,:,:)

TYPE(amrex_mfiter)                                  ::  mfi
TYPE(amrex_box)                                     ::  Box

REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  TransMat
INTEGER                                             ::  Their_DOF



IF ( Verbose_Flag ) CALL Run_Message('Receiving Newtonian Sources. Container : AMReX Multifab.')
CALL TimerStart(Timer_Poisson_SourceInput)



! Define Interpolation Matrix
Their_DOF = Input_NQ(1)*Input_NQ(2)*Input_NQ(3)


ALLOCATE(TransMat(1:Their_DOF, 1:Local_Quad_DOF))

TransMat =  Create_Translation_Matrix(  Input_NQ,          &
                                        Input_xL,          &
                                        Input_R_Quad,    &
                                        Input_T_Quad,    &
                                        Input_P_Quad,    &
                                        Their_DOF,         &
                                        [Num_R_Quad_Points, Num_T_Quad_Points, Num_P_Quad_Points ],            &
                                        [xLeftLimit, xRightLimit ],            &
                                        Int_R_Locations,      &
                                        Int_R_Locations,      &
                                        Int_R_Locations,      &
                                        Local_Quad_DOF            )



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

        lPF_SI_Flags(iPF_SI_MF_Ready) = .TRUE.
    END DO
END IF



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
        DO Local_Here = 1,Local_Quad_DOF

            My_PTR(re,te,pe,Local_Here) = DOT_PRODUCT( TransMat(:,Local_Here),          &
                                                       Their_PTR(re,te,pe,1:Their_DOF)  )

        END DO  ! Local_Here
        END DO  ! pe
        END DO  ! te
        END DO  ! re

    END DO      ! mfi
END DO          ! level



CALL TimerStop(Timer_Poisson_SourceInput)



END SUBROUTINE Poseidon_Input_Sources_Newtonian_AMReX






!+201+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources_Newtonian_AMReX_Caller( MF_Src_Input,           &
                                                          MF_Src_Input_nComps     )

TYPE(amrex_multifab),                   INTENT(IN)  ::  MF_SRC_Input(0:AMReX_Num_Levels-1)
INTEGER,                                INTENT(IN)  ::  MF_Src_Input_nComps

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




IF ( Verbose_Flag ) CALL Run_Message('Receiving XCFC Sources. Container : AMReX Multifab.')
CALL TimerStart(Timer_GR_SourceInput)

IF ( .NOT. lPF_SI_Flags(iPF_SI_MF_Ready) ) THEN
    DO level = 0,AMReX_Num_Levels-1

        CALL amrex_multifab_build(  MF_Source(level),           &
                                    MF_Src_Input(Level)%BA,     &
                                    MF_Src_Input(Level)%DM,     &
                                    MF_Source_nComps, 1         )

        lPF_SI_Flags(iPF_SI_MF_Ready) = .TRUE.
    END DO
END IF



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
        DO Local_Here = 1,Local_Quad_DOF

            My_PTR(re,te,pe,Local_Here) = DOT_PRODUCT( Translation_Matrix(:,Local_Here),    &
                                                       Their_PTR(re,te,pe,1:Their_DOF)      )

        END DO  ! Local_Here
        END DO  ! pe
        END DO  ! te
        END DO  ! re

    END DO      ! mfi
END DO          ! level



CALL TimerStop(Timer_Poisson_SourceInput)




END SUBROUTINE Poseidon_Input_Sources_Newtonian_AMReX_Caller






#else


SUBROUTINE Poseidon_Input_Sources_Newtonian_AMReX( )
END SUBROUTINE Poseidon_Input_Sources_Newtonian_AMReX


SUBROUTINE Poseidon_Input_Sources_Newtonian_AMReX_Caller( A_Difference )
INTEGER     :: A_Difference
END SUBROUTINE Poseidon_Input_Sources_Newtonian_AMReX_Caller



#endif



END MODULE Source_Input_Newtonian_AMReX_Module

