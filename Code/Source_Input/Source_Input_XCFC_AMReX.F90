  !################################################################################!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Source_Input_XCFC_AMReX_Module                                                    !##!
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

USE amrex_box_module,       ONLY: &
  amrex_box
USE amrex_multifab_module,  ONLY: &
  amrex_multifab, &
  amrex_multifab_build, &
  amrex_multifab_destroy, &
  amrex_mfiter, &
  amrex_mfiter_build, &
  amrex_mfiter_destroy
#endif


USE Poseidon_Parameters, &
            ONLY :  Domain_Dim,         &
                    Verbose_Flag

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

#ifdef POSEIDON_AMREX_FLAG
USE Variables_AMReX_Core, &
            ONLY :  MF_Source,          &
                    AMReX_Num_Levels,   &
                    MF_Source_nComps

USE Initialization_XCFC_with_AMReX_Module, &
            ONLY :  Initialization_XCFC_with_AMReX


#endif


USE Timer_Routines_Module, &
            ONLY :  TimerStart,     &
                    TimerSTop

USE Timer_Variables_Module, &
            ONLY :  Timer_GR_SourceInput


USE Variables_Interface, &
            ONLY :  Caller_Quad_DOF,            &
                    Translation_Matrix

USE Maps_X_Space, &
            ONLY :  Map_To_X_Space

USE Flags_Source_Input_Module, &
            ONLY :  lPF_SI_Flags,       &
                    iPF_SI_MF_Ready

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Flags,     &
                    iPF_Init_Method_Vars

use mpi



IMPLICIT NONE

CONTAINS


#ifdef POSEIDON_AMREX_FLAG
!+101+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources_XCFC_AMREX(   MF_Src_Input,          &
                                                MF_Src_Input_nComps,         &
                                                Num_Levels,            &
                                                Input_NQ,              &
                                                Input_R_Quad,          &
                                                Input_T_Quad,          &
                                                Input_P_Quad,          &
                                                Input_xL               )

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


INTEGER                                             ::  si
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




IF ( Verbose_Flag ) CALL Run_Message('Receiving XCFC Sources. Container : AMReX Multifab.')
CALL TimerStart(Timer_GR_SourceInput)

! Define Interpolation Matrix
Their_DOF = Input_NQ(1)*Input_NQ(2)*Input_NQ(3)


ALLOCATE(TransMat(1:Their_DOF, 1:Local_Quad_DOF))


TransMat = Create_Translation_Matrix(   Input_NQ,          &
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



        DO si = 1,2+DOMAIN_DIM
        DO Local_Here = 1,Local_Quad_DOF

            Here  = (si-1)*Their_DOF+1
            There = si*Their_DOF

            Index = (si-1)*Local_Quad_DOF+Local_Here

            My_PTR(re,te,pe,Index) = DOT_PRODUCT( TransMat(:,Local_Here), &
                                                  Their_PTR(re,te,pe,Here:There)    )

        END DO  ! Local_Here
        END DO  ! si
        END DO  ! pe
        END DO  ! te
        END DO  ! re
    END DO      ! mfi
END DO          ! level


CALL TimerStop(Timer_GR_SourceInput)


CALL Initialization_XCFC_with_AMReX()



END SUBROUTINE Poseidon_Input_Sources_XCFC_AMREX


!+102+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources_XCFC_AMREX_Caller(    MF_Src_Input,           &
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



        DO si = 1,2+DOMAIN_DIM
        DO Local_Here = 1,Local_Quad_DOF

            Here  = (si-1)*Caller_Quad_DOF+1
            There = si*Caller_Quad_DOF

            Index = (si-1)*Local_Quad_DOF+Local_Here

            My_PTR(re,te,pe,Index) = DOT_PRODUCT( Translation_Matrix(:,Local_Here), &
                                                  Their_PTR(re,te,pe,Here:There)    )

        END DO  ! Local_Here
        END DO  ! si


        END DO  ! pe
        END DO  ! te
        END DO  ! re
    END DO      ! mfi
END DO          ! level


CALL TimerStop(Timer_GR_SourceInput)


CALL Initialization_XCFC_with_AMReX()


END SUBROUTINE Poseidon_Input_Sources_XCFC_AMREX_Caller









!+201+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources_Part1_AMReX( MF_Src_Input,          &
                                          MF_Src_Input_nComps,         &
                                          Num_Levels,            &
                                          Input_NQ,              &
                                          Input_R_Quad,          &
                                          Input_T_Quad,          &
                                          Input_P_Quad,          &
                                          Input_xL               )

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


INTEGER                                             ::  si
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


IF ( Verbose_Flag ) CALL Run_Message('Receiving Part 1 XCFC Sources. Container : AMReX Multifab.')
CALL TimerStart(Timer_GR_SourceInput)

! Define Interpolation Matrix
Their_DOF = Input_NQ(1)*Input_NQ(2)*Input_NQ(3)


ALLOCATE(TransMat(1:Their_DOF, 1:Local_Quad_DOF))

TransMat = Create_Translation_Matrix(   Input_NQ,                       &
                                        Input_xL,                       &
                                        Input_R_Quad,                   &
                                        Input_T_Quad,                   &
                                        Input_P_Quad,                   &
                                        Their_DOF,                      &
                                        [Num_R_Quad_Points, Num_T_Quad_Points, Num_P_Quad_Points ],            &
                                        [xLeftLimit, xRightLimit ],            &
                                        Int_R_Locations,                &
                                        Int_R_Locations,                &
                                        Int_R_Locations,                &
                                        Local_Quad_DOF                  )


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



        Box = mfi%tilebox()

        iEL = Box%lo
        iEU = Box%hi


        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)



        si = iS_E
        DO Local_Here = 1,Local_Quad_DOF

            Here  = (si-1)*Their_DOF+1
            There = si*Their_DOF

            Index = (si-1)*Local_Quad_DOF+Local_Here

            My_PTR(re,te,pe,Index) = DOT_PRODUCT( TransMat(:,Local_Here), &
                                                  Their_PTR(re,te,pe,Here:There)    )

        END DO  ! Local_Here
        


        DO si = iS_S1,iS_S3
        DO Local_Here = 1,Local_Quad_DOF

            Here  = (si-1)*Their_DOF+1
            There = si*Their_DOF

            Index = (si-1)*Local_Quad_DOF+Local_Here

            My_PTR(re,te,pe,Index) = DOT_PRODUCT( TransMat(:,Local_Here), &
                                                  Their_PTR(re,te,pe,Here:There)    )

        END DO  ! Local_Here
        END DO  ! si





        END DO  ! pe
        END DO  ! te
        END DO  ! re
    END DO      ! mfi
END DO          ! level


CALL TimerStop(Timer_GR_SourceInput)


CALL Initialization_XCFC_with_AMReX()

END SUBROUTINE Poseidon_Input_Sources_Part1_AMReX


!+101+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources_Part1_AMReX_Caller( MF_Src_Input,       &
                                                      MF_Src_Input_nComps       )

TYPE(amrex_multifab),                   INTENT(IN)  ::  MF_Src_Input(0:AMReX_Num_Levels-1)
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


IF ( Verbose_Flag ) CALL Run_Message('Receiving Part 1 XCFC Sources. Container : AMReX Multifab.')
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
        


        DO si = iS_S1,iS_S3
        DO Local_Here = 1,Local_Quad_DOF

            Here  = (si-1)*Caller_Quad_DOF+1
            There = si*Caller_Quad_DOF

            Index = (si-1)*Local_Quad_DOF+Local_Here

            My_PTR(re,te,pe,Index) = DOT_PRODUCT( Translation_Matrix(:,Local_Here), &
                                                  Their_PTR(re,te,pe,Here:There)    )

        END DO  ! Local_Here
        END DO  ! si


        END DO  ! pe
        END DO  ! te
        END DO  ! re
    END DO      ! mfi
END DO          ! level


CALL TimerStop(Timer_GR_SourceInput)


CALL Initialization_XCFC_with_AMReX()


END SUBROUTINE Poseidon_Input_Sources_Part1_AMReX_Caller









#else


SUBROUTINE Poseidon_Input_Sources_XCFC_AMReX( )
END SUBROUTINE Poseidon_Input_Sources_XCFC_AMReX

SUBROUTINE Poseidon_Input_Sources_Part1_AMReX( )
END SUBROUTINE Poseidon_Input_Sources_Part1_AMReX



SUBROUTINE Poseidon_Input_Sources_XCFC_AMReX_Caller( A_Difference )
INTEGER     :: A_Difference
END SUBROUTINE Poseidon_Input_Sources_XCFC_AMReX_Caller

SUBROUTINE Poseidon_Input_Sources_Part1_AMReX_Caller( A_Difference )
INTEGER     :: A_Difference
END SUBROUTINE Poseidon_Input_Sources_Part1_AMReX_Caller

#endif




END MODULE Source_Input_XCFC_AMReX_Module