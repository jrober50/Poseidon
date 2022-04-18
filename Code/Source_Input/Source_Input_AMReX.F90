  !################################################################################!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Source_Input_AMReX                                                           !##!
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
USE amrex_boxarray_module,  ONLY: &
  amrex_boxarray,         &
  amrex_boxarray_build,   &
  amrex_boxarray_destroy
USE amrex_distromap_module, ONLY: &
  amrex_distromap,       &
  amrex_distromap_build, &
  amrex_distromap_destroy
USE amrex_multifab_module,  ONLY: &
  amrex_multifab, &
  amrex_multifab_build, &
  amrex_mfiter, &
  amrex_mfiter_build, &
  amrex_mfiter_destroy
#endif

USE Poseidon_Internal_Communication_Module, &
            ONLY :  Poseidon_CFA_Block_Share

USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM,         &
                    Poseidon_Remesh_Flag, &
                    Verbose_Flag


USE Parameters_Variable_Indices, &
            ONLY :  iS_E,               &
                    iS_S,               &
                    iS_S1,              &
                    iS_S2,              &
                    iS_S3


USE Variables_Mesh, &
            ONLY :  NUM_R_ELEMENTS,     &
                    NUM_T_ELEMENTS,     &
                    NUM_P_ELEMENTS,     &
                    drlocs,             &
                    dtlocs,             &
                    dplocs

USE Variables_Source, &
            ONLY :  Block_Source_E,     &
                    Block_Source_S,     &
                    Block_Source_Si

USE Poseidon_IO_Module, &
            ONLY :  OUTPUT_POSEIDON_SOURCES_1D


USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Num_T_Quad_Points,      &
                    Num_P_Quad_Points,      &
                    Num_TP_Quad_Points,     &
                    Int_R_Locations,        &
                    Int_T_Locations,        &
                    Int_P_Locations,        &
                    Int_R_Weights,          &
                    Int_T_Weights,          &
                    Int_P_Weights,          &
                    Int_TP_Weights



USE Functions_Math, &
            ONLY :  Lagrange_Poly

#ifdef POSEIDON_AMREX_FLAG
USE Variables_AMReX_Core, &
            ONLY :  MF_Source,          &
                    AMReX_Num_Levels

USE Poseidon_AMReX_Utilities_Module,    &
            ONLY :  AMReX2Poseidon,     &
                    UnpackSources_AMReX


USE Initialization_XCFC_with_AMReX_Module, &
            ONLY :  Initialization_XCFC_with_AMReX


#endif


USE Timer_Routines_Module, &
            ONLY :  TimerStart,     &
                    TimerSTop

USE Timer_Variables_Module, &
            ONLY :  Timer_GR_SourceInput


use mpi



IMPLICIT NONE

CONTAINS


#ifdef POSEIDON_AMREX_FLAG
!+101+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources_AMREX( MF_SRC_Input,          &
                                         MF_Src_nComps,         &
                                         Num_Levels,            &
                                         Input_NQ,              &
                                         Input_R_Quad,          &
                                         Input_T_Quad,          &
                                         Input_P_Quad,          &
                                         Input_xL               )

TYPE(amrex_multifab),                   INTENT(IN)  ::  MF_SRC_Input(0:Num_Levels-1)
INTEGER,                                INTENT(IN)  ::  MF_Src_nComps
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

INTEGER                                             ::  Local_R
INTEGER                                             ::  Local_T
INTEGER                                             ::  Local_P

INTEGER                                             ::  Input_T
INTEGER                                             ::  Input_P

INTEGER                                             ::  nComp
REAL(idp), CONTIGUOUS, POINTER                      ::  My_PTR(:,:,:,:)
REAL(idp), CONTIGUOUS, POINTER                      ::  Their_PTR(:,:,:,:)

TYPE(amrex_mfiter)                                  ::  mfi
TYPE(amrex_box)                                     ::  Box

REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  Translation_Matrix
INTEGER                                             ::  My_DOF
INTEGER                                             ::  Their_DOF

REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  R_Lag_Poly_Values
REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  T_Lag_Poly_Values
REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  P_Lag_Poly_Values



IF (Verbose_Flag) THEN
    PRINT*,"In Poseidon_Input_Sources_AMREX"
END IF
CALL TimerStart(Timer_GR_SourceInput)

! Define Interpolation Matrix
My_DOF    = Num_R_Quad_Points*Num_T_Quad_Points*Num_P_Quad_Points
Their_DOF = Input_NQ(1)*Input_NQ(2)*Input_NQ(3)


ALLOCATE(Translation_Matrix(1:Their_DOF, 1:My_DOF))
ALLOCATE( R_Lag_Poly_Values(1:Input_NQ(1),1:NUM_R_QUAD_POINTS) )
ALLOCATE( T_Lag_Poly_Values(1:Input_NQ(2),1:NUM_T_QUAD_POINTS) )
ALLOCATE( P_Lag_Poly_Values(1:Input_NQ(3),1:NUM_P_QUAD_POINTS) )




DO Local_R = 1,NUM_R_QUAD_POINTS
    R_Lag_Poly_Values(:,Local_R) = Lagrange_Poly( Int_R_Locations(Local_R), &
                                                  Input_NQ(1)-1,            &
                                                  Input_R_Quad              )

END DO

DO Local_T = 1,NUM_T_QUAD_POINTS
    T_Lag_Poly_Values(:,Local_T) = Lagrange_Poly( Int_T_Locations(Local_T), &
                                                  Input_NQ(2)-1,            &
                                                  Input_T_Quad              )

END DO

DO Local_P = 1,NUM_P_QUAD_POINTS
    P_Lag_Poly_Values(:,Local_P) = Lagrange_Poly( Int_P_Locations(Local_P), &
                                                  Input_NQ(3)-1,            &
                                                  Input_P_Quad              )

END DO



DO Local_P = 1,NUM_P_QUAD_POINTS
DO Local_T = 1,NUM_T_QUAD_POINTS
DO Local_R = 1,NUM_R_QUAD_POINTS

    Local_Here = (Local_P-1) * NUM_T_QUAD_POINTS * NUM_R_QUAD_POINTS        &
               + (Local_T-1) * NUM_R_QUAD_POINTS                            &
               + Local_R

    DO Input_P = 1,Input_NQ(3)
    DO Input_T = 1,Input_NQ(2)

            Here = (Input_P-1) * Input_NQ(2) * Input_NQ(1)   &
                 + (Input_T-1) * Input_NQ(1)

            There = Here + Input_NQ(1)

            Translation_Matrix(Here+1:There, Local_Here)  =                 &
                              R_Lag_Poly_Values(1:Input_NQ(1),Local_R)    &
                            * T_Lag_Poly_Values(Input_T,Local_T)            &
                            * P_Lag_Poly_Values(Input_P,Local_P)

    END DO  !   Input_T Loop
    END DO  !   Input_P Loop
END DO  !   Local_R Loop
END DO  !   Local_T Loop
END DO  !   Local_P looop






DO level = 0,AMReX_Num_Levels-1

    CALL amrex_multifab_build(  MF_Source(level),           &
                                MF_Src_Input(Level)%BA,     &
                                MF_Src_Input(Level)%DM,     &
                                MF_Src_nComps, 1                        )



    CALL amrex_mfiter_build(mfi, MF_Source(level), tiling = .true. )

    DO WHILE(mfi%next())
        Their_PTR => MF_Src_Input(level)%dataPtr(mfi)
        My_PTR    => MF_Source(level)%dataPtr(mfi)



        Box = mfi%tilebox()
        nComp =  MF_Source(level)%ncomp()

        iEL = Box%lo
        iEU = Box%hi


        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)



        DO si = 1,2+DOMAIN_DIM
        DO Local_Here = 1,My_DOF

            Here  = (si-1)*Their_DOF+1
            There = si*Their_DOF

            Index = (si-1)*My_DOF+Local_Here

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



IF ( Poseidon_Remesh_Flag ) THEN
    Call Initialization_XCFC_with_AMReX()
END IF



END SUBROUTINE Poseidon_Input_Sources_AMREX









!+101+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources1_AMREX( MF_Src_Input,          &
                                          MF_Src_nComps,         &
                                          Num_Levels,            &
                                          Input_NQ,              &
                                          Input_R_Quad,          &
                                          Input_T_Quad,          &
                                          Input_P_Quad,          &
                                          Input_xL               )

TYPE(amrex_multifab),                   INTENT(IN)  ::  MF_Src_Input(0:Num_Levels-1)
INTEGER,                                INTENT(IN)  ::  MF_Src_nComps
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

INTEGER                                             ::  Local_R
INTEGER                                             ::  Local_T
INTEGER                                             ::  Local_P

INTEGER                                             ::  Input_R
INTEGER                                             ::  Input_T
INTEGER                                             ::  Input_P

INTEGER                                             ::  nComp
REAL(idp), CONTIGUOUS, POINTER                      ::  My_PTR(:,:,:,:)
REAL(idp), CONTIGUOUS, POINTER                      ::  Their_PTR(:,:,:,:)

TYPE(amrex_mfiter)                                  ::  mfi
TYPE(amrex_box)                                     ::  Box

REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  Translation_Matrix
INTEGER                                             ::  My_DOF
INTEGER                                             ::  Their_DOF

REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  R_Lag_Poly_Values
REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  T_Lag_Poly_Values
REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  P_Lag_Poly_Values


IF (Verbose_Flag) THEN
    PRINT*,"In Poseidon_Input_Sources1_AMREX"
END IF
CALL TimerStart(Timer_GR_SourceInput)

! Define Interpolation Matrix
My_DOF    = Num_R_Quad_Points*Num_T_Quad_Points*Num_P_Quad_Points
Their_DOF = Input_NQ(1)*Input_NQ(2)*Input_NQ(3)


ALLOCATE(Translation_Matrix(1:Their_DOF, 1:My_DOF))
ALLOCATE( R_Lag_Poly_Values(1:Input_NQ(1),1:NUM_R_QUAD_POINTS) )
ALLOCATE( T_Lag_Poly_Values(1:Input_NQ(2),1:NUM_T_QUAD_POINTS) )
ALLOCATE( P_Lag_Poly_Values(1:Input_NQ(3),1:NUM_P_QUAD_POINTS) )




DO Local_R = 1,NUM_R_QUAD_POINTS
    R_Lag_Poly_Values(:,Local_R) = Lagrange_Poly( Int_R_Locations(Local_R), &
                                                  Input_NQ(1)-1,            &
                                                  Input_R_Quad              )

END DO

DO Local_T = 1,NUM_T_QUAD_POINTS
    T_Lag_Poly_Values(:,Local_T) = Lagrange_Poly( Int_T_Locations(Local_T), &
                                                  Input_NQ(2)-1,            &
                                                  Input_T_Quad              )

END DO

DO Local_P = 1,NUM_P_QUAD_POINTS
    P_Lag_Poly_Values(:,Local_P) = Lagrange_Poly( Int_P_Locations(Local_P), &
                                                  Input_NQ(3)-1,            &
                                                  Input_P_Quad              )

END DO



DO Local_P = 1,NUM_P_QUAD_POINTS
DO Local_T = 1,NUM_T_QUAD_POINTS
DO Local_R = 1,NUM_R_QUAD_POINTS

    Local_Here = (Local_P-1) * NUM_T_QUAD_POINTS * NUM_R_QUAD_POINTS        &
               + (Local_T-1) * NUM_R_QUAD_POINTS                            &
               + Local_R

    DO Input_P = 1,Input_NQ(3)
    DO Input_T = 1,Input_NQ(2)

            Here = (Input_P-1) * Input_NQ(2) * Input_NQ(1)   &
                 + (Input_T-1) * Input_NQ(1)

            There = Here + Input_NQ(1)

            Translation_Matrix(Here+1:There, Local_Here)  =                 &
                              R_Lag_Poly_Values(1:Input_NQ(1),Local_R)    &
                            * T_Lag_Poly_Values(Input_T,Local_T)            &
                            * P_Lag_Poly_Values(Input_P,Local_P)

    END DO  !   Input_T Loop
    END DO  !   Input_P Loop
END DO  !   Local_R Loop
END DO  !   Local_T Loop
END DO  !   Local_P looop






DO level = 0,AMReX_Num_Levels-1

    CALL amrex_multifab_build(  MF_Source(level),           &
                                MF_Src_Input(Level)%BA,     &
                                MF_Src_Input(Level)%DM,     &
                                MF_Src_nComps, 1                        )



    CALL amrex_mfiter_build(mfi, MF_Source(level), tiling = .true. )

    DO WHILE(mfi%next())
        Their_PTR => MF_Src_Input(level)%dataPtr(mfi)
        My_PTR    => MF_Source(level)%dataPtr(mfi)



        Box = mfi%tilebox()
        nComp =  MF_Source(level)%ncomp()

        iEL = Box%lo
        iEU = Box%hi


        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)



        si = iS_E
        DO Local_Here = 1,My_DOF

            Here  = (si-1)*Their_DOF+1
            There = si*Their_DOF

            Index = (si-1)*My_DOF+Local_Here

            My_PTR(re,te,pe,Index) = DOT_PRODUCT( Translation_Matrix(:,Local_Here), &
                                                  Their_PTR(re,te,pe,Here:There)    )

        END DO  ! Local_Here
        


        DO si = iS_S1,iS_S2
        DO Local_Here = 1,My_DOF

            Here  = (si-1)*Their_DOF+1
            There = si*Their_DOF

            Index = (si-1)*My_DOF+Local_Here

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



IF ( Poseidon_Remesh_Flag ) THEN
    Call Initialization_XCFC_with_AMReX()
END IF



END SUBROUTINE Poseidon_Input_Sources1_AMREX




!+101+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources2_AMREX( MF_SRC_Input,          &
                                         MF_Src_nComps,         &
                                         Num_Levels,            &
                                         Input_NQ,              &
                                         Input_R_Quad,          &
                                         Input_T_Quad,          &
                                         Input_P_Quad,          &
                                         Input_xL               )

TYPE(amrex_multifab),                   INTENT(IN)  ::  MF_SRC_Input(0:Num_Levels-1)
INTEGER,                                INTENT(IN)  ::  MF_Src_nComps
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

INTEGER                                             ::  Local_R
INTEGER                                             ::  Local_T
INTEGER                                             ::  Local_P

INTEGER                                             ::  Input_T
INTEGER                                             ::  Input_P

INTEGER                                             ::  nComp
REAL(idp), CONTIGUOUS, POINTER                      ::  My_PTR(:,:,:,:)
REAL(idp), CONTIGUOUS, POINTER                      ::  Their_PTR(:,:,:,:)

TYPE(amrex_mfiter)                                  ::  mfi
TYPE(amrex_box)                                     ::  Box

REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  Translation_Matrix
INTEGER                                             ::  My_DOF
INTEGER                                             ::  Their_DOF

REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  R_Lag_Poly_Values
REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  T_Lag_Poly_Values
REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  P_Lag_Poly_Values



IF (Verbose_Flag) THEN
    PRINT*,"In Poseidon_Input_Sources_AMREX"
END IF
CALL TimerStart(Timer_GR_SourceInput)

! Define Interpolation Matrix
My_DOF    = Num_R_Quad_Points*Num_T_Quad_Points*Num_P_Quad_Points
Their_DOF = Input_NQ(1)*Input_NQ(2)*Input_NQ(3)


ALLOCATE(Translation_Matrix(1:Their_DOF, 1:My_DOF))
ALLOCATE( R_Lag_Poly_Values(1:Input_NQ(1),1:NUM_R_QUAD_POINTS) )
ALLOCATE( T_Lag_Poly_Values(1:Input_NQ(2),1:NUM_T_QUAD_POINTS) )
ALLOCATE( P_Lag_Poly_Values(1:Input_NQ(3),1:NUM_P_QUAD_POINTS) )




DO Local_R = 1,NUM_R_QUAD_POINTS
    R_Lag_Poly_Values(:,Local_R) = Lagrange_Poly( Int_R_Locations(Local_R), &
                                                  Input_NQ(1)-1,            &
                                                  Input_R_Quad              )

END DO

DO Local_T = 1,NUM_T_QUAD_POINTS
    T_Lag_Poly_Values(:,Local_T) = Lagrange_Poly( Int_T_Locations(Local_T), &
                                                  Input_NQ(2)-1,            &
                                                  Input_T_Quad              )

END DO

DO Local_P = 1,NUM_P_QUAD_POINTS
    P_Lag_Poly_Values(:,Local_P) = Lagrange_Poly( Int_P_Locations(Local_P), &
                                                  Input_NQ(3)-1,            &
                                                  Input_P_Quad              )

END DO



DO Local_P = 1,NUM_P_QUAD_POINTS
DO Local_T = 1,NUM_T_QUAD_POINTS
DO Local_R = 1,NUM_R_QUAD_POINTS

    Local_Here = (Local_P-1) * NUM_T_QUAD_POINTS * NUM_R_QUAD_POINTS        &
               + (Local_T-1) * NUM_R_QUAD_POINTS                            &
               + Local_R

    DO Input_P = 1,Input_NQ(3)
    DO Input_T = 1,Input_NQ(2)

            Here = (Input_P-1) * Input_NQ(2) * Input_NQ(1)   &
                 + (Input_T-1) * Input_NQ(1)

            There = Here + Input_NQ(1)

            Translation_Matrix(Here+1:There, Local_Here)  =                 &
                              R_Lag_Poly_Values(1:Input_NQ(1),Local_R)    &
                            * T_Lag_Poly_Values(Input_T,Local_T)            &
                            * P_Lag_Poly_Values(Input_P,Local_P)

    END DO  !   Input_T Loop
    END DO  !   Input_P Loop
END DO  !   Local_R Loop
END DO  !   Local_T Loop
END DO  !   Local_P looop






DO level = 0,AMReX_Num_Levels-1

    CALL amrex_multifab_build(  MF_Source(level),           &
                                MF_Src_Input(Level)%BA,     &
                                MF_Src_Input(Level)%DM,     &
                                MF_Src_nComps, 1                        )



    CALL amrex_mfiter_build(mfi, MF_Source(level), tiling = .true. )

    DO WHILE(mfi%next())
        Their_PTR => MF_Src_Input(level)%dataPtr(mfi)
        My_PTR    => MF_Source(level)%dataPtr(mfi)



        Box = mfi%tilebox()
        nComp =  MF_Source(level)%ncomp()

        iEL = Box%lo
        iEU = Box%hi


        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)



        DO si = 1,2+DOMAIN_DIM
        DO Local_Here = 1,My_DOF

            Here  = (si-1)*Their_DOF+1
            There = si*Their_DOF

            Index = (si-1)*My_DOF+Local_Here

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



IF ( Poseidon_Remesh_Flag ) THEN
    Call Initialization_XCFC_with_AMReX()
END IF



END SUBROUTINE Poseidon_Input_Sources2_AMREX





#else


SUBROUTINE Poseidon_Input_Sources_AMREX( )
END SUBROUTINE Poseidon_Input_Sources_AMREX

SUBROUTINE Poseidon_XCFC_Input_Sources_AMREX( )
END SUBROUTINE Poseidon_XCFC_Input_Sources_AMREX


#endif




END MODULE Source_Input_AMReX
