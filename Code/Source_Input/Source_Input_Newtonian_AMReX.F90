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
                    Caller_PQ_xlocs

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

USE Flags_Source_Input_Module, &
            ONLY :  lPF_SI_Flags,       &
                    iPF_SI_MF_Ready


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

INTEGER                                             ::  Their_DOF



IF ( Verbose_Flag ) CALL Run_Message('Receiving Newtonian Sources. Container : AMReX Multifab.')
CALL TimerStart(Timer_Poisson_SourceInput)



! Define Interpolation Matrix
Their_DOF = Input_NQ(1)*Input_NQ(2)*Input_NQ(3)





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







END SUBROUTINE Poseidon_Input_Sources_Newtonian_AMReX_Caller






#else


SUBROUTINE Poseidon_Input_Sources_Newtonian_AMReX( )
END SUBROUTINE Poseidon_Input_Sources_Newtonian_AMReX


SUBROUTINE Poseidon_Input_Sources_Newtonian_AMReX_Caller( A_Difference )
INTEGER     :: A_Difference
END SUBROUTINE Poseidon_Input_Sources_Newtonian_AMReX_Caller



#endif



END MODULE Source_Input_Newtonian_AMReX_Module

