   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Return_Routines_CF                                           !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
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
            ONLY : idp

USE Poseidon_Parameters, &
            ONLY :  DEGREE

USE Parameters_Variable_Indices, &
            ONLY :  iU_CF

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    Num_T_Elements,         &
                    Num_P_Elements


USE Variables_Interface, &
            ONLY :  Caller_NQ,                      &
                    Caller_Quad_DOF,                &
                    Caller_xL,                      &
                    Caller_RQ_xlocs,                &
                    Caller_TQ_xlocs,                &
                    Caller_PQ_xlocs

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Num_Levels


#ifdef POSEIDON_AMREX_FLAG
use amrex_base_module

USE amrex_multifab_module,  &
            ONLY :  amrex_multifab

#endif

USE Return_Functions_Native, &
            ONLY :  Poseidon_Return_Native_Type_A,  &
                    Poseidon_Return_Native_Type_B

USE Return_Functions_AMReX, &
            ONLY :  Poseidon_Return_AMReX_Type_A,  &
                    Poseidon_Return_AMReX_Type_B


IMPLICIT NONE


CONTAINS
!+101+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor - Uses the results from the routine,                      !
!                                   Poseidon_XCFC_Run_Part1(), to calculate the value       !
!                                   of the Conformal Factor at the desired locations.       !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_CF_Native(   NE, NQ,                 &
                                        RQ_Input,               &
                                        TQ_Input,               &
                                        PQ_Input,               &
                                        Left_Limit,             &
                                        Right_Limit,            &
                                        Return_ConFactor        )


INTEGER,    DIMENSION(3),                                   INTENT(IN)  ::  NE, NQ
REAL(idp),  DIMENSION(NQ(1)),                               INTENT(IN)  ::  RQ_Input
REAL(idp),  DIMENSION(NQ(2)),                               INTENT(IN)  ::  TQ_Input
REAL(idp),  DIMENSION(NQ(3)),                               INTENT(IN)  ::  PQ_Input
REAL(idp),                                                  INTENT(IN)  ::  Left_Limit
REAL(idp),                                                  INTENT(IN)  ::  Right_Limit

REAL(idp),  DIMENSION(NQ(1)*NQ(2)*NQ(3),NE(1),NE(2),NE(3)), INTENT(OUT) ::  Return_Confactor

 INTEGER                                                                ::  iU

iU = iU_CF

CALL Poseidon_Return_Native_Type_A( iU,                     &
                                    NE,                     &
                                    NQ,                     &
                                    RQ_Input,               &
                                    TQ_Input,               &
                                    PQ_Input,               &
                                    Left_Limit,             &
                                    Right_Limit,            &
                                    Return_Confactor        )


END SUBROUTINE Poseidon_Return_CF_Native



!+101+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor - Uses the results from the routine,                      !
!                                   Poseidon_XCFC_Run_Part1(), to calculate the value       !
!                                   of the Conformal Factor at the desired locations.       !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_CF_Native_Caller( Return_ConFactor )


REAL(idp),  DIMENSION( Caller_Quad_DOF,                 &
                       Num_R_Elements,                  &
                       Num_T_Elements,                  &
                       Num_P_Elements),  INTENT(OUT)    ::  Return_Confactor

INTEGER                                                 ::  iU
INTEGER, DIMENSION(3)                                   ::  NE


iU = iU_CF
NE = [ Num_R_Elements, Num_T_Elements, Num_P_Elements ]

CALL Poseidon_Return_Native_Type_A( iU,                     &
                                    NE,                     &
                                    Caller_NQ,              &
                                    Caller_RQ_xlocs,        &
                                    Caller_TQ_xlocs,        &
                                    Caller_PQ_xlocs,        &
                                    Caller_xL(1),           &
                                    Caller_xL(2),           &
                                    Return_Confactor        )


END SUBROUTINE Poseidon_Return_CF_Native_Caller




#ifdef POSEIDON_AMREX_FLAG
!+101+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor_AMReX                                                     !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the desired locations, and fill an       !
!             AMReX multifab with the data.                                                 !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_CF_AMReX( NQ,                     &
                                    RQ_Input,               &
                                    TQ_Input,               &
                                    PQ_Input,               &
                                    Left_Limit,             &
                                    Right_Limit,            &
                                    nLevels,                &
                                    MF_Results              )


INTEGER,    DIMENSION(3),                       INTENT(IN)      ::  NQ
REAL(idp),  DIMENSION(NQ(1)),                   INTENT(IN)      ::  RQ_Input
REAL(idp),  DIMENSION(NQ(2)),                   INTENT(IN)      ::  TQ_Input
REAL(idp),  DIMENSION(NQ(3)),                   INTENT(IN)      ::  PQ_Input
REAL(idp),                                      INTENT(IN)      ::  Left_Limit
REAL(idp),                                      INTENT(IN)      ::  Right_Limit

INTEGER,                                        INTENT(IN)      ::  nLevels
TYPE(amrex_multifab),                           INTENT(INOUT)   ::  MF_Results(0:nLevels-1)



INTEGER                                                         ::  iU


iU = iU_CF

CALL Poseidon_Return_AMReX_Type_A(  iU,                     &
                                    NQ,                     &
                                    RQ_Input,               &
                                    TQ_Input,               &
                                    PQ_Input,               &
                                    Left_Limit,             &
                                    Right_Limit,            &
                                    nLevels,                &
                                    MF_Results              )

END SUBROUTINE Poseidon_Return_CF_AMReX



!+102+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor_AMReX_Caller                                              !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the locations set during initialization, !
!             and fill an AMReX multifab with the data.                                     !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_CF_AMReX_Caller( MF_Results )

TYPE(amrex_multifab),   INTENT(INOUT)           ::  MF_Results(0:AMReX_Num_Levels-1)



INTEGER                                                         ::  iU


iU = iU_CF

CALL Poseidon_Return_AMReX_Type_A(  iU,                     &
                                    Caller_NQ,              &
                                    Caller_RQ_xlocs,        &
                                    Caller_TQ_xlocs,        &
                                    Caller_PQ_xlocs,        &
                                    Caller_xL(1),           &
                                    Caller_xL(2),           &
                                    AMReX_Num_Levels,       &
                                    MF_Results              )

END SUBROUTINE Poseidon_Return_CF_AMReX_Caller

#else


!+101+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor_AMReX                                                     !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the desired locations, and fill an       !
!             AMReX multifab with the data.                                                 !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_CF_AMReX( )

END SUBROUTINE Poseidon_Return_CF_AMReX



!+102+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor_AMReX_Caller                                              !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the locations set during initialization, !
!             and fill an AMReX multifab with the data.                                     !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_CF_AMReX_Caller( A_Difference )

INTEGER, INTENT(IN)         :: A_Difference

END SUBROUTINE Poseidon_Return_CF_AMReX_Caller
#endif





END MODULE Poseidon_Return_Routines_CF
