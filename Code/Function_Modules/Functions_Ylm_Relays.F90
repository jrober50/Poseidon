   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Functions_Ylm_Relays_Module                                           !##!
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
            ONLY :  idp


USE Variables_Tables, &
            ONLY :  Ylm_Values,                 &
                    Ylm_dt_Values,              &
                    Ylm_dp_Values,              &
                    Ylm_CC_Values,              &
                    Ylm_Elem_Values,            &
                    Ylm_Elem_dt_Values,         &
                    Ylm_Elem_dp_Values,         &
                    Ylm_Elem_CC_Values





IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!     Ylm_Value_Relay                                              				!
!                                                                               !
!###############################################################################!
PURE ELEMENTAL FUNCTION Ylm_Value_Relay( lm, tpd, iTE, iPE )

INTEGER, INTENT(IN)                         ::  lm
INTEGER, INTENT(IN)                         ::  tpd
INTEGER, INTENT(IN)                         ::  iTE
INTEGER, INTENT(IN)                         ::  iPE

COMPLEX(idp)                                ::  Ylm_Value_Relay

#ifdef POSEIDON_AMREX_FLAG

    
    Ylm_Value_Relay = Ylm_Elem_Values( lm, tpd )

#else

    Ylm_Value_Relay = Ylm_Values( lm, tpd, iTE, iPE)

#endif

END FUNCTION Ylm_Value_Relay





!+102+##########################################################################!
!                                                                               !
!     Ylm_dt_Relay                                                              !
!                                                                               !
!###############################################################################!
PURE ELEMENTAL FUNCTION Ylm_dt_Relay( lm, tpd, iTE, iPE )

INTEGER, INTENT(IN)                         ::  lm
INTEGER, INTENT(IN)                         ::  tpd
INTEGER, INTENT(IN)                         ::  iTE
INTEGER, INTENT(IN)                         ::  iPE

COMPLEX(idp)                                ::  Ylm_dt_Relay

#ifdef POSEIDON_AMREX_FLAG

    
    Ylm_dt_Relay = Ylm_Elem_dt_Values( lm, tpd )

#else

    Ylm_dt_Relay = Ylm_dt_Values( lm, tpd, iTE, iPE)

#endif

END FUNCTION Ylm_dt_Relay




!+103+##########################################################################!
!                                                                               !
!     Ylm_dp_Relay                                                              !
!                                                                               !
!###############################################################################!
PURE ELEMENTAL FUNCTION Ylm_dp_Relay( lm, tpd, iTE, iPE )

INTEGER, INTENT(IN)                         ::  lm
INTEGER, INTENT(IN)                         ::  tpd
INTEGER, INTENT(IN)                         ::  iTE
INTEGER, INTENT(IN)                         ::  iPE

COMPLEX(idp)                                ::  Ylm_dp_Relay

#ifdef POSEIDON_AMREX_FLAG

    
    Ylm_dp_Relay = Ylm_Elem_dp_Values( lm, tpd )

#else

    Ylm_dp_Relay = Ylm_dp_Values( lm, tpd, iTE, iPE)

#endif

END FUNCTION Ylm_dp_Relay





!+201+##########################################################################!
!                                                                               !
!     Ylm_CC_Relay                                                              !
!                                                                               !
!###############################################################################!
PURE ELEMENTAL FUNCTION Ylm_CC_Relay( lm, tpd, iTE, iPE )

INTEGER, INTENT(IN)                         ::  lm
INTEGER, INTENT(IN)                         ::  tpd
INTEGER, INTENT(IN)                         ::  iTE
INTEGER, INTENT(IN)                         ::  iPE

COMPLEX(idp)                                ::  Ylm_CC_Relay

#ifdef POSEIDON_AMREX_FLAG

    
    Ylm_CC_Relay = Ylm_Elem_CC_Values( lm, tpd )

#else

    Ylm_CC_Relay = Ylm_CC_Values( lm, tpd, iTE, iPE)

#endif

END FUNCTION Ylm_CC_Relay





END MODULE Functions_Ylm_Relays_Module
