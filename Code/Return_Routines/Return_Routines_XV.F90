   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Return_Routines_XV                                           !##!
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

USE Parameters_Variable_Indices, &
            ONLY :  iVB_X,                      &
                    iU_X1,                      &
                    iU_X2,                      &
                    iU_X3

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    Num_T_Elements,         &
                    Num_P_Elements

USE Variables_Interface, &
            ONLY :  Caller_NQ,                      &
                    Caller_Quad_DOF,                     &
                    Caller_xL,                      &
                    Caller_RQ_xlocs,                &
                    Caller_TQ_xlocs,                &
                    Caller_PQ_xlocs


#ifdef POSEIDON_AMREX_FLAG
use amrex_base_module

USE amrex_box_module,   &
            ONLY :  amrex_box

USE amrex_boxarray_module, &
            ONLY :  amrex_boxarray

use amrex_fort_module, &
            ONLY :  amrex_spacedim

USE amrex_multifab_module,  &
            ONLY :  amrex_multifab,         &
                    amrex_multifab_build,   &
                    amrex_imultifab_build,  &
                    amrex_imultifab_destroy

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Num_Levels

USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask

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
!                                   Poseidon_XCFC_Run_Part2(), to calculate the value       !
!                                   of the Conformal Factor at the desired locations.       !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_XV_Native(   NE, NQ,                 &
                                    RQ_Input,           &
                                    TQ_Input,           &
                                    PQ_Input,           &
                                    Left_Limit,             &
                                    Right_Limit,            &
                                    Return_X            )

INTEGER,    DIMENSION(3),                                   INTENT(IN)  ::  NE, NQ
REAL(idp),  DIMENSION(NQ(1)),                               INTENT(IN)  ::  RQ_Input
REAL(idp),  DIMENSION(NQ(2)),                               INTENT(IN)  ::  TQ_Input
REAL(idp),  DIMENSION(NQ(3)),                               INTENT(IN)  ::  PQ_Input
REAL(idp),                                                  INTENT(IN)  ::  Left_Limit
REAL(idp),                                                  INTENT(IN)  ::  Right_Limit

REAL(idp),  DIMENSION(NQ(1)*NQ(2)*NQ(3),NE(1),NE(2),NE(3),1:3), INTENT(OUT) ::  Return_X

INTEGER                                                         ::  iU, iVB


iVB = iVB_X
DO iU = iU_X1, iU_X3

    CALL Poseidon_Return_Native_Type_B(  iU, iVB,                           &
                                        NE, NQ,                             &
                                        RQ_Input,                           &
                                        TQ_Input,                           &
                                        PQ_Input,                           &
                                        Left_Limit,                         &
                                        Right_Limit,                        &
                                        Return_X(:,:,:,:,iU-iU_X1+1)    )

END DO ! u Loop



END SUBROUTINE Poseidon_Return_XV_Native




!+102+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor - Uses the results from the routine,                      !
!                                   Poseidon_XCFC_Run_Part2(), to calculate the value       !
!                                   of the Conformal Factor at the desired locations.       !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_XV_Native_Caller( Return_X )


REAL(idp),  DIMENSION(  Caller_Quad_DOF,         &
                        Num_R_Elements,     &
                        Num_T_Elements,     &
                        Num_P_Elements,     &
                        1:3                 ),  INTENT(OUT)     ::  Return_X

INTEGER                                                                 ::  iU
INTEGER                                                                 ::  iVB
INTEGER, DIMENSION(3)                                                   ::  NE

NE = [ Num_R_Elements, Num_T_Elements, Num_P_Elements ]

iVB = iVB_X
DO iU = iU_X1, iU_X3

    CALL Poseidon_Return_Native_Type_B( iU, iVB,                            &
                                        NE,                                 &
                                        Caller_NQ,                          &
                                        Caller_RQ_xlocs,                    &
                                        Caller_TQ_xlocs,                    &
                                        Caller_PQ_xlocs,                    &
                                        Caller_xL(1),                       &
                                        Caller_xL(2),                       &
                                        Return_X(:,:,:,:,iU-iU_X1+1)    )

END DO ! u Loop



END SUBROUTINE Poseidon_Return_XV_Native_Caller





#ifdef POSEIDON_AMREX_FLAG
!+201+######################################################################################!
!                                                                                           !
!       Poseidon_Return_X_AMReX                                                     !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the desired locations, and fill an       !
!             AMReX multifab with the data.                                                 !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_XV_AMReX( NQ,                     &
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



INTEGER                                                     ::  iU
INTEGER                                                     ::  iVB

iVB = iVB_X

DO iU = iU_X1,iU_X3
    CALL Poseidon_Return_AMReX_Type_B(  iU,                     &
                                        iVB,                    &
                                        NQ,                     &
                                        RQ_Input,               &
                                        TQ_Input,               &
                                        PQ_Input,               &
                                        Left_Limit,             &
                                        Right_Limit,            &
                                        nLevels,                &
                                        MF_Results              )
END DO


END SUBROUTINE Poseidon_Return_XV_AMReX






!+302+######################################################################################!
!                                                                                           !
!       Poseidon_Return_X_AMReX_Caller                                                  !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the desired locations, and fill an       !
!             AMReX multifab with the data.                                                 !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_XV_AMReX_Caller( MF_Results )



TYPE(amrex_multifab),   INTENT(INOUT)               ::  MF_Results(0:AMReX_Num_Levels-1)



INTEGER                                             ::  iU
INTEGER                                             ::  iVB



iVB = iVB_X

DO iU = iU_X1,iU_X3
    CALL Poseidon_Return_AMReX_Type_B(  iU,                 &
                                        iVB,                &
                                        Caller_NQ,          &
                                        Caller_RQ_xlocs,    &
                                        Caller_TQ_xlocs,    &
                                        Caller_PQ_xlocs,    &
                                        Caller_xL(1),       &
                                        Caller_xL(2),       &
                                        AMReX_Num_Levels,   &
                                        MF_Results          )
END DO


END SUBROUTINE Poseidon_Return_XV_AMReX_Caller


#else

!+201+######################################################################################!
!                                                                                           !
!       Poseidon_Return_X_AMReX                                                     !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the desired locations, and fill an       !
!             AMReX multifab with the data.                                                 !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_XV_AMReX( )
END SUBROUTINE Poseidon_Return_XV_AMReX






!+302+######################################################################################!
!                                                                                           !
!       Poseidon_Return_X_AMReX_Caller                                                  !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the desired locations, and fill an       !
!             AMReX multifab with the data.                                                 !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_XV_AMReX_Caller(A_Difference )
INTEGER, INTENT(IN)         :: A_Difference
END SUBROUTINE Poseidon_Return_XV_AMReX_Caller

#endif



END MODULE Poseidon_Return_Routines_XV
