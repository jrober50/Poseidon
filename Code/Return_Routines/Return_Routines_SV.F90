   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Return_Routines_SV                                           !##!
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


USE Variables_Tables, &
            ONLY :  Ylm_Elem_Values,        &
                    Ylm_Elem_dt_Values,     &
                    Ylm_Elem_dp_Values,              &
                    Lagrange_Poly_Table,        &
                    Level_DX

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    Num_T_Elements,         &
                    Num_P_Elements

USE Variables_Derived, &
            ONLY :  LM_LENGTH

USE Variables_FP, &
            ONLY :  FP_Coeff_Vector_A,  &
                    FP_Coeff_Vector_B

USE Variables_Mesh, &
            ONLY :  rlocs,              &
                    tlocs,              &
                    plocs

USE Variables_AMReX_Source, &
            ONLY :  iLeaf,                &
                    iTrunk

USE Maps_Fixed_Point, &
            ONLY :  FP_Array_Map_TypeB

USE Maps_Quadrature, &
            ONLY :  Map_To_tpd

USE Maps_Domain, &
            ONLY :  Map_To_FEM_Node,        &
                    FEM_Elem_Map

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations

USE Functions_Math, &
            ONLY :  Lagrange_Poly

USE Initialization_Tables, &
            ONLY :  Initialize_Normed_Legendre_Tables_On_Level,     &
                    Initialize_Ylm_Tables_On_Elem


USE Variables_Interface, &
            ONLY :  Caller_nLevels,                 &
                    Caller_NQ,                      &
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
SUBROUTINE Poseidon_Return_SV_Native(   NE, NQ,                 &
                                    RQ_Input,           &
                                    TQ_Input,           &
                                    PQ_Input,           &
                                    Left_Limit,             &
                                    Right_Limit,            &
                                    Return_Shift            )

INTEGER,    DIMENSION(3),                                   INTENT(IN)  ::  NE, NQ
REAL(idp),  DIMENSION(NQ(1)),                               INTENT(IN)  ::  RQ_Input
REAL(idp),  DIMENSION(NQ(2)),                               INTENT(IN)  ::  TQ_Input
REAL(idp),  DIMENSION(NQ(3)),                               INTENT(IN)  ::  PQ_Input
REAL(idp),                                                  INTENT(IN)  ::  Left_Limit
REAL(idp),                                                  INTENT(IN)  ::  Right_Limit

REAL(idp),  DIMENSION(NQ(1)*NQ(2)*NQ(3),NE(1),NE(2),NE(3),1:3), INTENT(OUT) ::  Return_Shift


INTEGER                                                         ::  re, te, pe
INTEGER                                                         ::  rd, td, pd, tpd
INTEGER                                                         ::  d, lm, Here, iU, iVB

REAL(KIND = idp)                                                ::  Quad_Span
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  Local_Locations
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  LagP
REAL(KIND = idp), DIMENSION(1:NQ(1))                            ::  CUR_X_LOCS
COMPLEX(KIND = idp)                                             ::  TMP_U_Value
INTEGER                                                         ::  Current_Location

Quad_Span = Right_Limit - Left_Limit

Local_Locations = Initialize_LGL_Quadrature_Locations(DEGREE)
CUR_X_LOCS = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp

iVB = iVB_S
DO iU = iU_S1, iU_S3

    CALL Poseidon_Return_Native_Type_B(  iU, iVB,                           &
                                        NE, NQ,                             &
                                        RQ_Input,                           &
                                        TQ_Input,                           &
                                        PQ_Input,                           &
                                        Left_Limit,                         &
                                        Right_Limit,                        &
                                        Return_Shift(:,:,:,:,iU-iU_S1+1)    )

END DO ! u Loop



END SUBROUTINE Poseidon_Return_SV_Native




!+102+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor - Uses the results from the routine,                      !
!                                   Poseidon_XCFC_Run_Part2(), to calculate the value       !
!                                   of the Conformal Factor at the desired locations.       !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_SV_Native_Caller( Return_Shift )


REAL(idp),  DIMENSION(  Caller_Quad_DOF,         &
                        Num_R_Elements,     &
                        Num_T_Elements,     &
                        Num_P_Elements,     &
                        1:3                 ),  INTENT(OUT)     ::  Return_Shift

INTEGER                                                                 ::  iU
INTEGER                                                                 ::  iVB
INTEGER, DIMENSION(3)                                                   ::  NE

NE = [ Num_R_Elements, Num_T_Elements, Num_P_Elements ]

iVB = iVB_S
DO iU = iU_S1, iU_S3

    CALL Poseidon_Return_Native_Type_B( iU, iVB,                            &
                                        NE,                                 &
                                        Caller_NQ,                          &
                                        Caller_RQ_xlocs,                    &
                                        Caller_TQ_xlocs,                    &
                                        Caller_PQ_xlocs,                    &
                                        Caller_xL(1),                       &
                                        Caller_xL(2),                       &
                                        Return_Shift(:,:,:,:,iU-iU_S1+1)    )

END DO ! u Loop



END SUBROUTINE Poseidon_Return_SV_Native_Caller





#ifdef POSEIDON_AMREX_FLAG
!+201+######################################################################################!
!                                                                                           !
!       Poseidon_Return_Shift_AMReX                                                     !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the desired locations, and fill an       !
!             AMReX multifab with the data.                                                 !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_SV_AMReX( NQ,                     &
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

iVB = iVB_S

DO iU = iU_S1,iU_S3
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


END SUBROUTINE Poseidon_Return_SV_AMReX






!+302+######################################################################################!
!                                                                                           !
!       Poseidon_Return_Shift_AMReX_Caller                                                  !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the desired locations, and fill an       !
!             AMReX multifab with the data.                                                 !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_SV_AMReX_Caller( MF_Results )



TYPE(amrex_multifab),   INTENT(INOUT)               ::  MF_Results(0:Caller_nLevels-1)



INTEGER                                             ::  iU
INTEGER                                             ::  iVB



iVB = iVB_S

DO iU = iU_S1,iU_S3
    CALL Poseidon_Return_AMReX_Type_B(  iU,                 &
                                        iVB,                &
                                        Caller_NQ,          &
                                        Caller_RQ_xlocs,    &
                                        Caller_TQ_xlocs,    &
                                        Caller_PQ_xlocs,    &
                                        Caller_xL(1),       &
                                        Caller_xL(2),       &
                                        Caller_nLevels,     &
                                        MF_Results          )
END DO


END SUBROUTINE Poseidon_Return_SV_AMReX_Caller



#else
!+201+######################################################################################!
!                                                                                           !
!       Poseidon_Return_Shift_AMReX                                                     !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the desired locations, and fill an       !
!             AMReX multifab with the data.                                                 !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_SV_AMReX( )
END SUBROUTINE Poseidon_Return_SV_AMReX






!+302+######################################################################################!
!                                                                                           !
!       Poseidon_Return_Shift_AMReX_Caller                                                  !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the desired locations, and fill an       !
!             AMReX multifab with the data.                                                 !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_SV_AMReX_Caller( A_Difference )
INTEGER, INTENT(IN)         :: A_Difference
END SUBROUTINE Poseidon_Return_SV_AMReX_Caller


#endif


END MODULE Poseidon_Return_Routines_SV
