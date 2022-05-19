   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Return_Routines_LF                                           !##!
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
!                                   Poseidon_XCFC_Run_Part1(), to calculate the value       !
!                                   of the Conformal Factor at the desired locations.       !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_LF_Native(   NE, NQ,                 &
                                        RQ_Input,               &
                                        TQ_Input,               &
                                        PQ_Input,               &
                                        Left_Limit,             &
                                        Right_Limit,            &
                                        Return_Lapse        )


INTEGER,    DIMENSION(3),                                   INTENT(IN)  ::  NE, NQ
REAL(idp),  DIMENSION(NQ(1)),                               INTENT(IN)  ::  RQ_Input
REAL(idp),  DIMENSION(NQ(2)),                               INTENT(IN)  ::  TQ_Input
REAL(idp),  DIMENSION(NQ(3)),                               INTENT(IN)  ::  PQ_Input
REAL(idp),                                                  INTENT(IN)  ::  Left_Limit
REAL(idp),                                                  INTENT(IN)  ::  Right_Limit

REAL(idp),  DIMENSION(NQ(1)*NQ(2)*NQ(3),NE(1),NE(2),NE(3)), INTENT(OUT) ::  Return_Lapse

REAL(idp),  DIMENSION(NQ(1)*NQ(2)*NQ(3),NE(1),NE(2),NE(3))              ::  Tmp_CF



CALL Poseidon_Return_Native_Type_A( iU_LF,                  &
                                    NE,                     &
                                    NQ,                     &
                                    RQ_Input,               &
                                    TQ_Input,               &
                                    PQ_Input,               &
                                    Left_Limit,             &
                                    Right_Limit,            &
                                    Return_Lapse            )

CALL Poseidon_Return_Native_Type_A( iU_CF,                  &
                                    NE,                     &
                                    NQ,                     &
                                    RQ_Input,               &
                                    TQ_Input,               &
                                    PQ_Input,               &
                                    Left_Limit,             &
                                    Right_Limit,            &
                                    Tmp_CF                  )



Return_Lapse = Return_Lapse/Tmp_CF



END SUBROUTINE Poseidon_Return_LF_Native



!+101+######################################################################################!
!                                                                                           !
!       Poseidon_Return_Lapse - Uses the results from the routine,                      !
!                                   Poseidon_XCFC_Run_Part1(), to calculate the value       !
!                                   of the Conformal Factor at the desired locations.       !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_LF_Native_Caller( Return_Lapse)



REAL(idp),  DIMENSION(  Caller_Quad_DOF,                     &
                        Num_R_Elements,                 &
                        Num_T_Elements,                 &
                        Num_P_Elements),  INTENT(OUT)   ::  Return_Lapse

INTEGER, DIMENSION(3)                                   ::  NE

REAL(idp),  DIMENSION(  Caller_Quad_DOF,                     &
                        Num_R_Elements,                 &
                        Num_T_Elements,                 &
                        Num_P_Elements)                 ::  Tmp_CF


NE = [ Num_R_Elements, Num_T_Elements, Num_P_Elements ]

CALL Poseidon_Return_Native_Type_A( iU_LF,                  &
                                    NE,                     &
                                    Caller_NQ,              &
                                    Caller_RQ_xlocs,        &
                                    Caller_TQ_xlocs,        &
                                    Caller_PQ_xlocs,        &
                                    Caller_xL(1),           &
                                    Caller_xL(2),           &
                                    Return_Lapse            )


CALL Poseidon_Return_Native_Type_A( iU_Cf,                  &
                                    NE,                     &
                                    Caller_NQ,              &
                                    Caller_RQ_xlocs,        &
                                    Caller_TQ_xlocs,        &
                                    Caller_PQ_xlocs,        &
                                    Caller_xL(1),           &
                                    Caller_xL(2),           &
                                    Tmp_CF                  )


Return_Lapse = Return_Lapse/Tmp_CF

END SUBROUTINE Poseidon_Return_LF_Native_Caller




#ifdef POSEIDON_AMREX_FLAG
!+101+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor_AMReX                                                     !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the desired locations, and fill an       !
!             AMReX multifab with the data.                                                 !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_LF_AMReX( NQ,                     &
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


iU = iU_LF

CALL Poseidon_Return_AMReX_Type_A(  iU,                     &
                                    NQ,                     &
                                    RQ_Input,               &
                                    TQ_Input,               &
                                    PQ_Input,               &
                                    Left_Limit,             &
                                    Right_Limit,            &
                                    nLevels,                &
                                    MF_Results              )



END SUBROUTINE Poseidon_Return_LF_AMReX



!+102+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor_AMReX_Caller                                              !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the locations set during initialization, !
!             and fill an AMReX multifab with the data.                                     !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_LF_AMReX_Caller( MF_Results )

TYPE(amrex_multifab),   INTENT(INOUT)           ::  MF_Results(0:Caller_nLevels-1)



INTEGER                                                         ::  iU


iU = iU_LF

CALL Poseidon_Return_AMReX_Type_A(  iU,                     &
                                    Caller_NQ,                      &
                                    Caller_RQ_xlocs,                &
                                    Caller_TQ_xlocs,                &
                                    Caller_PQ_xlocs,                &
                                    Caller_xL(1),                   &
                                    Caller_xL(2),                   &
                                    Caller_nLevels,                 &
                                    MF_Results              )

END SUBROUTINE Poseidon_Return_LF_AMReX_Caller


#else
!+101+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor_AMReX                                                     !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the desired locations, and fill an       !
!             AMReX multifab with the data.                                                 !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_LF_AMReX( )
END SUBROUTINE Poseidon_Return_LF_AMReX



!+102+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor_AMReX_Caller                                              !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the locations set during initialization, !
!             and fill an AMReX multifab with the data.                                     !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_LF_AMReX_Caller(A_Difference )
INTEGER, INTENT(IN)         :: A_Difference
END SUBROUTINE Poseidon_Return_LF_AMReX_Caller



#endif










END MODULE Poseidon_Return_Routines_LF
