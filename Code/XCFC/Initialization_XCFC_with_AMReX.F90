   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Initialization_XCFC_with_AMReX_Module                                 !##!
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



USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Numbers_Module, &
            ONLY :  pi


USE Poseidon_Parameters, &
            ONLY :  Domain_Dim,             &
                    DEGREE,                 &
                    L_LIMIT,                &
                    Num_CFA_Eqs,            &
                    Num_CFA_Vars,           &
                    Verbose_Flag

USE Variables_Derived, &
            ONLY :  Num_R_Nodes,            &
                    Num_R_Nodesp1,          &
                    Var_Dim,                &
                    Elem_Var_Dim,           &
                    Block_Var_Dim,          &
                    LM_Length,              &
                    ULM_Length,             &
                    Prob_Dim,               &
                    Elem_Prob_Dim,          &
                    Elem_Prob_Dim_Sqr,      &
                    Block_Prob_Dim,         &
                    Num_Off_Diagonals,      &
                    Beta_Prob_Dim,          &
                    Beta_Elem_Prob_Dim

USE Variables_Mesh, &
            ONLY :  Num_R_Elements

USE Variables_Functions, &
            ONLY :  LM_Location,                    &
                    Calc_3D_Values_At_Location,     &
                    Calc_1D_CFA_Values,             &
                    Calc_Var_At_Loc_A,              &
                    Calc_Var_At_Loc_B

USE Variables_FP, &
            ONLY :  CFA_EQ_Flags,               &
                    CFA_EQ_Map,                 &
                    CFA_Var_Map,                &
                    CFA_Mat_Map,                &
                    Laplace_NNZ,                &
                    Beta_Diagonals,             &
                    Beta_Bandwidth,             &
                    Num_Matrices

USE Allocation_XCFC, &
            ONLY :  Allocate_XCFC

USE Allocation_Mesh, &
            ONLY :  Allocate_Mesh

USE FP_Functions_Results,   &
            ONLY :  Calc_FP_Values_At_Location, &
                    Calc_1D_CFA_Values_FP
                

USE FP_Intialize_Matrices, &
            ONLY :  Initialize_FP_Matrices

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_XCFC_Initialization,      &
                    Timer_XCFC_Matrix_Init

USE Variables_AMReX_Core, &
            ONLY :  MF_Source,          &
                    AMReX_Num_Levels,   &
                    iNumLeafElements

USE Poseidon_AMReX_BoxArraySize_Module, &
            ONLY : AMReX_BoxArraySize

USE Initialization_AMReX, &
            ONLY :  Initialize_AMReX_Maps

USE Initialization_Mesh_AMReX_Module, &
            ONLY :  Determine_AMReX_Mesh

USE Initialization_XCFC, &
            ONLY :  Create_Eq_Maps

USE IO_Setup_Report_Module, &
            ONLY :  PRINT_AMReX_Setup

USE FP_Functions_Results, &
            ONLY :  Calc_Var_At_Location_Type_A,    &
                    Calc_Var_At_Location_Type_B

IMPLICIT NONE


CONTAINS



 !+101+############################################################!
!                                                                   !
!              Initialization_XCFC_with_AMReX                       !
!                                                                   !
 !#################################################################!
SUBROUTINE  Initialization_XCFC_with_AMReX( CFA_EQ_Flags_Input )


INTEGER, DIMENSION(5), INTENT(IN), OPTIONAL             ::  CFA_EQ_Flags_Input



Integer                 :: Level
type(amrex_boxarray)    :: BA
type(amrex_box)         :: BX
type(c_ptr)             :: BX_ptr
INTEGER(c_int)          :: i, j

TYPE(amrex_mfiter)                              ::  mfi
TYPE(amrex_box)                                 ::  Box

TYPE(amrex_imultifab)                           ::  Level_Mask

INTEGER                                         ::  re, te, pe
INTEGER, DIMENSION(3)                           ::  iE
INTEGER, DIMENSION(3)                           ::  EHi, ELo
INTEGER                                         ::  nComp
INTEGER                                         ::  lvl

INTEGER                 :: lo(3)
INTEGER                 :: hi(3)
INTEGER                 :: Num_Boxes

! Determine Radial Base Variables from Multifab
PRINT*,"In Initialization_XCFC_with_AMReX"

IF ( PRESENT(CFA_EQ_Flags_Input) ) THEN
    CFA_EQ_Flags = CFA_EQ_Flags_Input
ELSE
    CFA_EQ_Flags = [1,1,1,0,0]
END IF

NUM_CFA_Eqs = SUM(CFA_EQ_Flags)
CALL Create_Eq_Maps()


CALL Initialize_AMReX_Maps()


CALL PRINT_AMReX_Setup()


Num_R_Elements = iNumLeafElements



! Determine Derived Varaibles
Num_R_Nodes         = DEGREE*NUM_R_ELEMENTS + 1
Num_R_Nodesp1       = Num_R_Nodes + 1


VAR_DIM             = LM_LENGTH*NUM_R_NODES
PROB_DIM            = NUM_CFA_VARS*VAR_DIM
Beta_Prob_Dim       = 3*Var_Dim
Beta_Diagonals = Beta_Elem_Prob_Dim
Beta_Bandwidth = 2*Beta_Diagonals+1

Laplace_NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1


CALL Allocate_Mesh()
CALL Determine_AMReX_Mesh()

! Allocate Arrays
CALL Allocate_XCFC()

! Construct Matrices
CALL TimerStart( Timer_XCFC_Matrix_Init )
CALL Initialize_FP_Matrices()
CALL TimerStop( Timer_XCFC_Matrix_Init )

Calc_3D_Values_At_Location  => Calc_FP_Values_At_Location
Calc_1D_CFA_Values          => Calc_1D_CFA_Values_FP

Calc_Var_At_Loc_A => Calc_Var_At_Location_Type_A
Calc_Var_At_Loc_B => Calc_Var_At_Location_Type_B

CALL TimerStop( Timer_XCFC_Initialization )







END SUBROUTINE  Initialization_XCFC_with_AMReX

#endif



END MODULE Initialization_XCFC_with_AMReX_Module
