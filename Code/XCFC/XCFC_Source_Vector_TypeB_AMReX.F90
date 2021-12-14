   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE XCFC_Source_Vector_TypeB_AMReX_Module                                 !##!
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

USE Poseidon_Numbers_Module, &
           ONLY :  pi,                         &
                   TwoPi

USE Poseidon_Units_Module, &
           ONLY :  GR_Source_Scalar

USE Poseidon_Parameters, &
           ONLY :  DEGREE,                     &
                   L_LIMIT,                    &
                   NUM_CFA_EQs,                &
                   NUM_CFA_VARs

USE Parameters_Variable_Indices, &
           ONLY :  iU_CF,                        &
                   iU_LF,                        &
                   iU_S1,                        &
                   iU_S2,                        &
                   iU_S3,                        &
                   iU_X1,                        &
                   iU_X2,                        &
                   iU_X3,                       &
                   iVB_S,                       &
                   iVB_X


USE Variables_Quadrature, &
           ONLY :  NUM_R_QUAD_POINTS,          &
                   NUM_T_QUAD_POINTS,          &
                   NUM_P_QUAD_POINTS,          &
                   NUM_TP_QUAD_POINTS,         &
                   INT_R_LOCATIONS,            &
                   INT_T_LOCATIONS,            &
                   INT_P_LOCATIONS

USE Variables_Mesh, &
           ONLY :  rlocs,                      &
                   tlocs,                      &
                   plocs
                

USE Variables_Tables, &
            ONLY :  Ylm_CC_Values,              &
                    Lagrange_Poly_Table

USE Variables_Derived, &
            ONLY :  LM_LENGTH

USE Variables_FP, &
            ONLY :  FP_Coeff_Vector_A,           &
                    FP_Coeff_Vector_B,           &
                    FP_Source_Vector_B

USE Functions_Jacobian, &
            ONLY :  Calc_Ahat

USE Poseidon_IO_Module, &
            ONLY :  Clock_In

USE FP_Functions_Mapping, &
            ONLY :  FP_FEM_Node_Map,            &
                    FP_tpd_Map,                 &
                    FP_Array_Map_TypeB

USE XCFC_Functions_Calc_Values_Module, &
            ONLY :  Calc_Int_Weights,               &
                    Calc_Val_On_Elem_TypeB,         &
                    Calc_Val_And_Drv_On_Elem_TypeA, &
                    Calc_Val_And_Drv_On_Elem_TypeB

USE XCFC_Functions_Physical_Source_Module, &
            ONLY :  Get_Physical_Source

USE XCFC_Source_Variables_Module, &
            ONLY :  Cur_R_Locs,         &
                    Cur_T_Locs,         &
                    R_Square,           &
                    Sin_Square,         &
                    RSin_Square,        &
                    TP_Sin_Val,         &
                    TP_Sin_Square,      &
                    TP_Cotan_Val,       &
                    TP_RSin_Square,     &
                    R_Int_Weights,      &
                    TP_Int_Weights,     &
                    Cur_Val_Psi,        &
                    Cur_Drv_Psi,        &
                    Cur_Val_AlphaPsi,   &
                    Cur_Drv_AlphaPsi,   &
                    Cur_Val_X,          &
                    Cur_Drv_X,          &
                    SourceTerm
                    
USE Variables_MPI, &
            ONLY :  myID_Poseidon,      &
                    nPROCS_Poseidon,    &
                    Poseidon_Comm_World

USE Poseidon_MPI_Utilities_Module, &
            ONLY :  STOP_MPI,               &
                    MPI_Master_Print,       &
                    MPI_All_Print


#ifdef POSEIDON_AMREX_FLAG
use amrex_base_module

USE amrex_box_module,   &
            ONLY:   amrex_box

USE amrex_boxarray_module, &
            ONLY:   amrex_boxarray


USE amrex_multifab_module,  &
            ONLY:   amrex_multifab,         &
                    amrex_multifab_build

USE Variables_AMReX_Multifabs, &
            ONLY :  MF_Source,  &
                    BA_Source,  &
                    DM_Source,  &
                    GM_Source,  &
                    nLevels,    &
                    Level_Ratio

USE Variables_AMReX_Source, &
            ONLY :  Source_PTR,         &
                    Mask_PTR,           &
                    iCoarse,            &
                    iFine

USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask

USE Poseidon_AMReX_Multilayer_Utilities_Module, &
            ONLY :  Find_Coarsest_Parent

#endif

USE MPI


IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!          XCFC_Calc_Source_Vector_TypeB_AMReX                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Calc_Source_Vector_TypeB_AMReX( iU, iVB )

INTEGER, INTENT(IN), DIMENSION(3)               ::  iU
INTEGER, INTENT(IN)                             ::  iVB


INTEGER                                         ::  lvl



FP_Source_Vector_B(:,iVB) = 0.0_idp
DO lvl = 0,nLevels-1


END DO ! lvl



END SUBROUTINE XCFC_Calc_Source_Vector_TypeB_AMReX










!+201+##########################################################################!
!                                                                               !
!          XCFC_Calc_Source_Vector_On_AMReX_Box_TypeB                           !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Calc_Source_Vector_On_AMReX_Box_TypeB(Level)

INTEGER, INTENT(IN)                             ::  Level

#ifdef POSEIDON_AMREX_FLAG

TYPE(amrex_mfiter)                              ::  mfi
TYPE(amrex_box)                                 ::  Box

TYPE(amrex_imultifab)                           ::  Level_Mask

INTEGER                                         ::  re, te, pe
INTEGER, DIMENSION(3)                           ::  iE
INTEGER, DIMENSION(3)                           ::  iEL, iEU
INTEGER                                         ::  nComp







!
!   MakeFineMask
!
IF ( Level < nLevels-1 ) THEN
    CALL AMReX_MakeFineMask(  Level_Mask,               &
                              MF_Source(Level)%ba,      &
                              MF_Source(Level)%dm,      &
                              MF_Source(Level+1)%ba,    &
                              iCoarse, iFine            )
ELSE
    ! Create Level_Mask all equal to 1
    CALL amrex_imultifab_build( Level_Mask,             &
                                MF_Source(Level)%ba,    &
                                MF_Source(Level)%dm,    &
                                1,                      &
                                0                       )
    CALL Level_Mask%SetVal(iCoarse)
END IF





!
!   Build mfiter
!
CALL amrex_mfiter_build(mfi, MF_Source(Level), tiling = .false. )
DO WHILE(mfi%next())
    Source_PTR => MF_Source(Level)%dataPtr(mfi)
    Mask_PTR   => Level_Mask%dataPtr(mfi)
    Box = mfi%tilebox()

    nComp =  MF_Source(Level)%ncomp()

    iEL = Box%lo
    iEU = Box%hi



    ! Create Tables
    Call XCFC_Create_Box_Tables( iEL, iEU)

    

    ! Move Through Box
    DO re = iEL(1),iEU(1)
    DO te = iEL(2),iEU(2)
    DO pe = iEL(3),iEU(3)
        IF ( Mask_PTR(RE,TE,PE,1) == 1 ) THEN
            iE = [re,te,pe]
!            CALL XCFC_Calc_Source_Vector_On_Element_TypeB( iU, iVB, iE, Level )
        END IF
    END DO ! pe
    END DO ! te
    END DO ! re

END DO
CALL amrex_mfiter_destroy(mfi)



#endif

END SUBROUTINE XCFC_Calc_Source_Vector_On_AMReX_Box_TypeB




END MODULE XCFC_Source_Vector_TypeB_AMReX_Module
