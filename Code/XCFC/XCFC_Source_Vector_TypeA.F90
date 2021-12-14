  !################################################################################!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE XCFC_Source_Vector_TypeA_Module                                              !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  !################################################################################!

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
           ONLY :  NUM_R_ELEMENTS,             &
                   NUM_T_ELEMENTS,             &
                   NUM_P_ELEMENTS,             &
                   rlocs,                      &
                   tlocs,                      &
                   plocs
                 
USE Variables_Source, &
           ONLY :  Block_Source_E,             &
                   Block_Source_S,             &
                   Block_Source_Si

USE Variables_Tables, &
           ONLY :   Ylm_CC_Values,              &
                    Ylm_Elem_CC_Values,         &
                    Level_dx,                   &
                    Lagrange_Poly_Table,        &
                    Lagpoly_MltiLayer_Table

USE Variables_Derived, &
           ONLY :  LM_LENGTH

USE Variables_FP, &
           ONLY :  FP_Coeff_Vector_A,            &
                   FP_Source_Vector_A

USE Functions_Jacobian, &
           ONLY :  Calc_Ahat

USE FP_Functions_Mapping, &
           ONLY :  FP_FEM_Node_Map,            &
                   FP_tpd_Map



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
                    Cur_Val_AlphaPsi,   &
                    Cur_Val_X,          &
                    Cur_Drv_X,          &
                    SourceTerm
                    
USE XCFC_Functions_Physical_Source_Module, &
            ONLY : Get_Physical_Source

USE XCFC_Functions_Calc_Values_Module, &
            ONLY :  Calc_Int_Weights,               &
                    Calc_Val_On_Elem_TypeA,         &
                    Calc_Val_And_Drv_On_Elem_TypeB

USE Initialization_Tables, &
            ONLY :  Initialize_Normed_Legendre_Tables_On_Level,     &
                    Initialize_Ylm_Tables_On_Elem


USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_XCFC_Lapse_SourceVector,  &
                    Timer_XCFC_ConFactor_SourceVector

USE MPI





#ifdef POSEIDON_AMREX_FLAG
use amrex_base_module

USE amrex_box_module,   &
            ONLY:   amrex_box

USE amrex_boxarray_module, &
            ONLY:   amrex_boxarray


USE amrex_multifab_module,  &
            ONLY:   amrex_multifab,         &
                    amrex_multifab_build,   &
                    amrex_imultifab_build

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






IMPLICIT NONE


CONTAINS

!+102+###########################################################################!
!                                                                                !
!           XCFC_Calc_Source_Vector_TypeA                                        !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Calc_Source_Vector_TypeA( iU, iEU, iEL )

INTEGER, INTENT(IN)                     ::  iU
INTEGER, INTENT(IN), DIMENSION(3)       ::  iEU
INTEGER, INTENT(IN), DIMENSION(3)       ::  iEL

INTEGER,             DIMENSION(3)       ::  iE
#ifndef POSEIDON_AMREX_FLAG
INTEGER                                 ::  re, te, pe
#endif


IF ( iU == iU_CF ) THEN
    CALL TimerStart( Timer_XCFC_ConFactor_SourceVector)
ELSEIF ( iU == iU_LF ) THEN
    CALL TimerStart( Timer_XCFC_Lapse_SourceVector)
END IF




#ifdef POSEIDON_AMREX_FLAG

    CALL XCFC_AMReX_Calc_Source_Vector_TypeA( iU )

#else
    
    FP_Source_Vector_A(:,:,iU) = 0.0_idp
    DO re = iEL(1),iEU(1)
    DO te = iEL(2),iEU(2)
    DO pe = iEL(3),iEU(3)
        iE = [re,te,pe]
        CALL XCFC_Calc_Source_Vector_On_Element_TypeA( iU, iE )
    END DO ! pe
    END DO ! te
    END DO ! re
#endif



IF ( iU == iU_CF ) THEN
    CALL TimerStop( Timer_XCFC_ConFactor_SourceVector)
ELSEIF ( iU == iU_LF ) THEN
    CALL TimerStop( Timer_XCFC_Lapse_SourceVector)
END IF



!PRINT*,FP_Source_Vector_A(:,:,iU)


END SUBROUTINE XCFC_Calc_Source_Vector_TypeA



!+102+###########################################################################!
!                                                                                !
!           XCFC_Calc_Source_Vector_TypeA                                        !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Calc_Source_Vector_On_Element_TypeA( iU, iE, Level_Option )

INTEGER, INTENT(IN)                             ::  iU
INTEGER, INTENT(IN), DIMENSION(3)               ::  iE
INTEGER, INTENT(IN), OPTIONAL                   ::  Level_Option


INTEGER                                         ::  rd, tpd, td, pd

REAL(KIND = idp)                                ::  DROT,     &
                                                    DTOT

INTEGER                                         ::  Level, i
INTEGER                                         ::  iCE(3)
INTEGER                                         ::  iRE(3)

IF (Present(Level_Option)) THEN
    Level = Level_Option
ELSE
    Level = 0
END IF


#ifdef POSEIDON_AMREX_FLAG


DO i = 1,3
    iCE(i) = Find_Coarsest_Parent(iE(i), Level)
    iRE(i) = 2.0_idp*MOD(iE(i),Level_Ratio(Level))
END DO


DROT = Level_dx(Level,1)/2.0_idp
DTOT = Level_dx(Level,2)/2.0_idp

CUR_R_LOCS(:) = DROT * (Int_R_Locations(:) + 1.0_idp + iE(1)*2.0_idp)
CUR_T_LOCS(:) = DTOT * (Int_T_Locations(:) + 1.0_idp + iE(2)*2.0_idp)

#else


DROT = 0.5_idp * (rlocs(iE(1)+1) - rlocs(iE(1)))
DTOT = 0.5_idp * (tlocs(iE(2)+1) - tlocs(iE(2)))

CUR_R_LOCS(:) = DROT * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(iE(1))
CUR_T_LOCS(:) = DTOT * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(iE(2))


#endif



R_SQUARE(:) = CUR_R_LOCS(:)*CUR_R_LOCS(:)
DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS
    tpd = FP_tpd_Map(td,pd)
    TP_Sin_Val(tpd)    = DSIN(CUR_T_LOCS(td))
    TP_Cotan_Val(tpd)  = 1.0_idp/DTAN(CUR_T_LOCS(td))
END DO
END DO
TP_Sin_Square(:) = TP_Sin_Val(:)*TP_Sin_Val


DO rd = 1,NUM_R_QUAD_POINTS
    TP_RSIN_SQUARE(:,rd) = R_SQUARE(rd)*TP_SIN_SQUARE(:)
END DO


!PRINT*,"Before Calc_Int_Weights",iE
CALL Calc_Int_Weights( DROT, DTOT,                  &
                       R_Square, TP_Sin_Val,        &
                       R_Int_Weights, TP_Int_Weights )

!PRINT*,"Before Calc_XCFC_CurVals_TypeA",iE
CALL Calc_XCFC_CurVals_TypeA( iE,       &
                              iU,       &
                              DROT,     &
                              DTOT,     &
                              Level     )

!PRINT*,"Before Create_XCFC_Vector_TypeA",iE
CALL Create_XCFC_Vector_TypeA( iE, iU, Level )



END SUBROUTINE XCFC_Calc_Source_Vector_On_Element_TypeA











!+201+##########################################################################!
!                                                                               !
!                  Create_XCFC_Vector_TypeA                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Create_XCFC_Vector_TypeA( iE, iU, Level )


INTEGER, INTENT(IN), DIMENSION(3)                           ::  iE
INTEGER, INTENT(IN)                                         ::  iU
INTEGER, INTENT(IN)                                         ::  Level


INTEGER                                                     ::  rd, d, lm_loc, td,pd,tpd
INTEGER                                                     ::  Current_i_Location

COMPLEX(KIND = idp)                                         ::  RHS_TMP
INTEGER                                                     ::  iCT


#ifdef POSEIDON_AMREX_FLAG

iCT = 2**(level+1) + mod(iE(1),2**level) - 2
DO lm_loc = 1,LM_LENGTH
DO d = 0,DEGREE
    
    RHS_TMP = 0.0_idp
    DO rd = 1,NUM_R_QUAD_POINTS

        RHS_TMP =  RHS_TMP                                          &
                 + SUM( SourceTerm( :, rd, iU )                     &
                       * Ylm_Elem_CC_Values( :, lm_loc )            &
                       * TP_Int_Weights(:)                     )    &
               * LagPoly_MltiLayer_Table( d, rd, 0, iCT)            &
               * R_Int_Weights(rd)

!        PRINT*,"________Source_Term_____________",lm_loc,d,rd
!        PRINT*,SourceTerm( :, rd, iU )
!        PRINT*,"_______Lagrange_Poly____________"
!        PRINT*,Ylm_CC_Values( :, lm_loc, iE(2), iE(3))
!        PRINT*,"_______TP_Int_Weights___________"
!        PRINT*,TP_INT_Weights(:)



!        PRINT*,lm_loc,d,rd,                                 &
!        SUM( SourceTerm( :, rd, iU )                        &
!                * Ylm_Elem_CC_Values( :, lm_loc )           &
!                * TP_Int_Weights(:)                     ),  &
!        LagPoly_MltiLayer_Table( d, rd, 0, iCT),            &
!        R_Int_Weights(rd)
        
    END DO  ! rd Loop
    


    Current_i_Location = FP_FEM_Node_Map(Find_Coarsest_Parent(iE(1),Level),d)
    FP_Source_Vector_A(Current_i_Location,lm_loc,iU)          &
        = FP_Source_Vector_A(Current_i_Location,lm_loc,iU)    &
        + RHS_TMP


!    PRINT*,Current_i_Location,FP_Source_Vector_A(Current_i_Location,lm_loc,iU)

END DO  ! d Loop
END DO  ! lm_loc Loop




#else

DO lm_loc = 1,LM_LENGTH
DO d = 0,DEGREE
    RHS_TMP = 0.0_idp
    DO rd = 1,NUM_R_QUAD_POINTS

        RHS_TMP =  RHS_TMP                                          &
                 + SUM( SourceTerm( :, rd, iU )                     &
                       * Ylm_CC_Values( :, lm_loc, iE(2), iE(3))    &
                       * TP_Int_Weights(:)                     )    &
               * Lagrange_Poly_Table( d, rd, 0)                     &
               * R_Int_Weights(rd)


!        PRINT*,"________Source_Term_____________",lm_loc,d,rd
!        PRINT*,SourceTerm( :, rd, iU )
!        PRINT*,"_______Lagrange_Poly____________"
!        PRINT*,Ylm_CC_Values( :, lm_loc, iE(2), iE(3))
!        PRINT*,"_______TP_Int_Weights___________"
!        PRINT*,TP_INT_Weights(:)

!        PRINT*,lm_loc,d,rd,                                 &
!        SUM( SourceTerm( :, rd, iU )                        &
!                * Ylm_CC_Values( :, lm_loc, iE(2), iE(3))   &
!                * TP_Int_Weights(:)                     ),  &
!        Lagrange_Poly_Table( d, rd, 0),                     &
!        R_Int_Weights(rd)

    END DO  ! rd Loop
    Current_i_Location = FP_FEM_Node_Map(iE(1),d)

    FP_Source_Vector_A(Current_i_Location,lm_loc,iU)          &
        = FP_Source_Vector_A(Current_i_Location,lm_loc,iU)    &
        + RHS_TMP

    
!    PRINT*,Current_i_Location,FP_Source_Vector_A(Current_i_Location,lm_loc,iU)

END DO  ! d Loop
END DO  ! lm_loc Loop


#endif




END SUBROUTINE Create_XCFC_Vector_TypeA










!+202+###########################################################################!
!                                                                                !
!                  Calc_Current_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_XCFC_CurVals_TypeA( iE, iU, DROT, DTOT, Level )

INTEGER, INTENT(IN), DIMENSION(3)                               ::  iE
INTEGER, INTENT(IN)                                             ::  iU
REAL(KIND = idp), INTENT(IN)                                    ::  DROT, DTOT
INTEGER, INTENT(IN)                                             ::  Level


REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3,3) ::  Ahat_Array
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3 )  ::  f
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points )    ::  AA_Array
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points)     ::  PhysSrc

INTEGER                                                         ::  tpd, rd
INTEGER                                                         ::  i, j

!PRINT*,"A",Level
CALL Calc_Int_Weights( DROT, DTOT,                      &
                       R_Square, TP_Sin_Val,               &
                       R_Int_Weights, TP_Int_Weights    )

!PRINT*,"B"
CALL Calc_Val_On_Elem_TypeA( iE, Cur_Val_Psi, iU_CF, Level )
CALL Calc_Val_On_Elem_TypeA( iE, Cur_Val_AlphaPsi, iU_LF, Level )

!PRINT*,"C"
CALL Calc_Val_And_Drv_On_Elem_TypeB(iE, DROT,      &
                                    CUR_Val_X(:,:,1),       &
                                    CUR_DRV_X(:,:,:,1),     &
                                    iU_X1, iVB_X,           &
                                    Level                   )

CALL Calc_Val_And_Drv_On_Elem_TypeB( iE, DROT,      &
                                    CUR_Val_X(:,:,2),       &
                                    CUR_DRV_X(:,:,:,2),     &
                                    iU_X2, iVB_X,           &
                                    Level                   )

CALL Calc_Val_And_Drv_On_Elem_TypeB( iE, DROT,      &
                                    CUR_Val_X(:,:,3),       &
                                    CUR_DRV_X(:,:,:,3),     &
                                    iU_X3, iVB_X,           &
                                    Level                   )




!PRINT*,"_____Psi______",iE
!PRINT*,Cur_Val_Psi
!PRINT*,"_____AlphaPsi______"
!PRINT*,Cur_Val_AlphaPsi
!PRINT*,"_____X1______"
!PRINT*,Cur_Val_X(:,:,1)
!PRINT*,"_____X2______"
!PRINT*,Cur_Val_X(:,:,2)
!PRINT*,"_____X3______"
!PRINT*,Cur_Val_X(:,:,3)

!PRINT*,"D"
CALL Calc_Ahat( Ahat_Array,                             &
                NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS,  &
                Cur_R_Locs, R_SQUARE,                   &
                TP_RSIN_SQUARE, TP_COTAN_VAL,           &
                CUR_VAL_X, CUR_DRV_X                  )




!PRINT*,"E"
DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS

    ! Calc Flat Metric Terms
    f(tpd,rd,1) = 1.0_idp
    f(tpd,rd,2) = R_Square(rd)
    f(tpd,rd,3) = R_Square(rd) * TP_SIN_SQUARE(tpd)
    
END DO ! tpd
END DO ! rd



!PRINT*,"F"
AA_Array = 0.0_idp
DO i = 1,3
DO j = 1,3

    AA_Array(:,:) = AA_Array(:,:)           &
        + f(:,:,i)*f(:,:,j) * (Ahat_Array(:,:,i,j))**2

END DO ! i
END DO ! j


!PRINT*,"G"
CALL Get_Physical_Source( PhysSrc, iU, iE )


!PRINT*,"H",iU
IF ( iU == iU_CF ) THEN
!    PRINT*,"1"
!    PRINT*,Cur_Val_Psi(:,:)
!    PRINT*,"2"
!    PRINT*,PhysSrc(:,:)
!    PRINT*,"3"
!    PRINT*,AA_Array(:,:)
    SourceTerm(:,:,iU) = -TwoPi * GR_Source_Scalar / Cur_Val_Psi(:,:)   &
                        * PhysSrc(:,:)                                  &
                      -1.0_idp / ( 8.0_idp * Cur_Val_Psi(:,:)**7)       &
                        * AA_Array(:,:)

ELSEIF ( iU == iU_LF ) THEN

!    PRINT*,"1"
!    PRINT*,Cur_Val_Psi(:,:)
!    PRINT*,"2"
!    PRINT*,PhysSrc(:,:)
!    PRINT*,"3"
!    PRINT*,AA_Array(:,:)
!    PRINT*,"4"
!    PRINT*,Cur_Val_AlphaPsi(:,:)
    SourceTerm(:,:,iU) = TwoPi * GR_Source_Scalar * Cur_Val_AlphaPsi(:,:)   &
                        / (Cur_Val_Psi(:,:)**2) * PhysSrc(:,:)              &
                    + 7.0_idp*Cur_Val_AlphaPsi(:,:)                         &
                        / ( 8.0_idp * Cur_Val_Psi(:,:)**8) * AA_Array(:,:)

ENDIF

!PRINT*,"I"


END SUBROUTINE Calc_XCFC_CurVals_TypeA






!+102+##########################################################################!
!                                                                               !
!          XCFC_AMReX_Calc_Source_Vector_TypeA                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_AMReX_Calc_Source_Vector_TypeA( iU )

INTEGER, INTENT(IN)                             ::  iU

#ifdef POSEIDON_AMREX_FLAG

TYPE(amrex_mfiter)                              ::  mfi
TYPE(amrex_box)                                 ::  Box

TYPE(amrex_imultifab)                           ::  Level_Mask

INTEGER                                         ::  re, te, pe
INTEGER, DIMENSION(3)                           ::  iE
INTEGER, DIMENSION(3)                           ::  iEL, iEU
INTEGER                                         ::  nComp
INTEGER                                         ::  lvl

PRINT*,"In XCFC_AMReX_Calc_Source_Vector_TypeA"

FP_Source_Vector_A(:,:,iU) = 0.0_idp

DO lvl = 0,nLevels-1
    IF ( lvl < nLevels-1 ) THEN
        CALL AMReX_MakeFineMask(  Level_Mask,               &
                                  MF_Source(lvl)%ba,        &
                                  MF_Source(lvl)%dm,        &
                                  MF_Source(lvl+1)%ba,      &
                                  iCoarse, iFine            )
    ELSE
        ! Create Level_Mask all equal to 1
        CALL amrex_imultifab_build( Level_Mask,             &
                                    MF_Source(lvl)%ba,      &
                                    MF_Source(lvl)%dm,      &
                                    1,                      &
                                    0                       )
        CALL Level_Mask%SetVal(iCoarse)
    END IF




    CALL amrex_mfiter_build(mfi, MF_Source(lvl), tiling = .false. )
    DO WHILE(mfi%next())

        Source_PTR => MF_Source(lvl)%dataPtr(mfi)
        Mask_PTR   => Level_Mask%dataPtr(mfi)


        Box = mfi%tilebox()

        nComp =  MF_Source(lvl)%ncomp()

        iEL = Box%lo
        iEU = Box%hi

        CALL Initialize_Normed_Legendre_Tables_on_Level( iEU, iEL, lvl )

        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)
            
            IF ( Mask_PTR(RE,TE,PE,1) == iCoarse ) THEN
                CALL Initialize_Ylm_Tables_On_Elem( te, pe, iEL, lvl )
                iE = [re,te,pe]
                CALL XCFC_Calc_Source_Vector_On_Element_TypeA( iU, iE, lvl )
            END IF
        END DO ! pe
        END DO ! te
        END DO ! re

    END DO
    CALL amrex_mfiter_destroy(mfi)
END DO ! lvl

#endif

END SUBROUTINE XCFC_AMReX_Calc_Source_Vector_TypeA



END MODULE XCFC_Source_Vector_TypeA_Module
