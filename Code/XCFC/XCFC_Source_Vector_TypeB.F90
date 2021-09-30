  !################################################################################!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE XCFC_Source_Vector_TypeB_Module                                              !##!
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

USE Units_Module, &
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
                    DM_Source

USE Variables_AMReX_Core,   &
            ONLY :  AMReX_Levels

USE Variables_AMReX_Source, &
            ONLY :  Source_PTR

#endif

USE MPI




IMPLICIT NONE


CONTAINS
!+101+###########################################################################!
!                                                                                !
!           XCFC_Calc_Source_Vector_TypeB                                        !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Calc_Source_Vector_TypeB( iU, iVB, iEU, iEL )

INTEGER, INTENT(IN), DIMENSION(3)       ::  iU
INTEGER, INTENT(IN)                     ::  iVB
INTEGER, INTENT(IN), DIMENSION(3)       ::  iEU
INTEGER, INTENT(IN), DIMENSION(3)       ::  iEL


#ifdef POSEIDON_AMREX_FLAG

    CALL XCFC_AMReX_Calc_Source_Vector_TypeB( iU, iVB )

#else

    CALL XCFC_Native_Calc_Source_Vector_TypeB( iU, iVB, iEU, iEL )

#endif



END SUBROUTINE XCFC_Calc_Source_Vector_TypeB



!+102+###########################################################################!
!                                                                                !
!           XCFC_Calc_Source_Vector_TypeB                                        !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Native_Calc_Source_Vector_TypeB( iU, iVB, iEU, iEL )

INTEGER, INTENT(IN), DIMENSION(3)               ::  iU
INTEGER, INTENT(IN)                             ::  iVB
INTEGER, INTENT(IN), DIMENSION(3)               ::  iEU
INTEGER, INTENT(IN), DIMENSION(3)               ::  iEL

INTEGER                                         ::  re, te, pe,     &
                                                    rd, tpd, td, pd

REAL(KIND = idp)                                ::  DROT,     &
                                                    DTOT,     &
                                                    DPOT

#ifndef POSEIDON_AMREX_FLAG
FP_Source_Vector_B(:,iVB) = 0.0_idp
#endif


DO re = iEL(1),iEU(1)
DO te = iEL(2),iEU(2)
DO pe = iEL(3),iEU(3)

    DROT = 0.5_idp * (rlocs(re+1) - rlocs(re))
    DPOT = 0.5_idp * (plocs(pe+1) - plocs(pe))
    DTOT = 0.5_idp * (tlocs(te+1) - tlocs(te))

    CUR_R_LOCS(:) = DROT * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(re)
    CUR_T_LOCS(:) = DTOT * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(te)

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

!            PRINT*,R_Square
!            PRINT*,TP_Sin_Val
    CALL Calc_Int_Weights( DROT, DTOT,                  &
                           R_Square, TP_Sin_Val,        &
                           R_Int_Weights, TP_Int_Weights )

    CALL Calc_XCFC_CurVals_TypeB( re, te, pe,       &
                                  iU, iVB,          &
                                  DROT,DTOT, DPOT   )

    CALL Create_XCFC_Vector_TypeB( re, te, pe, iU, iVB )



END DO ! te Loop
END DO ! pe Loop
END DO ! re Loop


END SUBROUTINE XCFC_Native_Calc_Source_Vector_TypeB





!+201+##########################################################################!
!                                                                               !
!                  Create_XCFC_Vector_TypeB                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Create_XCFC_Vector_TypeB( re, te, pe, iU, iVB )

INTEGER, INTENT(IN)                                             ::  re, te, pe
INTEGER, INTENT(IN), DIMENSION(3)                               ::  iU
INTEGER, INTENT(IN)                                             ::  iVB

INTEGER                                                     ::  ui, rd, d, lm_loc
INTEGER                                                     ::  Current_i_Location

COMPLEX(KIND = idp)                                         ::  RHS_TMP


DO ui = iU(1),iU(3)
DO lm_loc = 1,LM_LENGTH
DO d = 0,DEGREE


    RHS_TMP = 0.0_idp
    DO rd = 1,NUM_R_QUAD_POINTS

        RHS_TMP =  RHS_TMP                                         &
                 + SUM( SourceTerm( :, rd, ui )                    &
                       * Ylm_CC_Values( :, lm_loc, te, pe)          &
                       * TP_Int_Weights(:)                     )    &
               * Lagrange_Poly_Table(d, rd, 0)                      &
               * R_Int_Weights(rd)
!        PRINT*,re,te,pe,SourceTerm( :, rd, ui ) ,Ylm_CC_Values( :, lm_loc, te, pe),TP_Int_Weights(:)
!        if ( ui == iU(1) ) THEN
!        PRINT*,"lm ",lm_loc
!        PRINT*,"ST",SourceTerm( :, rd, ui )
!        PRINT*,"++++++++++++"
!        END IF

    END DO  ! rd Loop
    

    Current_i_Location = FP_Array_Map_TypeB(ui,iVB,re,d,lm_loc)
    FP_Source_Vector_B(Current_i_Location,iVB)          &
        = FP_Source_Vector_B(Current_i_Location,iVB)    &
        + RHS_TMP

!    PRINT*,re,te,pe,Current_i_Location,FP_Source_Vector_B(Current_i_Location,iVB)

END DO  ! d Loop
END DO  ! lm_loc Loop
END DO  ! ui Loop

END SUBROUTINE Create_XCFC_Vector_TypeB










!+202+###########################################################################!
!                                                                                !
!                  Calc_Current_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_XCFC_CurVals_TypeB( re, te, pe,         &
                                    iU, iVB,            &
                                    DROT, DTOT, DPOT    )

INTEGER, INTENT(IN)                                             ::  re, te, pe
INTEGER, INTENT(IN), DIMENSION(3)                               ::  iU
INTEGER, INTENT(IN)                                             ::  iVB
REAL(KIND = idp), INTENT(IN)                                    ::  DROT, DTOT, DPOT



REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3,3) ::  Ahat_Array
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3 )  ::  n_Array
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3 )  ::  f
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3)   ::  PhysSrc

INTEGER                                                         ::  tpd, rd
INTEGER                                                         ::  ui, i, ierr

IF ( iVB == iVB_X ) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO tpd = 1,NUM_TP_QUAD_POINTS

        ! Calc Flat Metric Terms
        f(tpd,rd,1) = 1.0_idp
        f(tpd,rd,2) = R_Square(rd)
        f(tpd,rd,3) = R_Square(rd) * TP_SIN_SQUARE(tpd)
    END DO ! tpd
    END DO ! rd

    DO ui = iU(1),iU(3)
        CALL Get_Physical_Source( PhysSrc(:,:,ui-5), ui-3, re, te, pe )
        
!        CALL MPI_Barrier(Poseidon_Comm_World, ierr)
!        DO i = 0,nPROCs_Poseidon-1
!            IF( myID_Poseidon == i ) THEN
!                PRINT*,"myID ",i
!                PRINT*,PhysSrc(:,:,ui-5)
!            END IF
!            CALL MPI_Barrier(Poseidon_Comm_World, ierr)
!        END DO
        SourceTerm(:,:,ui) = 8.0_idp * pi * GR_Source_Scalar * f(:,:,ui-5)*PhysSrc(:,:,ui-5)
    END DO


ELSE IF ( iVB == iVB_S ) THEN


    DO ui = iU(1),iU(3)
        CALL Get_Physical_Source( PhysSrc(:,:,ui-2), ui, re, te, pe )
    END DO

    CALL Calc_Val_And_Drv_On_Elem_TypeA( RE, TE, PE, DROT,    &
                                         Cur_Val_Psi,               &
                                         Cur_DRV_Psi,               &
                                         iU_CF                      )

    CALL Calc_Val_And_Drv_On_Elem_TypeA( RE, TE, PE, DROT,    &
                                         Cur_Val_AlphaPsi,          &
                                         Cur_DRV_AlphaPsi,          &
                                         iU_LF                      )
     
    CALL Calc_Val_And_Drv_On_Elem_TypeB( RE, TE, PE, DROT,    &
                                         CUR_Val_X(:,:,1),          &
                                         CUR_DRV_X(:,:,:,1),        &
                                         iU_X1, iVB_X               )

    CALL Calc_Val_And_Drv_On_Elem_TypeB( RE, TE, PE, DROT,    &
                                         CUR_Val_X(:,:,2),          &
                                         CUR_DRV_X(:,:,:,2),        &
                                         iU_X2, iVB_X               )

    CALL Calc_Val_And_Drv_On_Elem_TypeB( RE, TE, PE, DROT,    &
                                         CUR_Val_X(:,:,3),          &
                                         CUR_DRV_X(:,:,:,3),        &
                                         iU_X3, iVB_X               )


    CALL Calc_Ahat( Ahat_Array,                             &
                    NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS,  &
                    Cur_R_Locs, R_SQUARE,                   &
                    TP_RSIN_SQUARE, TP_COTAN_VAL,           &
                    CUR_VAL_X, CUR_DRV_X                  )



    DO i = 1,3
        n_Array(:,:,i) = Cur_DRV_AlphaPsi(:,:,i)/(Cur_VAL_Psi(:,:)**7)     &
                       - 7.0_idp * Cur_DRV_Psi(:,:,i)                      &
                         * Cur_Val_AlphaPsi(:,:)/(Cur_Val_psi(:,:)**8)
    END DO



    DO ui = 1,3

    SourceTerm(:,:,iU(ui)) = 16.0_idp * pi * GR_Source_Scalar               &   ! Physical Source
                        * Cur_Val_AlphaPsi(:,:)/(Cur_Val_Psi(:,:)**7)   &
                        * PhysSrc(:,:,ui)                               &
                      + 2.0_idp                                         &   ! Geometry Source
                        * ( N_Array(:,:,1)*Ahat_Array(:,:,ui,1)       &
                            + N_Array(:,:,2)*Ahat_Array(:,:,ui,2)     &
                            + N_Array(:,:,3)*Ahat_Array(:,:,ui,3)     )

    END DO




ELSE

    WRITE(*,'(A)')" Error in Poseidon Subroutine : Calc_XCFC_CurVals_TypeB"
    WRITE(*,'(A)')" Invalid input. "
    WRITE(*,'(A)')" Variable : iVB "
    WRITE(*,'(A,I2.2)')" Value    : ",iVB

END IF




END SUBROUTINE Calc_XCFC_CurVals_TypeB







!+101+##########################################################################!
!                                                                               !
!          XCFC_AMReX_Calc_Source_Vector_TypeB                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_AMReX_Calc_Source_Vector_TypeB( iU, iVB )

INTEGER, INTENT(IN), DIMENSION(3)               ::  iU
INTEGER, INTENT(IN)                             ::  iVB

#ifdef POSEIDON_AMREX_FLAG

TYPE(amrex_mfiter)                              ::  mfi
TYPE(amrex_box)                                 ::  Box

INTEGER, DIMENSION(3)                           ::  ELo, EHi
INTEGER                                         ::  nComp
INTEGER                                         ::  lvl

INTEGER                                         ::  ierr



FP_Source_Vector_B(:,iVB) = 0.0_idp
DO lvl = 0,AMReX_Levels-1
    CALL amrex_mfiter_build(mfi, MF_Source(lvl), tiling = .false. )
    DO WHILE(mfi%next())
!        Source_PTR => MF_Source(lvl)%dataPtr(mfi)
        
        Box = mfi%tilebox()

        nComp =  MF_Source(lvl)%ncomp()

        ELo = Box%lo
        EHi = Box%hi


        CALL XCFC_Native_Calc_Source_Vector_TypeB( iU, iVB, EHi-1, ELo-1 )
         

    END DO
    CALL amrex_mfiter_destroy(mfi)
END DO ! lvl


#endif


END SUBROUTINE XCFC_AMReX_Calc_Source_Vector_TypeB





END MODULE XCFC_Source_Vector_TypeB_Module
