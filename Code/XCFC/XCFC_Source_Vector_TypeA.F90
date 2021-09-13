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
           ONLY :  Ylm_CC_Values,              &
                   Lagrange_Poly_Table

USE Variables_Derived, &
           ONLY :  LM_LENGTH

USE Variables_FP, &
           ONLY :  FP_Coeff_Vector_A,            &
                   FP_Source_Vector_A

USE Functions_Jacobian, &
           ONLY :  Calc_Ahat

USE Poseidon_IO_Module, &
           ONLY :  Clock_In

USE FP_Functions_Mapping, &
           ONLY :  FP_FEM_Node_Map,            &
                   FP_tpd_Map

USE XCFC_Source_Functions_Module, &
            ONLY :  Calc_Int_Weights,               &
                    Calc_Val_On_Elem_TypeA,         &
                    Calc_Val_And_Drv_On_Elem_TypeB



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
                    


USE MPI




IMPLICIT NONE


CONTAINS




!+102+###########################################################################!
!                                                                                !
!           XCFC_Calc_Source_Vector_TypeA                                        !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Calc_Source_Vector_TypeA( iU )

INTEGER, INTENT(IN)                             ::  iU

INTEGER                                         ::  re, te, pe,     &
                                                    rd, tpd, td, pd

REAL(KIND = idp)                                ::  deltar_overtwo,     &
                                                    deltat_overtwo,     &
                                                    deltap_overtwo


FP_Source_Vector_A(:,:,iU) = 0.0_idp

DO re = 0,NUM_R_ELEMENTS-1

    DELTAR_OVERTWO = 0.5_idp *(rlocs(re + 1) - rlocs(re))
    CUR_R_LOCS(:) = deltar_overtwo * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(re)
    R_SQUARE(:) = CUR_R_LOCS(:)*CUR_R_LOCS(:)


    DO pe = 0,NUM_P_ELEMENTS-1

        deltap_overtwo = 0.5_idp * (plocs(pe + 1)-plocs(pe))

        DO te = 0,NUM_T_ELEMENTS-1

            deltat_overtwo = 0.5_idp*(tlocs(te + 1) - tlocs(te))
            CUR_T_LOCS(:) = deltat_overtwo * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(te)

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



            CALL Calc_XCFC_CurVals_TypeA( re, te, pe, iU,   &
                                          DELTAR_OVERTWO,   &
                                          DELTAT_OVERTWO,   &
                                          DELTAP_OVERTWO    )

            CALL Create_XCFC_Vector_TypeA( re, te, pe, iU )



        END DO ! te Loop
    END DO ! pe Loop
END DO ! re Loop


END SUBROUTINE XCFC_Calc_Source_Vector_TypeA





!+201+##########################################################################!
!                                                                               !
!                  Create_XCFC_Vector_TypeA                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Create_XCFC_Vector_TypeA( re, te, pe, iU )



INTEGER, INTENT(IN)                                         ::  re, te, pe, iU

INTEGER                                                     ::  rd, d, lm_loc
INTEGER                                                     ::  Current_i_Location

COMPLEX(KIND = idp)                                         ::  RHS_TMP

DO lm_loc = 1,LM_LENGTH
DO d = 0,DEGREE


    RHS_TMP = 0.0_idp
    DO rd = 1,NUM_R_QUAD_POINTS

        RHS_TMP =  RHS_TMP                                         &
                 + SUM( SourceTerm( :, rd, iU )                    &
                       * Ylm_CC_Values( :, lm_loc, te, pe)          &
                       * TP_Int_Weights(:)                     )    &
               * Lagrange_Poly_Table(d, rd, 0)                      &
               * R_Int_Weights(rd)

    END DO  ! rd Loop
    

    Current_i_Location = FP_FEM_Node_Map(re,d)
    FP_Source_Vector_A(Current_i_Location,lm_loc,iU)          &
        = FP_Source_Vector_A(Current_i_Location,lm_loc,iU)    &
        + RHS_TMP

END DO  ! d Loop
END DO  ! lm_loc Loop


END SUBROUTINE Create_XCFC_Vector_TypeA










!+202+###########################################################################!
!                                                                                !
!                  Calc_Current_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_XCFC_CurVals_TypeA( re, te, pe, iU,     &
                                    DROT, DTOT, DPOT    )

INTEGER, INTENT(IN)                                             ::  re, te, pe, iU
REAL(KIND = idp), INTENT(IN)                                    ::  DROT, DTOT, DPOT



REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3,3) ::  Ahat_Array
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3 )  ::  f
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points )    ::  AA_Array

INTEGER                                                         ::  tpd, rd
INTEGER                                                         ::  i, j


CALL Calc_Int_Weights( DROT, DTOT,                      &
                       R_Square, TP_Sin_Val,               &
                       R_Int_Weights, TP_Int_Weights    )

CALL Calc_Val_On_Elem_TypeA( RE, TE, PE, Cur_Val_Psi, iU_CF )
 
CALL Calc_Val_And_Drv_On_Elem_TypeB( RE, TE, PE, DROT,      &
                                    CUR_Val_X(:,:,1),       &
                                    CUR_DRV_X(:,:,:,1),     &
                                    iU_X1, iVB_X            )

CALL Calc_Val_And_Drv_On_Elem_TypeB( RE, TE, PE, DROT,      &
                                    CUR_Val_X(:,:,2),       &
                                    CUR_DRV_X(:,:,:,2),     &
                                    iU_X2, iVB_X            )

CALL Calc_Val_And_Drv_On_Elem_TypeB( RE, TE, PE, DROT,      &
                                    CUR_Val_X(:,:,3),       &
                                    CUR_DRV_X(:,:,:,3),     &
                                    iU_X3, iVB_X            )

CALL Calc_Ahat( Ahat_Array,                             &
                NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS,  &
                Cur_R_Locs, R_SQUARE,                   &
                TP_RSIN_SQUARE, TP_COTAN_VAL,           &
                CUR_VAL_X, CUR_DRV_X                  )



DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_T_QUAD_POINTS

    ! Calc Flat Metric Terms
    f(tpd,rd,1) = 1.0_idp
    f(tpd,rd,2) = R_Square(rd)
    f(tpd,rd,3) = R_Square(rd) * TP_SIN_SQUARE(tpd)
    
END DO ! tpd
END DO ! rd

AA_Array = 0.0_idp
DO i = 1,3
DO j = 1,3
    AA_Array(:,:) = AA_Array(:,:)           &
        + f(:,:,i)*f(:,:,j) * (Ahat_Array(:,:,i,j))**2

END DO ! i
END DO ! j


CALL Get_Source_Term( SourceTerm(:,:,iU), iU, re, te, pe, AA_Array)


END SUBROUTINE Calc_XCFC_CurVals_TypeA







!+701+###########################################################################!
!                                                                                !
!           Get_Physical_Source                                      !
!                                                                                !
!################################################################################!
SUBROUTINE Get_Source_Term( Source, iU, RE, TE, PE, AA )

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points), INTENT(OUT)    :: Source
INTEGER,                                                     INTENT(IN)     :: iU
INTEGER,                                                     INTENT(IN)     :: RE, TE, PE
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points), INTENT(IN)     :: AA

INTEGER                                                         :: rd, td, pd, tpd

IF ( iU == 1 ) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = FP_tpd_Map(td,pd)

        Source(tpd,rd) = -2.0_idp * pi * GR_Source_Scalar                   &
                            / Cur_Val_Psi(tpd,rd)                           &
                            * Block_Source_E(rd,td,pd,re,te,pe)             &
                         - 1.0_idp / ( 8.0_idp * Cur_Val_Psi(tpd,rd)**7)    &
                            * AA(tpd,rd)

    END DO ! pd
    END DO ! td
    END DO ! rd


ELSEIF ( iU == 2 ) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = FP_tpd_Map(td,pd)
        Source(tpd,rd) = 2.0_idp * pi * GR_Source_Scalar                        &
                            * Cur_Val_AlphaPsi(tpd,rd)                          &
                            / (Cur_Val_Psi(tpd,rd)**2)                          &
                            * (Block_Source_E(rd,td,pd,re,te,pe)                &
                                + 2.0_idp*Block_Source_S(rd,td,pd,re,te,pe) )   &
                        + (7.0_idp*Cur_Val_AlphaPsi(tpd,rd))                    &
                            / ( 8.0_idp * Cur_Val_Psi(tpd,rd)**8) * AA(tpd,rd)



    END DO ! pd
    END DO ! td
    END DO ! rd


ELSE IF ( iU == 3) THEN

   DO rd = 1,NUM_R_QUAD_POINTS
   DO td = 1,NUM_T_QUAD_POINTS
   DO pd = 1,NUM_P_QUAD_POINTS

       tpd = FP_tpd_Map(td,pd)
       Source(tpd,rd) = Block_Source_Si(rd,td,pd,re,te,pe,1)


   END DO ! pd
   END DO ! td
   END DO ! rd

ELSE IF ( iU == 4) THEN

   DO rd = 1,NUM_R_QUAD_POINTS
   DO td = 1,NUM_T_QUAD_POINTS
   DO pd = 1,NUM_P_QUAD_POINTS

       tpd = FP_tpd_Map(td,pd)
       Source(tpd,rd) = Block_Source_Si(rd,td,pd,re,te,pe,2)


   END DO ! pd
   END DO ! td
   END DO ! rd

ELSE IF ( iU == 5 ) THEN

   DO rd = 1,NUM_R_QUAD_POINTS
   DO td = 1,NUM_T_QUAD_POINTS
   DO pd = 1,NUM_P_QUAD_POINTS

       tpd = FP_tpd_Map(td,pd)
       Source(tpd,rd) = Block_Source_Si(rd,td,pd,re,te,pe,3)


   END DO ! pd
   END DO ! td
   END DO ! rd

END IF

END SUBROUTINE Get_Source_Term








END MODULE XCFC_Source_Vector_TypeA_Module
