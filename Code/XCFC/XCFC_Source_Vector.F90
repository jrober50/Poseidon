   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE XCFC_Source_Vector_Module                                             !##!
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

USE Variables_Quadrature, &
            ONLY :  NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    NUM_TP_QUAD_POINTS,         &
                    INT_R_LOCATIONS,            &
                    INT_T_LOCATIONS,            &
                    INT_P_LOCATIONS,            &
                    INT_R_WEIGHTS,              &
                    INT_T_WEIGHTS,              &
                    INT_P_WEIGHTS,              &
                    INT_TP_WEIGHTS

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
            ONLY :  Ylm_Values,                 &
                    Ylm_dt_Values,              &
                    Ylm_dp_Values,              &
                    Ylm_CC_Values,              &
                    Lagrange_Poly_Table

USE Variables_Derived, &
            ONLY :  NUM_R_NODES,                &
                    LM_LENGTH

USE Variables_FP, &
            ONLY :  FP_Coeff_Vector,            &
                    FP_Coeff_Vector_Beta,       &
                    FP_Source_Vector,           &
                    FP_Source_Vector_Beta,      &
                    FP_Source_Vector_X,         &
                    CFA_EQ_Map,                 &
                    CFA_EQ_Flags

USE Functions_Jacobian, &
            ONLY :  JCBN_kappa_FUNCTION_3D_ALL,     &
                    JCBN_BIGK_FUNCTION,         &
                    JCBN_Kappa_Array_3D

USE Poseidon_IO_Module, &
            ONLY :  Clock_In

USE Functions_Mapping, &
            ONLY :  Map_To_X_Space

USE FP_Functions_Mapping, &
            ONLY :  FP_FEM_Node_Map,            &
                    FP_Beta_Array_Map,          &
                    FP_Array_Map,               &
                    FP_LM_Map


USE MPI




IMPLICIT NONE

REAL(KIND = idp)                                        :: Time_S

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_R_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_T_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_P_LOCS

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: R_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: R_CUBED

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: SIN_VAL_B
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: SIN_VAL
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: SIN_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: TP_SIN_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COS_VAL
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COS_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CSC_VAL
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CSC_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COTAN_VAL

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: RSIN_SQUARE


REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: R_Int_Weights
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: TP_Int_Weights

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: CUR_VAL_PSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: CUR_VAL_ALPHAPSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_VAL_BETA

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_DRV_PSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_DRV_ALPHAPSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:)       :: CUR_DRV_BETA

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: Beta_DRV_Trace

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: SOURCE_TERMS

CONTAINS




!+101+###########################################################################!
!                                                                                !
!           Calc_XCFC_X_Source                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_XCFC_X_Source()



INTEGER                                                     ::  re, te, pe,     &
                                                                d, lm,          &
                                                                rd, tpd, td, pd

REAL(KIND = idp)                                                ::  TWOOVER_DELTAR,    &
                                                                    deltar_overtwo,     &
                                                                    deltat_overtwo,     &
                                                                    deltap_overtwo



FP_Source_Vector_X = 0.0_idp

DO re = 0,NUM_R_ELEMENTS-1

    DELTAR_OVERTWO = 0.5_idp *(rlocs(re + 1) - rlocs(re))
    TWOOVER_DELTAR = 1.0_idp/deltar_overtwo
    CUR_R_LOCS(:) = deltar_overtwo * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(re)


    R_SQUARE(:) = CUR_R_LOCS(:)*CUR_R_LOCS(:)



    DO pe = 0,NUM_P_ELEMENTS-1
        deltap_overtwo = 0.5_idp * (plocs(pe + 1)-plocs(pe))


        DO te = 0,NUM_T_ELEMENTS-1

            deltat_overtwo = 0.5_idp*(tlocs(te + 1) - tlocs(te))

            DO td = 1,NUM_T_QUAD_POINTS
            DO pd = 1,NUM_P_QUAD_POINTS
                tpd = (td-1)*NUM_P_QUAD_POINTS + pd
                SIN_VAL_B(tpd) =   SIN_VAL(td)
                TP_SIN_SQUARE(tpd) = SIN_SQUARE(td)
            END DO
            END DO


            CALL Create_XCFC_X_Vector( re, te, pe,       &
                                        deltat_overtwo,     &
                                        DELTAR_OVERTWO,   &
                                        SIN_VAL_B         )



        END DO ! te Loop
    END DO ! pe Loop


END DO ! re Loop




END SUBROUTINE Calc_XCFC_X_Source





!+204+###########################################################################!
!                                                                                !
!                  CREATE_3D_RHS_VECTOR                                          !
!                                                                                !
!################################################################################!
SUBROUTINE Create_XCFC_X_Vector( re, te, pe,            &
                                 DELTAR_OVERTWO,        &
                                 DeltaT_OverTwo,        &
                                 SIN_VAL            )



INTEGER, INTENT(IN)                                                     ::  re, te, pe

REAL(KIND = idp), INTENT(IN)                                            ::  DELTAR_OVERTWO
REAL(KIND = idp), INTENT(IN)                                            ::  DELTAT_OVERTWO
REAL(KIND = idp), DIMENSION(1:NUM_TP_QUAD_POINTS), INTENT(IN)           ::  SIN_VAL

INTEGER                                                                 ::  pd, td, rd, tpd,     &
                                                                            l, m, d,        &
                                                                            lm_loc, u,ui

INTEGER                                                                 ::  Current_i_Location

COMPLEX(KIND = idp), DIMENSION(1:5)                                     ::  RHS_TMP
COMPLEX(KIND = idp)                                                     ::  Test
COMPLEX(KIND = idp)                                                     ::  Common_Basis
REAL(KIND = idp)                                                        ::  Combined_Weights

COMPLEX(KIND = idp)                                                     ::  Inner, Middle


DO ui = 1,3
DO d = 0,DEGREE
DO lm_loc = 1,LM_LENGTH

    RHS_TMP = 0.0_idp
    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = (td-1)*NUM_P_QUAD_POINTS + pd
        RHS_TMP(ui) =  RHS_TMP(ui)                                          &
                        + 8.0_idp*pi                                        &
                        * Block_Source_Si(rd, td, pd, re, te, pe, ui)       &
                        * Ylm_CC_Values( tpd, lm_loc, te, pe)               &
                        * SIN_VAL(td)                                       &
                        * DELTAT_OVERTWO * INT_T_WEIGHTS(td)                &
                        * INT_P_WEIGHTS(pd)                                 &
                        * Lagrange_Poly_Table(d, rd, 0)                     &
                        * DELTAR_OVERTWO * R_SQUARE(rd) * INT_R_WEIGHTS(rd)

    END DO
    END DO
    END DO  ! rd Loop

    Current_i_Location = FP_Beta_Array_Map(re,d,ui,lm_loc)

    FP_Source_Vector_X(Current_i_Location)                &
        = FP_Source_Vector_X(Current_i_Location)          &
        + RHS_TMP(ui)

END DO  ! lm_loc Loop
END DO  ! d Loop
END DO ! ui


END SUBROUTINE Create_XCFC_X_Vector






















!+102+###########################################################################!
!                                                                                !
!           Calc_XCFC_ConFacto_Source                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_XCFC_ConFactor_Source()



INTEGER                                                     ::  re, te, pe,     &
                                                                d, lm,          &
                                                                rd, tpd, td, pd

REAL(KIND = idp)                                                ::  TWOOVER_DELTAR,    &
                                                                    deltar_overtwo,     &
                                                                    deltat_overtwo,     &
                                                                    deltap_overtwo

REAL(KIND = idp), DIMENSION(1:3,1:3)                            ::  Kappa_Array
REAL(KIND = idp), DIMENSION(1:3)                                ::  n_Array

REAL(KIND = idp)                                                ::  BigK_Value



FP_Source_Vector(:,:,1) = 0.0_idp


! E = psi^6 E



!DO re = 0,NUM_R_ELEMENTS-1
!
!    DELTAR_OVERTWO = 0.5_idp *(rlocs(re + 1) - rlocs(re))
!    TWOOVER_DELTAR = 1.0_idp/deltar_overtwo
!    CUR_R_LOCS(:) = deltar_overtwo * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(re)
!
!
!    R_SQUARE(:) = CUR_R_LOCS(:)*CUR_R_LOCS(:)
!    R_CUBED(:) = R_SQUARE(:)*CUR_R_LOCS(:)
!
!
!
!    DO pe = 0,NUM_P_ELEMENTS-1
!        deltap_overtwo = 0.5_idp * (plocs(pe + 1)-plocs(pe))
!        CUR_P_LOCS(:) = deltap_overtwo * (INT_P_LOCATIONS+1.0_idp) + plocs(pe)
!
!
!        DO te = 0,NUM_T_ELEMENTS-1
!
!!            IF ( re == 64 ) THEN
!!                PRINT*,re,pe,te,FINDLOC(Block_Source_E(:, :, 1, re, te, pe),0.0_idp,2)
!!                PRINT*,Block_Source_E(:, :, 1, re, te, pe)
!!            END IF
!
!            deltat_overtwo = 0.5_idp*(tlocs(te + 1) - tlocs(te))
!            CUR_T_LOCS(:) = deltat_overtwo * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(te)
!            
!
!
!
!            COTAN_VAL(:) = 1.0_idp/DTAN(CUR_T_LOCS(:))
!            CSC_VAL(:) = 1.0_idp/DSIN(CUR_T_LOCS(:))
!            SIN_VAL(:) = DSIN(CUR_T_LOCS(:))
!            COS_VAL(:) = DCOS(CUR_T_LOCS(:))
!
!            COS_SQUARE(:) = COS_VAL(:)*COS_VAL(:)
!            SIN_SQUARE(:) = SIN_VAL(:)*SIN_VAL(:)
!            CSC_SQUARE(:) = CSC_VAL(:)*CSC_VAL(:)
!
!            DO td = 1,NUM_T_QUAD_POINTS
!            DO pd = 1,NUM_P_QUAD_POINTS
!                tpd = (td-1)*NUM_P_QUAD_POINTS + pd
!                SIN_VAL_B(tpd) =   SIN_VAL(td)
!                TP_SIN_SQUARE(tpd) = SIN_SQUARE(td)
!            END DO
!            END DO
!
!            DO rd = 1,NUM_R_QUAD_POINTS
!
!                RSIN_SQUARE(:,rd) = R_SQUARE(rd)*SIN_SQUARE(:)
!
!            END DO
!
!
!
!
!            CALL Calc_FP_Current_Values(  re, te, pe,     &
!                                          DELTAR_OVERTWO, &
!                                          DELTAT_OVERTWO, &
!                                          DELTAP_OVERTWO  )
!            CALL Calc_FP_Source_Terms(    re, te, pe )
!            CALL Create_FP_Source_Vector( re, te, pe,       &
!                                          DELTAR_OVERTWO,   &
!                                          SIN_VAL_B         )
!
!
!
!        END DO ! te Loop
!    END DO ! pe Loop
!END DO ! re Loop

END SUBROUTINE Calc_XCFC_ConFactor_Source








!+103+###########################################################################!
!                                                                                !
!           Calc_XCFC_Lapse_Source                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_XCFC_Lapse_Source()


END SUBROUTINE Calc_XCFC_Lapse_Source


!+104+###########################################################################!
!                                                                                !
!           Calc_XCFC_Shift_Source                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_XCFC_Beta_Source()


END SUBROUTINE Calc_XCFC_Beta_Source





!+701+###########################################################################!
!                                                                                !
!           Allocate_Master_Build_Variables                                      !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_XCFC_Source_Variables()



ALLOCATE( CUR_R_LOCS(1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_T_LOCS(1:NUM_T_QUAD_POINTS) )
ALLOCATE( CUR_P_LOCS(1:NUM_P_QUAD_POINTS) )


ALLOCATE( R_SQUARE(1:NUM_R_QUAD_POINTS) )
ALLOCATE( R_CUBED(1:NUM_R_QUAD_POINTS) )

ALLOCATE( SIN_VAL_B( 1:NUM_TP_QUAD_POINTS ) )
ALLOCATE( SIN_VAL( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( SIN_SQUARE( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( TP_SIN_SQUARE( 1:NUM_TP_QUAD_POINTS ) )
ALLOCATE( COS_VAL( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( COS_SQUARE( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( CSC_VAL( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( CSC_SQUARE( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( COTAN_VAL( 1:NUM_T_QUAD_POINTS ) )

ALLOCATE( RSIN_SQUARE( 1:NUM_T_QUAD_POINTS, 1:NUM_R_QUAD_POINTS ) )


ALLOCATE( R_Int_Weights( 1:NUM_R_QUAD_POINTS ) )
ALLOCATE( TP_Int_Weights( 1:NUM_TP_QUAD_POINTS) )

ALLOCATE( CUR_VAL_PSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_VAL_ALPHAPSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_VAL_BETA(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )

ALLOCATE( CUR_DRV_PSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )
ALLOCATE( CUR_DRV_ALPHAPSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )
ALLOCATE( CUR_DRV_BETA(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3, 1:3) )

ALLOCATE( Beta_DRV_Trace(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS) )

ALLOCATE( Source_Terms( 1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:5) )


END SUBROUTINE Allocate_XCFC_Source_Variables



!+702+###########################################################################!
!                                                                                !
!           Deallocate_Master_Build_Variables                                    !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_XCFC_Source_Variables()

DEALLOCATE( CUR_R_LOCS )
DEALLOCATE( CUR_T_LOCS )
DEALLOCATE( CUR_P_LOCS )

DEALLOCATE( R_SQUARE )
DEALLOCATE( R_CUBED )

DEALLOCATE( SIN_VAL_B )
DEALLOCATE( SIN_VAL )
DEALLOCATE( SIN_SQUARE )
DEALLOCATE( TP_SIN_SQUARE )
DEALLOCATE( COS_VAL )
DEALLOCATE( COS_SQUARE )
DEALLOCATE( CSC_VAL )
DEALLOCATE( CSC_SQUARE )
DEALLOCATE( COTAN_VAL )

DEALLOCATE( RSIN_SQUARE )


DEALLOCATE( R_Int_Weights )
DEALLOCATE( TP_Int_Weights )

DEALLOCATE( CUR_VAL_PSI )
DEALLOCATE( CUR_VAL_ALPHAPSI )
DEALLOCATE( CUR_VAL_BETA )

DEALLOCATE( CUR_DRV_PSI )
DEALLOCATE( CUR_DRV_ALPHAPSI )
DEALLOCATE( CUR_DRV_BETA )

DEALLOCATE( Beta_DRV_Trace )

DEALLOCATE( Source_Terms )


END SUBROUTINE Deallocate_XCFC_Source_Variables

END MODULE XCFC_Source_Vector_Module
