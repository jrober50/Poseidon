   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_FP_Source_Vector_Module                                             !##!
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
USE Poseidon_Constants_Module, &
            ONLY :  idp, pi, fdp


USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM,                 &
                    DEGREE,                     &
                    L_LIMIT,                    &
                    NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS

USE Poseidon_Variables_Module, &
            ONLY :  NUM_R_ELEMENTS,             &
                    NUM_T_ELEMENTS,             &
                    NUM_P_ELEMENTS,             &
                    NUM_TP_QUAD_POINTS,         &
                    rlocs, tlocs, plocs,        &
                    NUM_R_NODES,                &
                    INT_R_LOCATIONS,            &
                    INT_T_LOCATIONS,            &
                    INT_P_LOCATIONS,            &
                    INT_R_WEIGHTS,              &
                    INT_T_WEIGHTS,              &
                    INT_P_WEIGHTS,              &
                    INT_TP_WEIGHTS,             &
                    Block_Source_E,             &
                    Block_Source_S,             &
                    Block_Source_Si,            &
                    LM_LENGTH

USE Jacobian_Internal_Functions_Module, &
            ONLY :  JCBN_kappa_FUNCTION_3D_ALL,     &
                    JCBN_BIGK_FUNCTION

USE Poseidon_IO_Module, &
                                ONLY :  Clock_In

USE Poseidon_Mapping_Functions_Module, &
                                ONLY :  Map_To_X_Space

USE SubJacobian_Functions_Module_3D, &
                                ONLY :  Calc_RHS_Terms

USE SubJacobian_Functions_Module_1D, &
                                ONLY :  Calc_RHS_Terms_1D





IMPLICIT NONE

REAL(KIND = idp)                                        :: Time_S

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_R_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_T_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_P_LOCS

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: R_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: R_CUBED

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: SIN_VAL
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: SIN_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: TP_SIN_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COS_VAL
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COS_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CSC_VAL
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CSC_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COTAN_VAL

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: RSIN_SQUARE

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)          :: PHI_EXP
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)          :: PHI_TWOEXP

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
!           Calc_Source_Vector                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_FP_Source_Vector()


REAL(KIND = idp),DIMENSION(1:3)                             ::  Timer

INTEGER                                                     ::  re, te, pe,     &
                                                                d, lm,          &
                                                                rd, tpd, td, pd

REAL(KIND = idp)                                                ::  TWOOVER_DELTAR,    &
                                                                    deltar_overtwo,     &
                                                                    deltat_overtwo,     &
                                                                    deltap_overtwo

REAL(KIND = idp), DIMENSION(1:3,1:3)                            ::  JCBN_kappa_Array
REAL(KIND = idp), DIMENSION(1:3)                                ::  JCBN_n_ARRAY

REAL(KIND = idp)                                                ::  JCBN_BIGK_VALUE

Timer = 0.0_idp

DO pe = 0,NUM_P_ELEMENTS-1
    deltap_overtwo = 0.5_idp * (plocs(pe + 1)-plocs(pe))
    CUR_P_LOCS(:) = deltap_overtwo * (INT_P_LOCATIONS+1.0_idp) + plocs(pe)

    PHI_EXP(:) = EXP( CMPLX(0, -CUR_P_LOCS(:), KIND = idp) )
    PHI_TWOEXP(:) = EXP( CMPLX(0, -2.0_idp*CUR_P_LOCS(:), KIND = idp) )

    DO te = 0,NUM_T_ELEMENTS-1

        deltat_overtwo = 0.5_idp*(tlocs(te + 1) - tlocs(te))
        CUR_T_LOCS(:) = deltat_overtwo * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(te)
        

        COTAN_VAL(:) = 1.0_idp/DTAN(CUR_T_LOCS(:))
        CSC_VAL(:) = 1.0_idp/DSIN(CUR_T_LOCS(:))
        SIN_VAL(:) = DSIN(CUR_T_LOCS(:))
        COS_VAL(:) = DCOS(CUR_T_LOCS(:))

        COS_SQUARE(:) = COS_VAL(:)*COS_VAL(:)
        SIN_SQUARE(:) = SIN_VAL(:)*SIN_VAL(:)
        CSC_SQUARE(:) = CSC_VAL(:)*CSC_VAL(:)

        DO td = 1,NUM_T_QUAD_POINTS
            DO pd = 1,NUM_P_QUAD_POINTS
                tpd = (td-1)*NUM_P_QUAD_POINTS + pd
                TP_SIN_SQUARE(tpd) = SIN_SQUARE(td)
            END DO
        END DO

        DO re = 0,NUM_R_ELEMENTS-1

            deltar_overtwo = 0.5_idp *(rlocs(re + 1) - rlocs(re))
            TWOOVER_DELTAR = 1.0_idp/deltar_overtwo
            CUR_R_LOCS(:) = deltar_overtwo * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(re)


            R_SQUARE(:) = CUR_R_LOCS(:)*CUR_R_LOCS(:)
            DO rd = 1,NUM_R_QUAD_POINTS

                RSIN_SQUARE(:,rd) = R_SQUARE(rd)*SIN_SQUARE(:)

            END DO


!            CALL Calc_FP_Current_Values(Global_re , Local_te  , Local_pe,               &
!                                        DELTAR_OVERTWO, DELTAT_OVERTWO, DELTAP_OVERTWO  )


!            CALL Calc_FP_Source_Terms( Local_re, Local_te, Local_pe )

!            CALL CREATE_FP_Source_Vector( Local_re, Local_te, Local_pe, DELTAR_OVERTWO )


        END DO ! re Loop
    END DO ! te Loop
END DO ! pe Loop



END SUBROUTINE Calc_FP_Source_Vector







!+202+###########################################################################!
!                                                                                !
!                  Calc_Current_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_FP_Current_Values( re, te, pe,                                  &
                                    DELTAR_OVERTWO,                             &
                                    DELTAT_OVERTWO,                             &
                                    DELTAP_OVERTWO                              )

INTEGER, INTENT(IN)                                             ::  re, te, pe


REAL(KIND = idp), INTENT(IN)                                    ::  DELTAR_OVERTWO,     &
                                                                    DELTAT_OVERTWO,     &
                                                                    DELTAP_OVERTWO



COMPLEX(KIND = idp), DIMENSION(1:5)                             ::  Tmp_U_Value,        &
                                                                    Tmp_U_R_DRV_Value,  &
                                                                    Tmp_U_T_DRV_Value,  &
                                                                    Tmp_U_P_DRV_Value



INTEGER                                                         ::  tpd, td, pd, rd,    &
                                                                    lm, d



                          !                                                 !
                         !!                                                 !!
                        !!!          Initialize Local Quadratures           !!!
                         !!                                                 !!
                          !                                                 !
R_Int_Weights(:) = DELTAR_OVERTWO * R_SQUARE(:) * INT_R_WEIGHTS(:)

DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS
        TP_Int_Weights( (td-1)*NUM_P_QUAD_POINTS + pd ) = SIN_VAL(td)                           &
                                                        * DELTAT_OVERTWO * INT_T_WEIGHTS(td)    &
                                                        * DELTAP_OVERTWO * INT_P_WEIGHTS(pd)
    END DO
END DO





    
DO lm = 0,LM_LENGTH
    DO d = 0,DEGREE

        DO tpd = 1,NUM_TP_QUAD_POINTS
            DO rd = 1,NUM_R_QUAD_POINTS

!                CUR_PSI_VALUE(tpd,rd) = SUM( FP_C)

            END DO ! rd
        END DO ! tpd
    END DO ! d
END DO ! lm



END SUBROUTINE Calc_FP_Current_Values













!+701+###########################################################################!
!                                                                                !
!           Allocate_Master_Build_Variables                                      !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_FP_Source_Variables()



ALLOCATE( CUR_R_LOCS(1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_T_LOCS(1:NUM_T_QUAD_POINTS) )
ALLOCATE( CUR_P_LOCS(1:NUM_P_QUAD_POINTS) )


ALLOCATE( R_SQUARE(1:NUM_R_QUAD_POINTS) )
ALLOCATE( R_CUBED(1:NUM_R_QUAD_POINTS) )

ALLOCATE( SIN_VAL( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( SIN_SQUARE( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( TP_SIN_SQUARE( 1:NUM_TP_QUAD_POINTS ) )
ALLOCATE( COS_VAL( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( COS_SQUARE( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( CSC_VAL( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( CSC_SQUARE( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( COTAN_VAL( 1:NUM_T_QUAD_POINTS ) )

ALLOCATE( RSIN_SQUARE( 1:NUM_T_QUAD_POINTS, 1:NUM_R_QUAD_POINTS ) )

ALLOCATE( PHI_EXP( 1:NUM_P_QUAD_POINTS ) )
ALLOCATE( PHI_TWOEXP( 1:NUM_P_QUAD_POINTS ) )

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


END SUBROUTINE Allocate_FP_Source_Variables



!+702+###########################################################################!
!                                                                                !
!           Deallocate_Master_Build_Variables                                    !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_FP_Source_Variables()

DEALLOCATE( CUR_R_LOCS )
DEALLOCATE( CUR_T_LOCS )
DEALLOCATE( CUR_P_LOCS )

DEALLOCATE( R_SQUARE )
DEALLOCATE( R_CUBED )

DEALLOCATE( SIN_VAL )
DEALLOCATE( SIN_SQUARE )
DEALLOCATE( TP_SIN_SQUARE )
DEALLOCATE( COS_VAL )
DEALLOCATE( COS_SQUARE )
DEALLOCATE( CSC_VAL )
DEALLOCATE( CSC_SQUARE )
DEALLOCATE( COTAN_VAL )

DEALLOCATE( RSIN_SQUARE )

DEALLOCATE( PHI_EXP )
DEALLOCATE( PHI_TWOEXP )

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


END SUBROUTINE Deallocate_FP_Source_Variables

END MODULE Poseidon_FP_Source_Vector_Module
