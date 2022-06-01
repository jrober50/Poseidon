   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE FP_Source_Terms_Module                                             !##!
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
            ONLY :  TwoPi,          &
                    Pi

USE Poseidon_Units_Module, &
            ONLY :  GR_Source_Scalar

USE Poseidon_Parameters, &
            ONLY :  DEGREE,                     &
                    L_LIMIT,                    &
                    NUM_CFA_EQs,                &
                    NUM_CFA_VARs

USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,                      &
                    iU_LF,                      &
                    iU_S1,                      &
                    iU_S2,                      &
                    iU_S3


USE Variables_Quadrature, &
            ONLY :  NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    NUM_TP_QUAD_POINTS

USE Variables_Source, &
            ONLY :  Block_Source_E,             &
                    Block_Source_S,             &
                    Block_Source_Si

USE Maps_Quadrature, &
            ONLY :  Quad_Map

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
#endif


USE MPI




IMPLICIT NONE





CONTAINS


!+101+##########################################################################!
!                                                                               !
!           Calc_Source_Terms                                                   !
!                                                                               !
!###############################################################################!
SUBROUTINE Calc_Source_Terms( Source_Terms,                      &
                                     re, te, pe,                        &
                                     td, pd, tpd, rd,                   &
                                     PSI_POWER, ALPHAPSI_POWER,         &
                                     BigK_Value, n_Array, Kappa_Array    )

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:NUM_TP_QUAD_POINTS,   &
                                            1:NUM_R_QUAD_POINTS,    &
                                            1:5                     )   ::  Source_Terms

INTEGER, INTENT(IN)                                                     ::  re, te, pe
INTEGER, INTENT(IN)                                                     ::  td, pd, tpd, rd

REAL(KIND = idp), INTENT(IN), DIMENSION(1:11)                           ::  PSI_POWER
REAL(KIND = idp), INTENT(IN), DIMENSION(1:4)                            ::  ALPHAPSI_POWER


REAL(KIND = idp), INTENT(IN)                                            ::  BigK_VALUE
REAL(KIND = idp), INTENT(IN), DIMENSION(1:3)                            ::  n_ARRAY
REAL(KIND = idp), INTENT(IN), DIMENSION(1:3,1:3)                        ::  Kappa_Array


#ifdef POSEIDON_AMREX_FLAG
CALL Calc_Source_Terms_Native( Source_Terms,                       &
                                re, te, pe,                        &
                                td, pd, tpd, rd,                   &
                                PSI_POWER, ALPHAPSI_POWER,         &
                                BigK_Value, n_Array, Kappa_Array   )
#else
CALL Calc_Source_Terms_Native( Source_Terms,                       &
                                re, te, pe,                        &
                                td, pd, tpd, rd,                   &
                                PSI_POWER, ALPHAPSI_POWER,         &
                                BigK_Value, n_Array, Kappa_Array   )
#endif




END SUBROUTINE Calc_Source_Terms




!+101+##########################################################################!
!                                                                               !
!           Calc_Source_Terms_Native                                            !
!                                                                               !
!###############################################################################!
SUBROUTINE Calc_Source_Terms_Native( Source_Terms,                      &
                                     re, te, pe,                        &
                                     td, pd, tpd, rd,                   &
                                     PSI_POWER, ALPHAPSI_POWER,         &
                                     BigK_Value, n_Array, Kappa_Array    )

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:NUM_TP_QUAD_POINTS,   &
                                            1:NUM_R_QUAD_POINTS,    &
                                            1:5                     )   ::  Source_Terms

INTEGER, INTENT(IN)                                                     ::  re, te, pe
INTEGER, INTENT(IN)                                                     ::  td, pd, tpd, rd

REAL(KIND = idp), INTENT(IN), DIMENSION(1:11)                           ::  PSI_POWER
REAL(KIND = idp), INTENT(IN), DIMENSION(1:4)                            ::  ALPHAPSI_POWER


REAL(KIND = idp), INTENT(IN)                                            ::  BigK_VALUE
REAL(KIND = idp), INTENT(IN), DIMENSION(1:3)                            ::  n_ARRAY
REAL(KIND = idp), INTENT(IN), DIMENSION(1:3,1:3)                        ::  Kappa_Array

REAL(KIND = idp)                                                        ::  Beta_Source_Prefix
INTEGER                                                                 ::  Here

Here = Quad_Map(rd, td, pd)


Source_Terms(tpd, rd, iU_CF) = - TwoPi                                      &
                                * GR_Source_Scalar                          &
                                * Block_Source_E(Here, re, te, pe)    &
                                * PSI_POWER(5)                              &
                             - PSI_POWER(7)                                 &
                                / ( 16.0_idp * ALPHAPSI_POWER(2) )          &
                                * BigK_Value
     




Source_Terms(tpd, rd, iU_LF) = TwoPi                                                        &
                                * ALPHAPSI_POWER(1)                                         &
                                * PSI_POWER(4)                                              &
                                * GR_Source_Scalar                                          &
                                    * ( Block_Source_E(Here, re, te, pe)              &
                                        + 2.0_idp                                           &
                                            * Block_Source_S(Here, re, te, pe)  )     &
                             + 7.0_idp                                                      &
                                * PSI_POWER(6)                                              &
                                / ( 16.0_idp * ALPHAPSI_POWER(1) )                          &
                                * BigK_Value





Beta_Source_Prefix = 16.0_idp * pi * ALPHAPSI_POWER(1) * PSI_POWER(3) * GR_Source_Scalar

Source_Terms(tpd, rd, iU_S1) = Beta_Source_Prefix * Block_Source_Si(Here, re, te, pe, 1)      &
                             + Kappa_Array(1,1) * n_Array(1)                                        &
                             + Kappa_Array(2,1) * n_Array(2)                                        &
                             + Kappa_Array(3,1) * n_Array(3)


Source_Terms(tpd, rd, iU_S2) = Beta_Source_Prefix * Block_Source_Si(Here, re, te, pe, 2)      &
                             + Kappa_Array(1,2) * n_Array(1)                                        &
                             + Kappa_Array(2,2) * n_Array(2)                                        &
                             + Kappa_Array(3,2) * n_Array(3)


Source_Terms(tpd, rd, iU_S3) = Beta_Source_Prefix * Block_Source_Si(Here, re, te, pe, 3)      &
                             + Kappa_Array(1,3) * n_Array(1)                                        &
                             + Kappa_Array(2,3) * n_Array(2)                                        &
                             + Kappa_Array(3,3) * n_Array(3)

!PRINT*,"~~~~~~~~~~~"
!PRINT*,Beta_Source_Prefix * Block_Source_Si(Here, re, te, pe, 3),      &
!                            Kappa_Array(1,3) * n_Array(1),                                        &
!                            Kappa_Array(2,3) * n_Array(2),                                        &
!                            Kappa_Array(3,3) * n_Array(3)
!PRINT*,Kappa_Array(1,3), n_Array(1)

END SUBROUTINE Calc_Source_Terms_Native







#ifdef POSEIDON_AMREX_FLAG
!+102+##########################################################################!
!                                                                               !
!           Calc_Source_Terms_AMReX                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE Calc_Source_Terms_AMReX( Source_Terms,                      &
                                     re, te, pe,                        &
                                     td, pd, tpd, rd,                   &
                                     PSI_POWER, ALPHAPSI_POWER,         &
                                     BigK_Value, n_Array, Kappa_Array    )

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:NUM_TP_QUAD_POINTS,   &
                                            1:NUM_R_QUAD_POINTS,    &
                                            1:5                     )   ::  Source_Terms

INTEGER, INTENT(IN)                                                     ::  re, te, pe
INTEGER, INTENT(IN)                                                     ::  td, pd, tpd, rd

REAL(KIND = idp), INTENT(IN), DIMENSION(1:11)                           ::  PSI_POWER
REAL(KIND = idp), INTENT(IN), DIMENSION(1:4)                            ::  ALPHAPSI_POWER


REAL(KIND = idp), INTENT(IN)                                            ::  BigK_VALUE
REAL(KIND = idp), INTENT(IN), DIMENSION(1:3)                            ::  n_ARRAY
REAL(KIND = idp), INTENT(IN), DIMENSION(1:3,1:3)                        ::  Kappa_Array

REAL(KIND = idp)                                                        ::  Beta_Source_Prefix
INTEGER                                                                 ::  Here





Source_Terms(tpd, rd, 1) = - TwoPi                                      &
                            * GR_Source_Scalar                          &
                            * Block_Source_E(Here, re, te, pe)    &
                            * PSI_POWER(5)                              &
                         - PSI_POWER(7)                                 &
                            / ( 16.0_idp * ALPHAPSI_POWER(2) )          &
                            * BigK_Value
     




Source_Terms(tpd, rd, 2) = TwoPi                                                        &
                            * ALPHAPSI_POWER(1)                                         &
                            * PSI_POWER(4)                                              &
                            * GR_Source_Scalar                                          &
                                * ( Block_Source_E(Here, re, te, pe)              &
                                    + 2.0_idp                                           &
                                        * Block_Source_S(Here, re, te, pe)  )     &
                         + 7.0_idp                                                      &
                            * PSI_POWER(6)                                              &
                            / ( 16.0_idp * ALPHAPSI_POWER(1) )                          &
                            * BigK_Value





Beta_Source_Prefix = 16.0_idp * pi * ALPHAPSI_POWER(1) * PSI_POWER(3) * GR_Source_Scalar

Source_Terms(tpd, rd, 3) = Beta_Source_Prefix * Block_Source_Si(Here, re, te, pe, 1)      &
                         + Kappa_Array(1,1) * n_Array(1)                                        &
                         + Kappa_Array(2,1) * n_Array(2)                                        &
                         + Kappa_Array(3,1) * n_Array(3)


Source_Terms(tpd, rd, 4) = Beta_Source_Prefix * Block_Source_Si(Here, re, te, pe, 2)      &
                         + Kappa_Array(1,2) * n_Array(1)                                        &
                         + Kappa_Array(2,2) * n_Array(2)                                        &
                         + Kappa_Array(3,2) * n_Array(3)


Source_Terms(tpd, rd, 5) = Beta_Source_Prefix * Block_Source_Si(Here, re, te, pe, 3)      &
                         + Kappa_Array(1,3) * n_Array(1)                                        &
                         + Kappa_Array(2,3) * n_Array(2)                                        &
                         + Kappa_Array(3,3) * n_Array(3)




END SUBROUTINE Calc_Source_Terms_AMReX
#endif




END MODULE FP_Source_Terms_Module
