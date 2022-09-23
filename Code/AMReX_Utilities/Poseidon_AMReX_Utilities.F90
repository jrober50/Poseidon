  !##################################################################!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################!
!##!                                                                !##!
!##!                                                                !##!
MODULE Poseidon_AMReX_Utilities_Module                              !##!
!##!                                                                !##!
!##!________________________________________________________________!##!
!##!                                                                !##!
!##!    Contains:                                                   !##!
!##!                                                                !##!
!##!                                                                !##!
!######################################################################!
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  !##################################################################!

#ifdef POSEIDON_AMREX_FLAG


USE Poseidon_Kinds_Module, &
               ONLY :  idp

use amrex_base_module

USE amrex_box_module, &
            ONLY:   amrex_box
USE amrex_multifab_module, &
            ONLY:   amrex_multifab,         &
                    amrex_multifab_build,   &
                    amrex_multifab_destroy, &
                    amrex_mfiter,           &
                    amrex_mfiter_build,     &
                    amrex_mfiter_destroy


IMPLICIT NONE


CONTAINS

 !+101+####################################################!
!                                                           !
!       Initialize_Derived                                  !
!                                                           !
 !#########################################################!
SUBROUTINE Build_MF_From_MF( MF_In,         &
                             MF_In_nComps,  &
                             MF_Num_Levels, &
                             MF_Out         )


TYPE(amrex_multifab),                   INTENT(IN)  ::  MF_In(0:MF_Num_Levels-1)
INTEGER,                                INTENT(IN)  ::  MF_In_nComps
INTEGER,                                INTENT(IN)  ::  MF_Num_Levels

TYPE(amrex_multifab),                   INTENT(OUT) ::  MF_Out(0:MF_Num_Levels-1)

INTEGER                                             ::  Level

DO Level = 0,MF_Num_Levels-1

    CALL amrex_multifab_build(  MF_Out(Level),      &
                                MF_In(Level)%BA,    &
                                MF_In(Level)%DM,    &
                                MF_In_nComps, 1     )

END DO


END SUBROUTINE Build_MF_From_MF




 !+101+####################################################!
!                                                           !
!       Initialize_Derived                                  !
!                                                           !
 !#########################################################!
SUBROUTINE Destroy_MF( MF_In,           &
                       MF_Num_Levels    )


TYPE(amrex_multifab),                   INTENT(OUT) ::  MF_In(0:MF_Num_Levels-1)
INTEGER,                                INTENT(IN)  ::  MF_Num_Levels

INTEGER                                             ::  Level

DO Level = 0,MF_Num_Levels-1

    CALL amrex_multifab_destroy( MF_In(Level) )

END DO


END SUBROUTINE Destroy_MF





#endif


END MODULE Poseidon_AMReX_Utilities_Module
