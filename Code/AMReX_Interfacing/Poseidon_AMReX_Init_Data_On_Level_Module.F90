   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Init_Data_On_Level_Module                                    !##!
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
USE ISO_C_BINDING

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
  amrex_multifab_build
#endif


USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Variables_Source, &
            ONLY :  Block_Source_E,             &
                    Block_Source_S,             &
                    Block_Source_Si

USE Variables_Quadrature, &
            ONLY :  NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    Num_Quad_DOF

USE Poseidon_AMReX_Multilayer_Utilities_Module, &
            ONLY :  Find_Coarsest_Parent


IMPLICIT NONE


CONTAINS



 !+101+################################################################!
!                                                                       !
!          Poseidon_Init_Data_On_Level                                  !
!                                                                       !
 !#####################################################################!
SUBROUTINE Poseidon_Init_Data_On_Level( Level,          &
                                        BLo, BHi,       &
                                        SLo, SHi,       &
                                        nComps,         &
                                        Src             )

INTEGER,                INTENT(IN)      ::  Level
INTEGER,                INTENT(IN)      ::  BLo(3), BHi(3)
INTEGER,                INTENT(IN)      ::  SLo(3), SHi(3)
INTEGER,                INTENT(IN)      ::  nComps
REAL(idp),              INTENT(INOUT)   ::  Src(SLo(1):SHi(1),  &
                                                SLo(2):SHi(2),  &
                                                SLo(3):SHi(3),  &
                                                nComps          )

INTEGER                                 ::  re, te, pe
INTEGER                                 ::  rd, td, pd
INTEGER                                 ::  Here
INTEGER, DIMENSION(3)                   ::  CP
INTEGER                                 ::  Num_DOF

Num_DOF = nComps/5


DO pe = BLo(3),BHi(3)
DO te = BLo(2),BHi(2)
DO re = BLo(1),BHi(1)

    CP(1) = Find_Coarsest_Parent(re,Level)
    CP(2) = Find_Coarsest_Parent(te,Level)
    CP(3) = Find_Coarsest_Parent(pe,Level)

    DO pd = 1,Num_P_Quad_Points
    DO td = 1,Num_T_Quad_Points
    DO rd = 1,Num_R_Quad_Points

        here = (pd-1)*Num_T_Quad_Points*Num_R_Quad_Points   &
             + (td-1)*Num_R_Quad_Points                     &
             + rd

        Src(re,te,pe,0*Num_DOF+Here) = Block_Source_E( Here,CP(1),CP(2),CP(3))
        Src(re,te,pe,1*Num_DOF+Here) = Block_Source_S( Here,CP(1),CP(2),CP(3))
        Src(re,te,pe,2*Num_DOF+Here) = Block_Source_Si(Here,CP(1),CP(2),CP(3),1)
        Src(re,te,pe,3*Num_DOF+Here) = Block_Source_Si(Here,CP(1),CP(2),CP(3),2)
        Src(re,te,pe,4*Num_DOF+Here) = Block_Source_Si(Here,CP(1),CP(2),CP(3),3)

    END DO ! rd
    END DO ! td
    END DO ! pd
    
END DO ! re
END DO ! te
END DO ! pe


END SUBROUTINE Poseidon_Init_Data_On_Level









END MODULE Poseidon_Init_Data_On_Level_Module
