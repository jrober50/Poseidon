   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Maps_Quadrature                                                       !##!
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
USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Num_T_Quad_Points,      &
                    Num_P_Quad_Points




IMPLICIT NONE


INTERFACE Quad_Map
    MODULE PROCEDURE Quad_Map_Short
    MODULE PROCEDURE Quad_Map_Long
    MODULE PROCEDURE Quad_Map_Long_Array
END INTERFACE Quad_Map


CONTAINS



!+101+##########################################################################!
!                                                                               !
!                                                     				!
!                                                                               !
!###############################################################################!
PURE INTEGER FUNCTION Quad_Map_Short(rd, td, pd)

INTEGER, INTENT(IN)     :: rd, td, pd

!Quad_Map_Short =  (rd-1) * Num_P_Quad_Points * Num_T_Quad_Points   &
!                + (td-1) * Num_P_Quad_Points                       &
!                + pd

Quad_Map_Short =  (pd-1) * Num_R_Quad_Points * Num_T_Quad_Points   &
                + (td-1) * Num_R_Quad_Points                       &
                + rd


END FUNCTION Quad_Map_Short



!+101+##########################################################################!
!                                                                               !
!                                                                     !
!                                                                               !
!###############################################################################!
PURE INTEGER FUNCTION Quad_Map_Long( rd, td, pd,     &   ! Quadrature point in question
                                     nr, nt, np      )   ! # of Quad points in each
                                                         ! direction

INTEGER, INTENT(IN)     :: rd, td, pd
INTEGER, INTENT(IN)     :: nr, nt, np

!Quad_Map_Long = (rd-1) * np * nt     &
!              + (td-1) * np          &
!              + pd

Quad_Map_Long = (pd-1) * nr * nt     &
              + (td-1) * nr          &
              + rd

END FUNCTION Quad_Map_Long


!+101+##########################################################################!
!                                                                               !
!                                                                     !
!                                                                               !
!###############################################################################!
PURE INTEGER FUNCTION Quad_Map_Long_Array( rd, td, pd,  &   ! Current Quad Point
                                           nq           )   ! # of Quad points in each
                                                         ! direction, in array.

INTEGER, INTENT(IN)     :: rd, td, pd
INTEGER, INTENT(IN)     :: nq(3)

!Quad_Map_Long_Array = (rd-1) * nq(3) * nq(2)   &
!                    + (td-1) * nq(3)           &
!                    + pd

Quad_Map_Long_Array = (pd-1) * nq(1) * nq(2)   &
                    + (td-1) * nq(1)           &
                    + rd

END FUNCTION Quad_Map_Long_Array







!+103+###########################################################################!
!                                                                                !
!                  FP_Vector_Map                                                  !
!                                                                                !
!################################################################################!
PURE FUNCTION Map_To_tpd( td, pd )

INTEGER                                     :: Map_To_tpd

INTEGER, INTENT(IN)                         :: td, pd


Map_To_tpd = (td-1)*Num_P_Quad_Points + pd


END FUNCTION Map_To_tpd





END MODULE Maps_Quadrature
