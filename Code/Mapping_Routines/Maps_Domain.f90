   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Maps_Domain                                                     	     !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!    +101+   Map_To_FEM_Node                                              !##!
!##!    +102+   Map_To_lm                                                    !##!
!##!    +103+   Map_To_tpd                                                   !##!
!##!                                                                         !##!
!###############################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !##########################################################################!


!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE Poseidon_Parameters, &
            ONLY :  Degree
                    
USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,          &
                    Num_T_Quad_Points,          &
                    Num_P_Quad_Points

USE Variables_AMReX_Core, &
            ONLY :  Findloc_Table,      &
                    FEM_Elem_Table,     &
                    Table_Offsets


IMPLICIT NONE


CONTAINS



!+101+###########################################################################!
!                                                                                !
!                  FP_Vector_Map                                                  !
!                                                                                !
!################################################################################!
PURE FUNCTION Map_To_FEM_Node( re, d )

INTEGER                                     :: Map_To_FEM_Node

INTEGER, INTENT(IN)                         :: re, d


Map_To_FEM_Node = re*DEGREE + d + 1


END FUNCTION Map_To_FEM_Node






!+102+###########################################################################!
!                                                                                !
!                  FP_Vector_Map                                                  !
!                                                                                !
!################################################################################!
PURE FUNCTION Map_To_lm( l, m )

INTEGER                                     :: Map_To_lm

INTEGER, INTENT(IN)                         :: l, m


Map_To_lm = l*(l+1) + m + 1


END FUNCTION Map_To_lm




 !+701+################################################################!
!                                                                       !
!                   FEM_Elem_Map                                        !
!                                                                       !
 !#####################################################################!
INTEGER FUNCTION FEM_Elem_Map( AMReX_Elem_Num, AMReX_Level )

INTEGER, INTENT(IN)             ::  AMReX_Elem_Num
INTEGER, INTENT(IN)             ::  AMReX_Level

INTEGER                         ::  i
INTEGER                         ::  Here
INTEGER                         ::  There
INTEGER                         ::  Index


Here  = Table_Offsets(AMReX_Level)
There = Table_Offsets(AMReX_Level+1)-1

!Index = Findloc(Findloc_Table(Here:There), AMReX_Elem_Num, DIM=1)


DO i = Here,There
    IF ( FindLoc_Table(i) == AMReX_Elem_Num ) THEN
        Index = i-Here+1
        EXIT
    ELSE
        Index = -1
    ENDIF
END DO



! Since the arrays start their indexing at 0,
! and the Fortran standard is to start at 1,
! Index-1 must be used as the array index.
FEM_Elem_Map = FEM_Elem_Table(Index+Table_Offsets(AMReX_Level)-1)
                                            

END FUNCTION FEM_Elem_Map



END MODULE Maps_Domain
