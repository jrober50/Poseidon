   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Allocation_Poisson                                                    !##!
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
USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Variables_Poisson, &
            ONLY :  Source_Vector,      &
                    Coefficient_Vector, &
                    Source_Terms,       &
                    First_Column_Storage,   &
                    Last_Column_Storage,    &
                    STF_NNZ,            &
                    STF_Mat_Integrals,  &
                    STF_Elem_Val,       &
                    STF_Row_Ind,        &
                    STF_Col_Ptr

USE Poseidon_Parameters, &
            ONLY :  Degree,             &
                    L_Limit

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,     &
                    Num_T_Elements,     &
                    Num_P_Elements


USE Variables_Quadrature, &
            ONLY :  Num_Quad_DOF

USE Variables_Derived, &
            ONLY :  LM_Length,          &
                    Num_R_Nodes

IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!                                                     				!
!                                                                               !
!###############################################################################!
SUBROUTINE Allocate_Poseidon_Poisson_Variables

ALLOCATE( Source_Vector(0:NUM_R_NODES-1, 1:LM_Length) )
ALLOCATE( Coefficient_Vector(0:NUM_R_NODES-1, -L_LIMIT:L_LIMIT, 0:L_LIMIT) )
ALLOCATE( Source_Terms( 1:Num_Quad_DOF,         &
                        0:NUM_R_ELEMENTS-1,     &
                        0:NUM_T_ELEMENTS-1,     &
                        0:NUM_P_ELEMENTS-1)     )


ALLOCATE( First_Column_Storage(0:Degree, 0:L_Limit) )
ALLOCATE(  Last_Column_Storage(0:Degree, 0:L_Limit) )

Source_Vector      = 0.0_idp
Coefficient_Vector = 0.0_idp


STF_NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1

    !!! Allocate CCS Matrix Arrays !!!
ALLOCATE( STF_MAT_Integrals(0:3, 0:DEGREE, 0:DEGREE) )
ALLOCATE( STF_ELEM_VAL(0:STF_NNZ-1,0:L_LIMIT)        )
ALLOCATE( STF_ROW_IND(0:STF_NNZ-1)                   )
ALLOCATE( STF_COL_PTR(0:NUM_R_NODES)                 )





END SUBROUTINE Allocate_Poseidon_Poisson_Variables





!+101+##########################################################################!
!                                                                               !
!                                                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Deallocate_Poseidon_Poisson_Variables

DEALLOCATE( Source_Vector       )

DEALLOCATE( Coefficient_Vector  )

DEALLOCATE( First_Column_Storage )
DEALLOCATE(  Last_Column_Storage )

DEALLOCATE( STF_MAT_Integrals   )
DEALLOCATE( STF_ELEM_VAL        )
DEALLOCATE( STF_ROW_IND         )
DEALLOCATE( STF_COL_PTR         )

END SUBROUTINE Deallocate_Poseidon_Poisson_Variables








END MODULE Allocation_Poisson
