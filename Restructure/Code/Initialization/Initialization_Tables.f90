   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Initialization_Tables                                                        !##!
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

USE Parameters_Poseidon, &
        ONLY :  Domain_Dim,     &
                Degree

USE Variables_Quadrature, &
        ONLY :  Num_R_Quad_Points

USE Variables_Tables, &
        ONLY :  M_Values

USE Tables_Module,  &
        ONLY :  Initalize_Ylm_Tables,    &
                Initalize_Lagrange_Poly_Tables

IMPLICIT NONE


CONTAINS
 !+101+################################################################################!
!                                                                                       !
!       Initialize_Tables                                                               !
!                                                                                       !
 !#####################################################################################!
SUBROUTINE Initalize_Tables()


IF ( DOMAIN_DIM == 1 ) THEN
     M_VALUES = 0
ELSE IF ( DOMAIN_DIM == 2 ) THEN
     M_VALUES = 0
ELSE IF ( DOMAIN_DIM == 3 ) THEN
     M_VALUES = (/(l,l=0,L_LIMIT,1)/)
END IF


CALL Initialize_Ylm_Tables()


CALL Initialize_Lagrange_Poly_Tables( DEGREE, NUM_R_QUAD_POINTS )


END SUBROUTINE Initalize_Tables





END MODULE Initialization_Tables


