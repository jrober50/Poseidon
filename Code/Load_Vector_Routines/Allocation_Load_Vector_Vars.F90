   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Load_Vector_Allocation_Module                                         !##!
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

USE Load_Vector_Variables_Module

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,          &
                    Num_T_Quad_Points,          &
                    Num_P_Quad_Points,          &
                    Num_TP_Quad_Points

USE Flags_Core_Module, &
            ONLY :  iPF_Core_Flags,         &
                    iPF_Core_Method_Mode,   &
                    iPF_Core_Method_Newtonian


IMPLICIT NONE



CONTAINS

!+701+###########################################################################!
!                                                                                !
!           Allocate_Master_Build_Variables                                      !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_Load_Vector_Construction_Variables()

ALLOCATE( Cur_R_Locs(1:Num_R_Quad_Points) )
ALLOCATE( Cur_T_Locs(1:Num_T_Quad_Points) )
ALLOCATE( Cur_P_Locs(1:Num_P_Quad_Points) )


ALLOCATE( R_Square(1:Num_R_Quad_Points) )

ALLOCATE( TP_Sin_Val( 1:Num_TP_Quad_Points ) )
ALLOCATE( TP_Sin_Square( 1:Num_TP_Quad_Points ) )
ALLOCATE( TP_RSin_Square( 1:Num_TP_Quad_Points, 1:Num_R_Quad_Points ) )


ALLOCATE( R_Int_Weights( 1:Num_R_Quad_Points ) )
ALLOCATE( TP_Int_Weights( 1:Num_TP_Quad_Points) )


IF ( iPF_Core_Flags(iPF_Core_Method_Mode) .NE. iPF_Core_Method_Newtonian ) THEN

    ALLOCATE( TP_Cotan_Val( 1:Num_TP_Quad_Points ) )

    ALLOCATE( Cur_Val_Psi(1:Num_TP_Quad_Points, 1:Num_R_Quad_Points) )
    ALLOCATE( Cur_Val_AlphaPsi(1:Num_TP_Quad_Points, 1:Num_R_Quad_Points) )
    ALLOCATE( Cur_Val_BETA(1:Num_TP_Quad_Points, 1:Num_R_Quad_Points, 1:3) )
    ALLOCATE( Cur_Val_X(1:Num_TP_Quad_Points, 1:Num_R_Quad_Points, 1:3) )

    ALLOCATE( Cur_Drv_Psi(1:Num_TP_Quad_Points, 1:Num_R_Quad_Points, 1:3) )
    ALLOCATE( Cur_Drv_AlphaPsi(1:Num_TP_Quad_Points, 1:Num_R_Quad_Points, 1:3) )
    ALLOCATE( Cur_Drv_BETA(1:Num_TP_Quad_Points, 1:Num_R_Quad_Points, 1:3, 1:3) )
    ALLOCATE( Cur_Drv_X(1:Num_TP_Quad_Points, 1:Num_R_Quad_Points, 1:3, 1:3) )
    
END IF


IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_Newtonian ) THEN
    ALLOCATE( SourceTerm( 1:Num_TP_Quad_Points, 1:Num_R_Quad_Points, 1 ) )
ELSE
    ALLOCATE( SourceTerm( 1:Num_TP_Quad_Points, 1:Num_R_Quad_Points, 1:8 ) )
END IF

END SUBROUTINE Allocate_Load_Vector_Construction_Variables



!+702+###########################################################################!
!                                                                                !
!           Deallocate_Master_Build_Variables                                    !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_Load_Vector_Construction_Variables()

DEALLOCATE( Cur_R_Locs )
DEALLOCATE( Cur_T_Locs )
DEALLOCATE( Cur_P_Locs )

DEALLOCATE( R_Square )

DEALLOCATE( TP_Sin_Val )
DEALLOCATE( TP_Sin_Square )
DEALLOCATE( TP_RSin_Square )

DEALLOCATE( R_Int_Weights )
DEALLOCATE( TP_Int_Weights )

IF ( iPF_Core_Flags(iPF_Core_Method_Mode) .NE. iPF_Core_Method_Newtonian ) THEN


    DEALLOCATE( TP_Cotan_Val )

    DEALLOCATE( Cur_Val_Psi )
    DEALLOCATE( Cur_Val_AlphaPsi )
    DEALLOCATE( Cur_Val_Beta )
    DEALLOCATE( Cur_Val_X )

    DEALLOCATE( Cur_Drv_Psi )
    DEALLOCATE( Cur_Drv_AlphaPsi )
    DEALLOCATE( Cur_Drv_Beta )
    DEALLOCATE( Cur_Drv_X )

END IF

DEALLOCATE( SourceTerm )

END SUBROUTINE Deallocate_Load_Vector_Construction_Variables




END MODULE Load_Vector_Allocation_Module

