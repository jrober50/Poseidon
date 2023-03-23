   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Allocation_Tables                                                            !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
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
            ONLY : idp

USE Poseidon_Message_Routines_Module, &
            ONLY :  Init_Message

USE Poseidon_Parameters, &
            ONLY :  Degree,                 &
                    L_Limit,                &
                    Max_Iterations,         &
                    Verbose_Flag

USE Variables_Derived, &
            ONLY :  LM_Length,              &
                    LM_Short_Length

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    Num_T_Elements,         &
                    Num_P_Elements

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Num_T_Quad_Points,      &
                    Num_P_Quad_Points,      &
                    Num_TP_Quad_Points

USE Variables_Tables, &
            ONLY :  Slm_Elem_Values,            &
                    Slm_Elem_dt_Values,         &
                    Slm_Elem_dp_Values,         &
                    Plm_Values,                 &
                    Plm_dt_Values,              &
                    Nlm_Values,                 &
                    Am_Values,                  &
                    Am_dp_Values,               &
                    Lagrange_Poly_Table,        &
                    LPT_LPT,                    &
                    M_Values,                   &
                    Level_dx,                   &
                    Level_Ratios

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Max_Grid_Size,        &
                    AMReX_Num_Levels


USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Tables_Flags,    &
                    iPF_Init_Tables_Alloc


IMPLICIT NONE


CONTAINS




!+101+##########################################################################!
!                                                                               !
!                            Allocate_Mesh                                      !
!                                                                               !
!###############################################################################!
SUBROUTINE Allocate_Tables()

IF ( Verbose_Flag ) CALL Init_Message('Allocating Table Variables.')

#ifdef POSEIDON_AMREX_FLAG

ALLOCATE( Level_dx( 0:AMReX_Num_Levels-1, 3 ) )
ALLOCATE( Level_Ratios(0:AMReX_Num_Levels) )


ALLOCATE( Plm_Values(       1:Num_T_Quad_Points,        &
                            1:LM_Short_Length,          &
                            0:AMReX_Max_Grid_Size(2)-1) )

ALLOCATE( Plm_dt_Values(    1:Num_T_Quad_Points,        &
                            1:LM_Short_Length,          &
                            0:AMReX_Max_Grid_Size(2)-1) )
                       
ALLOCATE( Nlm_Values(       1:LM_Short_Length)          )

ALLOCATE( Am_Values(        1:Num_P_Quad_Points,        &
                            -L_Limit:L_Limit,           &
                            0:AMReX_Max_Grid_Size(3)-1) )

ALLOCATE( Am_dp_Values(     1:Num_P_Quad_Points,        &
                            -L_Limit:L_Limit,           &
                            0:AMReX_Max_Grid_Size(3)-1) )

ALLOCATE( Slm_Elem_Values(      1:LM_Length,                &
                                1:Num_TP_Quad_Points    )   )
                                
ALLOCATE( Slm_Elem_dt_Values(   1:LM_Length,                &
                                1:Num_TP_Quad_Points    )   )
                                
ALLOCATE( Slm_Elem_dp_Values(   1:LM_Length,                &
                                1:Num_TP_Quad_Points    )   )

#else
                            
                            
ALLOCATE( Plm_Values(       1:Num_T_Quad_Points,        &
                            1:LM_Short_Length,          &
                            0:Num_T_Elements-1)         )

ALLOCATE( Plm_dt_Values(    1:Num_T_Quad_Points,        &
                            1:LM_Short_Length,          &
                            0:Num_T_Elements-1)         )
                            
ALLOCATE( Nlm_Values(       1:LM_Short_Length)          )

ALLOCATE( Am_Values(        1:Num_P_Quad_Points,        &
                            -L_Limit:L_Limit,           &
                            0:Num_P_Elements-1)         )

ALLOCATE( Am_dp_Values(     1:Num_P_Quad_Points,        &
                            -L_Limit:L_Limit,           &
                            0:Num_P_Elements-1)         )
                            
ALLOCATE( Slm_Elem_Values(      1:LM_Length,                &
                                1:Num_TP_Quad_Points    )   )
                                
ALLOCATE( Slm_Elem_dt_Values(   1:LM_Length,                &
                                1:Num_TP_Quad_Points    )   )
                                
ALLOCATE( Slm_Elem_dp_Values(   1:LM_Length,                &
                                1:Num_TP_Quad_Points    )   )

#endif

ALLOCATE( M_VALUES(0:L_LIMIT) )
ALLOCATE( Lagrange_Poly_Table(0:DEGREE, 1:NUM_R_QUAD_POINTS, 0:2)   )
ALLOCATE( LPT_LPT( 1:NUM_R_QUAD_POINTS,0:DEGREE,0:DEGREE,0:1,0:2)       )


lPF_Init_Tables_Flags(iPF_Init_Tables_Alloc) = .TRUE.

END SUBROUTINE Allocate_Tables









!+102+##########################################################################!
!                                                                               !
!                           Deallocate_Mesh                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Deallocate_Tables()

DEALLOCATE( Slm_Elem_Values )
DEALLOCATE( Slm_Elem_dt_Values )
DEALLOCATE( Slm_Elem_dp_Values )

DEALLOCATE( Nlm_Values )

DEALLOCATE( Am_Values )
DEALLOCATE( Am_dp_Values )

DEALLOCATE( Plm_Values )
DEALLOCATE( Plm_dt_Values )

#ifdef POSEIDON_AMREX_FLAG
DEALLOCATE( Level_dx )
DEALLOCATE( Level_Ratios )
#endif

DEALLOCATE( M_Values )


DEALLOCATE( Lagrange_Poly_Table )
DEALLOCATE( LPT_LPT )

END SUBROUTINE Deallocate_Tables





END MODULE Allocation_Tables


