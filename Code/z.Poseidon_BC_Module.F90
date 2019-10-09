   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_BC_Module                                                           !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   CFA_3D_Apply_BCs_Part1                                              !##!
!##!    +102+   CFA_3D_Apply_BCs_Part2                                              !##!
!##!    +103+   CFA_3D_Dirichlet_BCs_Part1                                          !##!
!##!    +104+   CFA_3D_Dirichlet_BCs_Part2                                          !##!
!##!    +105+   CFA_3D_Neumann_BCs                                                  !##!
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
                                ONLY :  idp,                &
                                        pi,                 &
                                        TwoPi,              &
                                        OneThird,           &
                                        TwoThirds,          &
                                        FourThirds,         &
                                        OneThirtySecond

USE Poseidon_Parameters, &
                                ONLY :  DOMAIN_DIM,                 &
                                        DEGREE,                     &
                                        L_LIMIT,                    &
                                        NUM_SHELLS,                 &
                                        NUM_SUBSHELLS,              &
                                        NUM_SUBSHELLS_PER_SHELL,    &
                                        NUM_CFA_VARS,               &
                                        nPROCS_POSEIDON,            &
                                        NUM_R_QUAD_POINTS,          &
                                        NUM_T_QUAD_POINTS,          &
                                        NUM_P_QUAD_POINTS,          &
                                        NUM_QUAD_DOF,               &
                                        NUM_R_ELEMS_PER_BLOCK,      &
                                        NUM_T_ELEMS_PER_BLOCK,      &
                                        NUM_P_ELEMS_PER_BLOCK,      &
                                        NUM_R_ELEMS_PER_SUBSHELL,   &
                                        NUM_BLOCK_THETA_ROWS,       &
                                        NUM_BLOCK_PHI_COLUMNS,      &
                                        NUM_R_ELEMS_PER_SHELL

USE Poseidon_Variables_Module, &
                                ONLY :  NUM_R_ELEMENTS,             &
                                        NUM_T_ELEMENTS,             &
                                        NUM_P_ELEMENTS,             &
                                        rlocs, tlocs, plocs,        &
                                        NUM_R_NODES,                &
                                        INT_R_LOCATIONS,            &
                                        INT_T_LOCATIONS,            &
                                        INT_P_LOCATIONS,            &
                                        INT_R_WEIGHTS,              &
                                        INT_T_WEIGHTS,              &
                                        INT_P_WEIGHTS,              &
                                        VAR_DIM,                    &
                                        Block_Source_E,             &
                                        Block_Source_S,             &
                                        Block_Source_Si,            &
                                        Ylm_Table_Block,            &
                                        Ylm_Values,                 &
                                        Ylm_CC_Values,              &
                                        Ylm_dt_Values,              &
                                        Ylm_dp_Values,              &
                                        Ylm_dtt_Values,             &
                                        Ylm_dtp_Values,             &
                                        Ylm_dpp_Values,             &
                                        Lagrange_Poly_Table,        &
                                        LPT_LPT,                    &
                                        Coefficient_Vector,         &
                                        nPROCS_SHELL,               &
                                        myID_Poseidon,              &
                                        myID_Shell,                 &
                                        myID_SubShell,              &
                                        myID_PETSc,                 &
                                        INNER_CFA_BC_VALUES,        &
                                        OUTER_CFA_BC_VALUES,        &
                                        INNER_CFA_BC_TYPE,          &
                                        OUTER_CFA_BC_TYPE,          &
                                        PROB_DIM,                   &
                                        Block_PROB_DIM,             &
                                        Coefficient_Vector,         &
                                        NUM_OFF_DIAGONALS,          &
                                        LM_LENGTH,                  &
                                        M_VALUES,                   &
                                        ULM_LENGTH,                 &
                                        Local_Length,               &
                                        POSEIDON_COMM_WORLD,        &
                                        POSEIDON_COMM_SHELL,        &
                                        POSEIDON_COMM_PETSC,        &
                                        Block_STF_MAT,              &
                                        ELEM_PROB_DIM,              &
                                        ELEM_PROB_DIM_SQR,          &
                                        SUBSHELL_PROB_DIM,          &
                                        BLOCK_ELEM_STF_MATVEC,      &
                                        myShell,                    &
                                        Block_RHS_Vector,           &
                                        BLOCK_NONZEROS,             &
                                        Matrix_Location,            &
                                        LM_Location

USE Additional_Functions_Module, &
                                ONLY :  Lagrange_Poly,              &
                                        Spherical_Harmonic,         &
                                        Initialize_LGL_Quadrature,  &
                                        Map_To_X_Space,             &
                                        CFA_Matrix_Map,             &
                                        CFA_ALL_Matrix_Map,         &
                                        Calc_3D_Values_At_Location

IMPLICIT NONE


CONTAINS


!+101+###########################################################################!
!                                                                                !
!                  CFA_3D_Apply_BCs                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_3D_Apply_BCs_Part1(  )

CALL CFA_3D_Dirichlet_BCs_Part1()
CALL CFA_3D_Neumann_BCs()

END SUBROUTINE CFA_3D_Apply_BCs_Part1



!+102+###########################################################################!
!                                                                                !
!                  CFA_3D_Apply_BCs                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_3D_Apply_BCs_Part2(  )

CALL CFA_3D_Dirichlet_BCs_Part2()
CALL CFA_3D_Neumann_BCs()

END SUBROUTINE CFA_3D_Apply_BCs_Part2


!+103+###########################################################################!
!                                                                                !
!                  CFA_3D_Dirichlet_BCs                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_3D_Dirichlet_BCs_Part1( )

INTEGER                                             ::  i, l, m, ui, d

INTEGER                                             ::  Value_Location

!*!
!*!     Steps 1
!*!
DO ui = 1,NUM_CFA_VARS


    !*!
    !*!     Innner BCs
    !*!
    IF ( INNER_CFA_BC_TYPE(ui) == 'D' ) THEN

        DO l = 0,L_LIMIT
            DO m = -M_VALUES(l),M_VALUES(l)

                Value_Location =  Matrix_Location( ui, l, m, 0, 0 )
                Coefficient_Vector( Value_Location ) = 0.0_idp

            END DO
        END DO

        Value_Location =  Matrix_Location( ui, 0, 0, 0, 0 )
        Coefficient_Vector( Value_Location ) = 2.0_idp * sqrt(pi) * INNER_CFA_BC_VALUES(ui)
    END IF





    !*!
    !*!     Outer BCs
    !*!
    IF ( OUTER_CFA_BC_TYPE(ui) == 'D' ) THEN

        DO l = 1,L_LIMIT
            DO m = -M_VALUES(l),M_VALUES(l)

                Value_Location =  Matrix_Location( ui, l, m, NUM_R_ELEMENTS-1, DEGREE )
                Coefficient_Vector( Value_Location ) = 0.0_idp

            END DO
        END DO

        Value_Location =  Matrix_Location( ui, 0, 0, NUM_R_ELEMENTS-1, DEGREE  )
        Coefficient_Vector( Value_Location ) = 2.0_idp * sqrt(pi) * OUTER_CFA_BC_VALUES(ui)

    END IF

END DO


END SUBROUTINE CFA_3D_Dirichlet_BCs_Part1






!+104+###########################################################################!
!                                                                                !
!                  CFA_Dirichlet_BCs                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_3D_Dirichlet_BCs_Part2( )



INTEGER                                             ::  i, l, m, ui, d

INTEGER                                             ::  Value_Location



!*!
!*!     Step 3
!*!
DO ui = 1,NUM_CFA_VARS



    !*!
    !*!     Innner BCs
    !*!

    IF ( INNER_CFA_BC_TYPE(ui) == 'D' .AND. myID_PETSc == 0 ) THEN

        DO l = 0,L_LIMIT
            DO m = -M_VALUES(l),M_VALUES(l)

                Value_Location =  MAtrix_Location( ui, l, m, 0, 0 )
                Block_RHS_VECTOR(Value_Location ) = 0.0_idp

            END DO
        END DO
            !*!
            !*!     Modify the Stiffness Matrix !
            !*!


        !Value_Location =  MAtrix_Location( ui, 0, 0, 0, 0 )
        DO i = 0,ELEM_PROB_DIM-1

          BLOCK_ELEM_STF_MATVEC( i*ELEM_PROB_DIM+ui-1, 0)=0.0_idp

        END DO
!        BLOCK_ELEM_STF_MATVEC(Value_Location*ELEM_PROB_DIM:(Value_Location+1)*ELEM_PROB_DIM-1, 0) = 0.0_idp
!        BLOCK_ELEM_STF_MATVEC(Value_Location*ELEM_PROB_DIM+ Value_Location, 0) = 1.0_idp
        BLOCK_ELEM_STF_MATVEC((ui-1)*ELEM_PROB_DIM:ui*ELEM_PROB_DIM-1, 0) = 0.0_idp
        BLOCK_ELEM_STF_MATVEC((ui-1) + (ui-1)*ELEM_PROB_DIM, 0) = 1.0_idp

    END IF






    !*!
    !*!     Outer BCs
    !*!

    IF ( OUTER_CFA_BC_TYPE(ui) == 'D' .AND. myID_PETSC == NUM_SUBSHELLS-1 ) THEN

        DO l = 0,L_LIMIT
            DO m = -M_VALUES(l),M_VALUES(l)

                Value_Location =  Matrix_Location( ui, l, m, NUM_R_ELEMS_PER_BLOCK-1, DEGREE )
                Block_RHS_VECTOR(Value_Location ) = 0.0_idp

             END DO
         END DO



                !*!
                !*!     Modify the Stiffness Matrix !
                !*!
          Value_Location =  Matrix_Location( ui, 0, 0, 0, DEGREE )

          DO i = 0,ELEM_PROB_DIM-1

              ! Clear the Column !
              BLOCK_ELEM_STF_MATVEC( i*ELEM_PROB_DIM + Value_Location, NUM_R_ELEMS_PER_BLOCK-1) = 0.0_idp

          END DO
          ! Clear the Row !
          BLOCK_ELEM_STF_MATVEC( Value_Location*ELEM_PROB_DIM:(Value_Location+1)*ELEM_PROB_DIM-1,       &
                                 NUM_R_ELEMS_PER_BLOCK-1) = 0.0_idp

          BLOCK_ELEM_STF_MATVEC( Value_Location*ELEM_PROB_DIM+Value_Location, NUM_R_ELEMS_PER_BLOCK-1) = 1.0_idp

    END IF

END DO


END SUBROUTINE CFA_3D_Dirichlet_BCs_Part2







!+105+###########################################################################!
!                                                                                !
!                  CFA_Neumann_BCs                                               !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_3D_Neumann_BCs( )


! Nothing needs to be done. Score!


END SUBROUTINE CFA_3D_Neumann_BCs



END MODULE Poseidon_BC_Module
