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



USE Poseidon_Kinds_Module, &
            ONLY :  idp
                    

USE Poseidon_Numbers_Module, &
            ONLY :  pi

USE Units_Module, &
            ONLY :  C_Square,           &
                    GR_Source_Scalar

USE Variables_Functions, &
            ONLY :  Potential_Solution,     &
                    Matrix_Location


USE Poseidon_Parameters, &
            ONLY :  DEGREE,                     &
                    L_LIMIT,                    &
                    NUM_CFA_VARS

USE Variables_MPI,  &
            ONLY :  Num_SubShells,              &
                    Num_R_Elems_Per_Block,      &
                    myID_PETSc

USE Variables_Mesh, &
            ONLY :  Num_R_Elements

USE Variables_NR, &
            ONLY :  NR_Coeff_Vector,         &
                    BLOCK_ELEM_STF_MATVEC,      &
                    Block_RHS_Vector

USE Variables_BC, &
            ONLY :  INNER_CFA_BC_VALUES,        &
                    OUTER_CFA_BC_VALUES,        &
                    INNER_CFA_BC_TYPE,          &
                    OUTER_CFA_BC_TYPE
            
USE Variables_Derived, &
            ONLY :  Prob_Dim,                   &
                    Elem_Prob_Dim
                    
USE Variables_Tables, &
            ONLY :  M_Values

USE Functions_Math, &
            ONLY :  Lagrange_Poly

USE Functions_Mapping, &
            ONLY :  Map_To_X_Space

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature,       &
                    Initialize_LGL_Quadrature


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

INTEGER                                             ::  l, m, ui

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
                NR_Coeff_Vector( Value_Location ) = 0.0_idp

            END DO
        END DO

        Value_Location =  Matrix_Location( ui, 0, 0, 0, 0 )
        NR_Coeff_Vector( Value_Location ) = 2.0_idp * sqrt(pi) * INNER_CFA_BC_VALUES(ui)
    END IF





    !*!
    !*!     Outer BCs
    !*!
    IF ( OUTER_CFA_BC_TYPE(ui) == 'D' ) THEN

        DO l = 1,L_LIMIT
            DO m = -M_VALUES(l),M_VALUES(l)

                Value_Location =  Matrix_Location( ui, l, m, NUM_R_ELEMENTS-1, DEGREE )
                NR_Coeff_Vector( Value_Location ) = 0.0_idp

            END DO
        END DO

        Value_Location =  Matrix_Location( ui, 0, 0, NUM_R_ELEMENTS-1, DEGREE  )
        NR_Coeff_Vector( Value_Location ) = 2.0_idp * sqrt(pi) * OUTER_CFA_BC_VALUES(ui)

    END IF

END DO


END SUBROUTINE CFA_3D_Dirichlet_BCs_Part1






!+104+###########################################################################!
!                                                                                !
!                  CFA_Dirichlet_BCs                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_3D_Dirichlet_BCs_Part2( )



INTEGER                                             ::  i, l, m, ui

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

                Value_Location =  MAtrix_Location( ui, l, m, 0, 0 )
                DO i = 0,ELEM_PROB_DIM-1
                    ! Clear the Column !
                    BLOCK_ELEM_STF_MATVEC( i*ELEM_PROB_DIM+value_location, 0)=0.0_idp

                END DO

                ! Clear the Row !
                BLOCK_ELEM_STF_MATVEC(Value_Location*ELEM_PROB_DIM:(Value_Location+1)*ELEM_PROB_DIM-1, 0) = 0.0_idp
                BLOCK_ELEM_STF_MATVEC(Value_Location*ELEM_PROB_DIM+Value_Location, 0) = 1.0_idp


            END DO
        END DO
            !*!
            !*!     Modify the Stiffness Matrix !
            !*!




    END IF






    !*!
    !*!     Outer BCs
    !*!

    IF ( OUTER_CFA_BC_TYPE(ui) == 'D' .AND. myID_PETSC == NUM_SUBSHELLS-1 ) THEN

        DO l = 0,L_LIMIT
            DO m = -M_VALUES(l),M_VALUES(l)

                Value_Location =  Matrix_Location( ui, l, m, NUM_R_ELEMS_PER_BLOCK-1, DEGREE )
                Block_RHS_VECTOR(Value_Location ) = 0.0_idp

                    !*!
                    !*!     Modify the Stiffness Matrix !
                    !*!
                Value_Location =  Matrix_Location( ui, l, m, 0, DEGREE )

                DO i = 0,ELEM_PROB_DIM-1

                    ! Clear the Column !
                    BLOCK_ELEM_STF_MATVEC( i*ELEM_PROB_DIM + Value_Location, NUM_R_ELEMS_PER_BLOCK-1) = 0.0_idp

                END DO
                ! Clear the Row !
                BLOCK_ELEM_STF_MATVEC( Value_Location*ELEM_PROB_DIM:(Value_Location+1)*ELEM_PROB_DIM-1,       &
                                     NUM_R_ELEMS_PER_BLOCK-1) = 0.0_idp

                BLOCK_ELEM_STF_MATVEC( Value_Location*ELEM_PROB_DIM+Value_Location, NUM_R_ELEMS_PER_BLOCK-1) = 1.0_idp

             END DO
         END DO

    END IF

END DO ! ui Loop


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
