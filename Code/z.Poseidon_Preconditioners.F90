   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Preconditioner_Module                                               !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
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
                                ONLY :  idp, pi

USE Poseidon_Parameters, &
                                ONLY :  DEGREE,                     &
                                        L_LIMIT,                    &
                                        NUM_CFA_VARS,               &
                                        NUM_R_ELEMS_PER_BLOCK


USE Poseidon_Variables_Module, &
                                ONLY :  NUM_R_ELEMENTS,             &
                                        rlocs,                      &
                                        INT_R_LOCATIONS,            &
                                        NUM_R_NODES,                &
                                        BLOCK_ELEM_STF_MATVEC,      &
                                        Block_RHS_Vector,           &
                                        Lagrange_Poly_Table,        &
                                        LM_LENGTH,                  &
                                        ULM_LENGTH,                 &
                                        ELEM_PROB_DIM,              &
                                        BLOCK_PROB_DIM


USE Poseidon_Mapping_Functions_Module, &
                                ONLY :  CFA_ALL_Matrix_Map

USE Poseidon_Additional_Functions_Module, &
                        ONLY :  Initialize_LGL_Quadrature_Locations

IMPLICIT NONE

CONTAINS



!+101+###########################################################################!
!                                                                                !
!           Jacobi_Type_A_PC                                                     !
!                                                                                !
!################################################################################!
SUBROUTINE Jacobi_Type_A_PC()

INTEGER                                                     ::  re, d, i, F, lm_loc
REAL(KIND = idp), DIMENSION(0:Block_PROB_DIM-1)             ::  Modifier

INTEGER                                                     ::  Start, Finish, Here

REAL(KIND = idp)                                            ::  Delta_R_Over_Two
REAL(KIND = idp), DIMENSION(0:DEGREE)                       ::  Node_X_Locs,        &
                                                                Node_R_Locs


Node_X_Locs = Initialize_LGL_Quadrature_Locations(DEGREE)


! Create Modifiers
DO re = 0,NUM_R_ELEMS_PER_BLOCK-1

    Delta_R_Over_Two = (rlocs(re+1) - rlocs(re))/2.0_idp
    Node_R_Locs(:) = Delta_R_Over_Two * ( Node_X_Locs(:) + 1.0_idp) + rlocs(re)
    IF ( Node_R_Locs(0) == 0.0_idp ) THEN
        Node_R_Locs(0) = 1.0_idp
    END IF

    DO d = 0,DEGREE
        Start  = CFA_ALL_Matrix_Map(1, 0, re, d)
        Finish = Start + ULM_LENGTH - 1

        Modifier(Start:Finish) = 1.0_idp/Node_R_Locs(d)
!        Modifier(Start:Finish) = 1.0_idp/(Node_R_Locs(d)*Node_R_Locs(d))
!        Modifier(Start:Finish) = 1.0_idp/(Node_R_Locs(d)*Node_R_Locs(d)*Node_R_Locs(d))

    END DO
END DO

! Modify RHS
Block_RHS_Vector(:) = Block_RHS_Vector(:)*Modifier(:)



! Modify Jacobian
DO re = 0,NUM_R_ELEMS_PER_BLOCK-1
    DO d = 0,DEGREE
        DO F = 1,NUM_CFA_VARS
            DO lm_loc = 0,LM_LENGTH-1

                Start  = (d*ULM_LENGTH + (F-1)*LM_LENGTH + lm_loc)*ELEM_PROB_DIM
                Finish = Start + ELEM_PROB_DIM-1
                Here = CFA_ALL_Matrix_Map(F, lm_loc, re, d)
                BLOCK_ELEM_STF_MATVEC(Start:Finish,re) = BLOCK_ELEM_STF_MATVEC(Start:Finish,re) &
                                                       * Modifier(Here)

            END DO ! lm_loc
        END DO ! F
    END DO  ! d
END DO  ! re



END SUBROUTINE Jacobi_Type_A_PC


!+102+###########################################################################!
!                                                                                !
!           Jacobi_Type_B_PC                                                     !
!                                                                                !
!################################################################################!
SUBROUTINE Jacobi_Type_B_PC()

INTEGER                                                     ::  re, d, i, F, lm_loc
REAL(KIND = idp), DIMENSION(0:Block_PROB_DIM-1)             ::  Modifier

INTEGER                                                     ::  Start, Finish, Here

REAL(KIND = idp)                                            ::  Delta_R_Over_Two
REAL(KIND = idp), DIMENSION(0:DEGREE)                       ::  Node_X_Locs,        &
                                                                Node_R_Locs


Node_X_Locs = Initialize_LGL_Quadrature_Locations(DEGREE)


! Create Modifiers
Modifier = 0.0_idp
DO re = 0,NUM_R_ELEMS_PER_BLOCK-1
    DO i = 0,ELEM_PROB_DIM-1

        Start = re*(ELEM_PROB_DIM-ULM_LENGTH) + i
        Here = i*ELEM_PROB_DIM + i
        Modifier(Start) = Modifier(Start)                               &
                        + 1.0_idp/Block_Elem_STF_MATVEC(Here, re)

    END DO
END DO


! Modify RHS
Block_RHS_Vector(:) = Block_RHS_Vector(:)*Modifier(:)



! Modify Jacobian
DO re = 0,NUM_R_ELEMS_PER_BLOCK-1
    DO d = 0,DEGREE
        DO F = 1,NUM_CFA_VARS
            DO lm_loc = 0,LM_LENGTH-1

                Start  = (d*ULM_LENGTH + (F-1)*LM_LENGTH + lm_loc)*ELEM_PROB_DIM
                Finish = Start + ELEM_PROB_DIM-1
                Here = CFA_ALL_Matrix_Map(F, lm_loc, re, d)
                BLOCK_ELEM_STF_MATVEC(Start:Finish,re) = BLOCK_ELEM_STF_MATVEC(Start:Finish,re) &
                                                       * Modifier(Here)

            END DO ! lm_loc
        END DO ! F
    END DO  ! d
END DO  ! re



END SUBROUTINE Jacobi_Type_B_PC

END MODULE Poseidon_Preconditioner_Module
