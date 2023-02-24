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


USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Parameters, &
            ONLY :  DEGREE,                     &
                    Num_Vars

USE Variables_Mesh, &
            ONLY :  rlocs,                      &
                    Num_R_Elements

USE Variables_Derived, &
            ONLY :  LM_LENGTH,                  &
                    ULM_LENGTH,                 &
                    ELEM_PROB_DIM,              &
                    Elem_Prob_Dim_Sqr,          &
                    BLOCK_PROB_DIM

USE Maps_Fixed_Point, &
            ONLY :  FP_Array_Map

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations

IMPLICIT NONE

CONTAINS



!+101+###########################################################################!
!                                                                                !
!           Jacobi_Type_A_PC    - 1/r                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE Jacobi_Type_A_PC( A_Mat, b_Vec)

REAL(idp), DIMENSION(0:ELEM_PROB_DIM_SQR-1 ,0:Num_R_Elements-1), INTENT(INOUT) :: A_Mat
REAL(idp), DIMENSION(0:Block_PROB_DIM-1), INTENT(INOUT)     ::  b_Vec

INTEGER                                                     ::  re, d, F, lm_loc
REAL(KIND = idp), DIMENSION(0:Block_PROB_DIM-1)             ::  Modifier

INTEGER                                                     ::  Start, Finish, Here

REAL(KIND = idp)                                            ::  Delta_R_Over_Two
REAL(KIND = idp), DIMENSION(0:DEGREE)                       ::  Node_X_Locs,        &
                                                                Node_R_Locs


Node_X_Locs = Initialize_LGL_Quadrature_Locations(DEGREE)


! Create Modifiers
DO re = 0,Num_R_Elements-1

    Delta_R_Over_Two = (rlocs(re+1) - rlocs(re))/2.0_idp
    Node_R_Locs(:) = Delta_R_Over_Two * ( Node_X_Locs(:) + 1.0_idp) + rlocs(re)
    IF ( Node_R_Locs(0) .LE. 0.0_idp ) THEN
        Node_R_Locs(0) = 1.0_idp
    END IF

    DO d = 0,DEGREE
        Start  = FP_Array_Map(re, d, 1, 0)
        Finish = Start + ULM_LENGTH - 1

        Modifier(Start:Finish) = 1.0_idp/Node_R_Locs(d)
!        Modifier(Start:Finish) = 1.0_idp/(Node_R_Locs(d)*Node_R_Locs(d))
!        Modifier(Start:Finish) = 1.0_idp/(Node_R_Locs(d)*Node_R_Locs(d)*Node_R_Locs(d))

    END DO
END DO

! Modify RHS
b_Vec(:) = b_Vec(:)*Modifier(:)



! Modify Jacobian
DO re = 0,Num_R_Elements-1
    DO d = 0,DEGREE
        DO F = 1,Num_Vars
            DO lm_loc = 0,LM_LENGTH-1

                Start  = (d*ULM_LENGTH + (F-1)*LM_LENGTH + lm_loc)*ELEM_PROB_DIM
                Finish = Start + ELEM_PROB_DIM-1
                Here = FP_Array_Map(re, d, F, lm_loc)

!                PRINT*,"Start: ",start," Finish: ",Finish," Here : ",Here
                A_Mat(Start:Finish,re) = A_Mat(Start:Finish,re) * Modifier(Here)

            END DO ! lm_loc
        END DO ! F
    END DO  ! d
END DO  ! re



END SUBROUTINE Jacobi_Type_A_PC






!+102+###########################################################################!
!                                                                                !
!           Jacobi_Type_B_PC  - Diagonal                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE Jacobi_Type_B_PC( A_Mat, b_Vec )

REAL(idp), DIMENSION(0:ELEM_PROB_DIM_SQR-1 ,0:Num_R_Elements-1), INTENT(INOUT) :: A_Mat
REAL(idp), DIMENSION(0:Block_PROB_DIM-1), INTENT(INOUT)     ::  b_Vec

INTEGER                                                     ::  re, d, i, F, lm_loc
REAL(KIND = idp), DIMENSION(0:Block_PROB_DIM-1)             ::  Modifier

INTEGER                                                     ::  Start, Finish, Here

REAL(KIND = idp), DIMENSION(0:DEGREE)                       ::  Node_X_Locs


Node_X_Locs = Initialize_LGL_Quadrature_Locations(DEGREE)


! Create Modifiers
Modifier = 0.0_idp
DO re = 0,Num_R_Elements-1
    DO i = 0,ELEM_PROB_DIM-1

        Start = re*(ELEM_PROB_DIM-ULM_LENGTH) + i
        Here = i*ELEM_PROB_DIM + i
        Modifier(Start) = 1.0_idp/A_Mat(Here, re)

    END DO
END DO
!PRINT*,Modifier
!PRINT*,"++++++++++++++++++++++++++++++"

! Modify RHS
b_Vec(:) = b_Vec(:)*Modifier(:)



! Modify Jacobian
DO re = 0,Num_R_Elements-1
    DO d = 0,DEGREE
        DO F = 1,Num_Vars
            DO lm_loc = 0,LM_LENGTH-1

                Start  = (d*ULM_LENGTH + (F-1)*LM_LENGTH + lm_loc)*ELEM_PROB_DIM
                Finish = Start + ELEM_PROB_DIM-1
                Here = FP_Array_Map(re, d, F, lm_loc)

                A_Mat(Start:Finish,re) = A_Mat(Start:Finish,re) * Modifier(Here)

            END DO ! lm_loc
        END DO ! F
    END DO  ! d
END DO  ! re



END SUBROUTINE Jacobi_Type_B_PC






END MODULE Poseidon_Preconditioner_Module

