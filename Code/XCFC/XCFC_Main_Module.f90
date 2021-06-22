   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE XCFC_Main_Module                                                         !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the top level subroutines needed to inialize, run, and close       !##!
!##!        Poseidon.                                                               !##!
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

USE Poseidon_Parameters, &
            ONLY :  DEGREE,                     &
                    L_LIMIT,                    &
                    CUR_ITERATION,              &
                    MAX_ITERATIONS,             &
                    Convergence_Type,           &
                    Convergence_Type_Names,     &
                    Convergence_Criteria,       &
                    Convergence_Flag,           &
                    NUM_CFA_EQs,                &
                    Verbose_Flag

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,             &
                    Num_T_Elements,             &
                    Num_P_Elements,             &
                    R_Inner,                    &
                    R_Outer

USE Variables_Derived, &
            ONLY :  LM_Length,                  &
                    Num_R_Nodes,                &
                    Var_Dim,                    &
                    Beta_Prob_Dim,              &
                    Prob_Dim

USE Variables_FP,  &
            ONLY :  Matrix_Format,              &
                    Linear_Solver,              &
                    FP_Source_Vector,           &
                    FP_Coeff_Vector,            &
                    FP_Update_Vector,           &
                    FP_Laplace_Vector,          &
                    FP_Residual_Vector,         &
                    FP_Laplace_Vector_Beta,     &
                    FP_Residual_Vector_Beta,    &
                    Laplace_Matrix_Full,        &
                    Laplace_Matrix_VAL,         &
                    Laplace_Matrix_ROW,         &
                    Laplace_Matrix_COL,         &
                    Laplace_Factored_VAL,       &
                    Laplace_Factored_ROW,       &
                    Laplace_Factored_COL,       &
                    Laplace_Matrix_Beta,        &
                    FP_Source_Vector_Beta,      &
                    FP_Coeff_Vector_Beta,       &
                    FP_Coeff_Vector_X,          &
                    CFA_EQ_Flags,               &
                    CFA_Eq_Map,                 &
                    CFA_MAT_Map,                &
                    Laplace_NNZ,                &
                    Num_Matrices,               &
                    MCF_Flag,                   &
                    FP_Anderson_M


USE Variables_Quadrature, &
            ONLY :  NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    NUM_TP_QUAD_POINTS,         &
                    INT_R_LOCATIONS,            &
                    INT_T_LOCATIONS,            &
                    INT_P_LOCATIONS,            &
                    INT_R_WEIGHTS,              &
                    INT_T_WEIGHTS,              &
                    INT_P_WEIGHTS,              &
                    INT_TP_WEIGHTS

USE Variables_Tables, &
            ONLY :  Ylm_Values,                 &
                    Ylm_dt_Values,              &
                    Ylm_dp_Values,              &
                    Ylm_CC_Values,              &
                    Lagrange_Poly_Table



USE IO_Print_Results, &
            ONLY :  Print_Results


USE XCFC_Source_Vector_Module, &
            ONLY :  Allocate_XCFC_Source_Variables,   &
                    Deallocate_XCFC_Source_Variables,   &
                    Calc_XCFC_X_Source,                 &
                    Calc_XCFC_ConFactor_Source,         &
                    Calc_XCFC_Lapse_Source,             &
                    Calc_XCFC_Beta_Source


USE FP_System_Solvers_Module,   &
            ONLY :  Solve_FP_System,                &
                    Solve_FP_System_Beta,           &
                    Solve_XCFC_System_X

USE FP_Functions_Mapping, &
            ONLY :  FP_FEM_Node_Map,            &
                    FP_Beta_Array_Map,          &
                    FP_Array_Map,               &
                    FP_tpd_Map


IMPLICIT NONE

INTEGER                                                 ::  M
INTEGER                                                 ::  LWORK


CONTAINS


!+101+##########################################################################!
!                                                                               !
!                       XCFC_Method                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Method()




CALL Allocate_XCFC_Source_Variables()


IF ( Verbose_Flag ) THEN
    PRINT*,"Begining XCFC Fixed Point Iterative Solve."
END IF


! Define  E*, Si*


! Solve for X
CALL Solve_X_System()


! Calculate A_ij A^ij/8


! Solve for Conformal Factor
CALL Solve_ConFactor_System()


! Define S*Block_Source_Si(rd, td, pd, re, te, pe, ui)

! Solve for Lapse Function
CALL Solve_Lapse_System()

! Solve for Shift Vector
CALL Solve_Shift_System()




END SUBROUTINE XCFC_Method





!+201+##########################################################################!
!                                                                               !
!                       Solve_X_System                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE Solve_X_System()


LOGICAL                                                 ::  CONVERGED
LOGICAL                                                 ::  CONVERGED_Residual
INTEGER                                                 ::  i, k, lm_loc
INTEGER                                                 ::  here, there
INTEGER                                                 ::  ui, re, d, l



INTEGER                                                 ::  mk
INTEGER                                                 ::  INFO

REAL(KIND = idp), DIMENSION(1:4)                        :: timer

COMPLEX(idp), DIMENSION(1:NUM_R_NODES)                  ::  WORK_VEC
COMPLEX(idp), DIMENSION(1:NUM_R_NODES,1:NUM_R_NODES)    ::  WORK_MAT


COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: Resid_Vector
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: FVector
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: GVector
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: BVector
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: UVector
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: GVectorM
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: FVectorM
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: AMatrix
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: Work
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: Alpha

LOGICAL                                                 :: PR = .FALSE.

M = FP_Anderson_M
LWORK = 2*M

ALLOCATE( Resid_Vector(1:Beta_Prob_Dim) )
ALLOCATE( FVector(1:Beta_Prob_Dim,1:M) )
ALLOCATE( GVectorM(1:Beta_Prob_Dim) )
ALLOCATE( FVectorM(1:Beta_Prob_Dim) )
ALLOCATE( BVector(1:Beta_Prob_Dim) )
ALLOCATE( GVector(1:Beta_Prob_Dim,1:M) )
ALLOCATE( Work(1:LWORK) )
ALLOCATE( Alpha(1:M) )
ALLOCATE( AMatrix(1:Beta_Prob_Dim,1:M))
ALLOCATE( UVector(1:Beta_Prob_Dim) )





IF ( Verbose_Flag ) THEN
    PRINT*,"Begining X system. "
END IF



CALL Calc_XCFC_X_Source()

DO WHILE ( .NOT. CONVERGED  .AND. Cur_Iteration < Max_Iterations)



    Cur_Iteration = Cur_Iteration+1

    mk = MIN(Cur_Iteration, M)





    IF ( Cur_Iteration .NE. 1 ) THEN
        UVector = GVectorM
    ELSE
        UVector = FP_Coeff_Vector_X
    END IF



    !
    !   Solve Systems
    !
!    Call Solve_XCFC_System_X()



    GVector(:,mk) = FP_Coeff_Vector_X(:)
    FVector(:,mk) = GVector(:,mk) - UVector(:)



    IF ( mk == 1 ) THEN
        GVectorM = GVector(:,mk)
    ELSE
        BVector = -FVector(:,mk)
        AMatrix(:,1:mk-1) = FVector(:,1:mk-1) - SPREAD( FVector(:,mk), DIM=2, NCOPIES = mk-1)

        CALL ZGELS('N',Prob_Dim,mk-1,1,              &
                    AMatrix(:,1:mk-1), Prob_Dim,     &
                    BVector, Prob_Dim,                &
                    WORK, LWORK, INFO )

        Alpha(1:mk-1) = BVector(1:mk-1)
        Alpha(mk)     = 1.0_idp - SUM( Alpha(1:mk-1) )
        GVectorM = 0.0_idp
        DO i = 1,mk
            GVectorM = GVectorM + Alpha(i)*GVector(:,i)
        END DO
    END IF
    FVectorM = GVectorM - UVector(:)






    IF ( Verbose_Flag ) THEN
        PRINT*,"L_Inf(FVectorM) = ",MAXVAL(ABS(FVectorM))
    END IF

    IF ( ALL( ABS( FVectorM ) <= Convergence_Criteria ) ) THEN
        IF ( Verbose_Flag ) THEN
            PRINT*,"The Method has converged. The absolute update is less than the tolerance set. "
        END IF
        Converged = .TRUE.
        
    END IF




    IF ( mk == M .AND. .NOT. Converged ) THEN
        GVector = CSHIFT( GVector, SHIFT = +1, DIM = 2)
        FVector = CSHIFT( FVector, SHIFT = +1, DIM = 2)
    END IF



    IF ( Cur_Iteration == Max_Iterations ) THEN
        PRINT*,"FP_Accelerated has reached the maximum number of allowed iterations. "
        Converged = .TRUE.
    END IF




    !
    !   Calculate Source Vector with u_k
    !
    FP_Coeff_Vector_X = GVectorM
    CALL Calc_XCFC_X_Source()




    IF ( PR ) THEN
        WRITE(*,'(A,I3.3,A)')'Iteration ',Cur_Iteration,' Results'
        CALL Print_Results()
        PRINT*," "
    END IF





    IF ( Verbose_Flag .EQV. .TRUE. ) THEN
        WRITE(*,'(A,1X,I3.3,/)') "End of Iteration",Cur_Iteration
    END IF
!    STOP


END DO ! Converged Loop





DEALLOCATE( Resid_Vector )
DEALLOCATE( FVector )
DEALLOCATE( GVectorM)
DEALLOCATE( FVectorM )
DEALLOCATE( BVector )
DEALLOCATE( GVector )
DEALLOCATE( Work )
DEALLOCATE( Alpha )
DEALLOCATE( AMatrix )
DEALLOCATE( UVector )




END  SUBROUTINE Solve_X_System








!+201+##########################################################################!
!                                                                               !
!                       Solve_ConFactor_System                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE Solve_ConFactor_System()

LOGICAL                                                 ::  CONVERGED
LOGICAL                                                 ::  CONVERGED_Residual
INTEGER                                                 ::  i, k, lm_loc
INTEGER                                                 ::  here, there
INTEGER                                                 ::  ui, re, d, l



INTEGER                                                 ::  mk
INTEGER                                                 ::  INFO

REAL(KIND = idp), DIMENSION(1:4)                        :: timer

COMPLEX(idp), DIMENSION(1:NUM_R_NODES)                  ::  WORK_VEC
COMPLEX(idp), DIMENSION(1:NUM_R_NODES,1:NUM_R_NODES)    ::  WORK_MAT


COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: Resid_Vector
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: FVector
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: GVector
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: BVector
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: UVector
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: GVectorM
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: FVectorM
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: AMatrix
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: Work
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: Alpha

LOGICAL                                                 :: PR = .FALSE.

M = FP_Anderson_M
LWORK = 2*M

ALLOCATE( Resid_Vector(1:Var_Dim) )
ALLOCATE( FVector(1:Var_Dim,1:M) )
ALLOCATE( GVectorM(1:Var_Dim) )
ALLOCATE( FVectorM(1:Var_Dim) )
ALLOCATE( BVector(1:Var_Dim) )
ALLOCATE( GVector(1:Var_Dim,1:M) )
ALLOCATE( Work(1:LWORK) )
ALLOCATE( Alpha(1:M) )
ALLOCATE( AMatrix(1:Var_Dim,1:M))
ALLOCATE( UVector(1:Var_Dim) )





IF ( Verbose_Flag ) THEN
    PRINT*,"Begining Conformal Factor system. "
END IF



CALL Calc_XCFC_ConFactor_Source()

DO WHILE ( .NOT. CONVERGED  .AND. Cur_Iteration < Max_Iterations)



    Cur_Iteration = Cur_Iteration+1

    mk = MIN(Cur_Iteration, M)





    IF ( Cur_Iteration .NE. 1 ) THEN
        UVector = GVectorM
    ELSE
        DO lm_loc = 1,LM_Length
            here = (ui-1)*Var_Dim + (lm_loc-1)*Num_R_Nodes + 1
            there = (ui-1)*Var_Dim + lm_loc*Num_R_Nodes
            UVector(here:there) = FP_Coeff_Vector(:,lm_loc,1)
        END DO
    END IF



    !
    !   Solve Systems
    !
!    Call Solve_XCFC_System_ConFactor()



    DO lm_loc = 1,LM_Length
        here = (ui-1)*Var_Dim + (lm_loc-1)*Num_R_Nodes + 1
        there = (ui-1)*Var_Dim + lm_loc*Num_R_Nodes
        GVector(here:there,mk) = FP_Coeff_Vector(:,lm_loc,1)
    END DO
    FVector(:,mk) = GVector(:,mk) - UVector(:)



    IF ( mk == 1 ) THEN
        GVectorM = GVector(:,mk)
    ELSE
        BVector = -FVector(:,mk)
        AMatrix(:,1:mk-1) = FVector(:,1:mk-1) - SPREAD( FVector(:,mk), DIM=2, NCOPIES = mk-1)

        CALL ZGELS('N',Prob_Dim,mk-1,1,              &
                    AMatrix(:,1:mk-1), Prob_Dim,     &
                    BVector, Prob_Dim,                &
                    WORK, LWORK, INFO )

        Alpha(1:mk-1) = BVector(1:mk-1)
        Alpha(mk)     = 1.0_idp - SUM( Alpha(1:mk-1) )
        GVectorM = 0.0_idp
        DO i = 1,mk
            GVectorM = GVectorM + Alpha(i)*GVector(:,i)
        END DO
    END IF
    FVectorM = GVectorM - UVector(:)






    IF ( Verbose_Flag ) THEN
        PRINT*,"L_Inf(FVectorM) = ",MAXVAL(ABS(FVectorM))
    END IF

    IF ( ALL( ABS( FVectorM ) <= Convergence_Criteria ) ) THEN
        IF ( Verbose_Flag ) THEN
            PRINT*,"The Method has converged. The absolute update is less than the tolerance set. "
        END IF
        Converged = .TRUE.
        
    END IF




    IF ( mk == M .AND. .NOT. Converged ) THEN
        GVector = CSHIFT( GVector, SHIFT = +1, DIM = 2)
        FVector = CSHIFT( FVector, SHIFT = +1, DIM = 2)
    END IF



    IF ( Cur_Iteration == Max_Iterations ) THEN
        PRINT*,"FP_Accelerated has reached the maximum number of allowed iterations. "
        Converged = .TRUE.
    END IF




    !
    !   Calculate Source Vector with u_k
    !
    DO lm_loc = 1,LM_Length
        here = (ui-1)*Var_Dim + (lm_loc-1)*Num_R_Nodes + 1
        there = (ui-1)*Var_Dim + lm_loc*Num_R_Nodes
        FP_Coeff_Vector(:,lm_loc,1) = GVectorM(here:There)
    END DO

    CALL Calc_XCFC_ConFactor_Source()




    IF ( PR ) THEN
        WRITE(*,'(A,I3.3,A)')'Iteration ',Cur_Iteration,' Results'
        CALL Print_Results()
        PRINT*," "
    END IF





    IF ( Verbose_Flag .EQV. .TRUE. ) THEN
        WRITE(*,'(A,1X,I3.3,/)') "End of Iteration",Cur_Iteration
    END IF
!    STOP


END DO ! Converged Loop





DEALLOCATE( Resid_Vector )
DEALLOCATE( FVector )
DEALLOCATE( GVectorM)
DEALLOCATE( FVectorM )
DEALLOCATE( BVector )
DEALLOCATE( GVector )
DEALLOCATE( Work )
DEALLOCATE( Alpha )
DEALLOCATE( AMatrix )
DEALLOCATE( UVector )

END  SUBROUTINE Solve_ConFactor_System




!+203+##########################################################################!
!                                                                               !
!                       Solve_Lapse_System                                      !
!                                                                               !
!###############################################################################!
SUBROUTINE Solve_Lapse_System()


LOGICAL                                                 ::  CONVERGED
LOGICAL                                                 ::  CONVERGED_Residual
INTEGER                                                 ::  i, k, lm_loc
INTEGER                                                 ::  here, there
INTEGER                                                 ::  ui, re, d, l



INTEGER                                                 ::  mk
INTEGER                                                 ::  INFO

REAL(KIND = idp), DIMENSION(1:4)                        :: timer

COMPLEX(idp), DIMENSION(1:NUM_R_NODES)                  ::  WORK_VEC
COMPLEX(idp), DIMENSION(1:NUM_R_NODES,1:NUM_R_NODES)    ::  WORK_MAT


COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: Resid_Vector
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: FVector
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: GVector
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: BVector
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: UVector
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: GVectorM
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: FVectorM
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: AMatrix
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: Work
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: Alpha

LOGICAL                                                 :: PR = .FALSE.

M = FP_Anderson_M
LWORK = 2*M

ALLOCATE( Resid_Vector(1:Var_Dim) )
ALLOCATE( FVector(1:Var_Dim,1:M) )
ALLOCATE( GVectorM(1:Var_Dim) )
ALLOCATE( FVectorM(1:Var_Dim) )
ALLOCATE( BVector(1:Var_Dim) )
ALLOCATE( GVector(1:Var_Dim,1:M) )
ALLOCATE( Work(1:LWORK) )
ALLOCATE( Alpha(1:M) )
ALLOCATE( AMatrix(1:Var_Dim,1:M))
ALLOCATE( UVector(1:Var_Dim) )





IF ( Verbose_Flag ) THEN
    PRINT*,"Begining Lapse Function system. "
END IF



CALL Calc_XCFC_Lapse_Source()

DO WHILE ( .NOT. CONVERGED  .AND. Cur_Iteration < Max_Iterations)



    Cur_Iteration = Cur_Iteration+1

    mk = MIN(Cur_Iteration, M)



    IF ( Cur_Iteration .NE. 1 ) THEN
        UVector = GVectorM
    ELSE
        DO lm_loc = 1,LM_Length
            here = (ui-1)*Var_Dim + (lm_loc-1)*Num_R_Nodes + 1
            there = (ui-1)*Var_Dim + lm_loc*Num_R_Nodes
            UVector(here:there) = FP_Coeff_Vector(:,lm_loc,2)
        END DO
    END IF



    !
    !   Solve Systems
    !
!    Call Solve_XCFC_System_Lapse()



    DO lm_loc = 1,LM_Length
        here = (ui-1)*Var_Dim + (lm_loc-1)*Num_R_Nodes + 1
        there = (ui-1)*Var_Dim + lm_loc*Num_R_Nodes
        GVector(here:there,mk) = FP_Coeff_Vector(:,lm_loc,2)
    END DO
    FVector(:,mk) = GVector(:,mk) - UVector(:)



    IF ( mk == 1 ) THEN
        GVectorM = GVector(:,mk)
    ELSE
        BVector = -FVector(:,mk)
        AMatrix(:,1:mk-1) = FVector(:,1:mk-1) - SPREAD( FVector(:,mk), DIM=2, NCOPIES = mk-1)

        CALL ZGELS('N',Prob_Dim,mk-1,1,              &
                    AMatrix(:,1:mk-1), Prob_Dim,     &
                    BVector, Prob_Dim,                &
                    WORK, LWORK, INFO )

        Alpha(1:mk-1) = BVector(1:mk-1)
        Alpha(mk)     = 1.0_idp - SUM( Alpha(1:mk-1) )
        GVectorM = 0.0_idp
        DO i = 1,mk
            GVectorM = GVectorM + Alpha(i)*GVector(:,i)
        END DO
    END IF
    FVectorM = GVectorM - UVector(:)






    IF ( Verbose_Flag ) THEN
        PRINT*,"L_Inf(FVectorM) = ",MAXVAL(ABS(FVectorM))
    END IF

    IF ( ALL( ABS( FVectorM ) <= Convergence_Criteria ) ) THEN
        IF ( Verbose_Flag ) THEN
            PRINT*,"The Method has converged. The absolute update is less than the tolerance set. "
        END IF
        Converged = .TRUE.
        
    END IF




    IF ( mk == M .AND. .NOT. Converged ) THEN
        GVector = CSHIFT( GVector, SHIFT = +1, DIM = 2)
        FVector = CSHIFT( FVector, SHIFT = +1, DIM = 2)
    END IF



    IF ( Cur_Iteration == Max_Iterations ) THEN
        PRINT*,"FP_Accelerated has reached the maximum number of allowed iterations. "
        Converged = .TRUE.
    END IF




    !
    !   Calculate Source Vector with u_k
    !
    DO lm_loc = 1,LM_Length
        here = (ui-1)*Var_Dim + (lm_loc-1)*Num_R_Nodes + 1
        there = (ui-1)*Var_Dim + lm_loc*Num_R_Nodes
        FP_Coeff_Vector(:,lm_loc,2) = GVectorM(here:There)
    END DO

    CALL Calc_XCFC_Lapse_Source()





    IF ( Verbose_Flag .EQV. .TRUE. ) THEN
        WRITE(*,'(A,1X,I3.3,/)') "End of Iteration",Cur_Iteration
    END IF
!    STOP


END DO ! Converged Loop





DEALLOCATE( Resid_Vector )
DEALLOCATE( FVector )
DEALLOCATE( GVectorM)
DEALLOCATE( FVectorM )
DEALLOCATE( BVector )
DEALLOCATE( GVector )
DEALLOCATE( Work )
DEALLOCATE( Alpha )
DEALLOCATE( AMatrix )
DEALLOCATE( UVector )

END  SUBROUTINE Solve_Lapse_System





!+204+##########################################################################!
!                                                                               !
!                       Solve_X_System                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE Solve_Shift_System()



LOGICAL                                                 ::  CONVERGED
LOGICAL                                                 ::  CONVERGED_Residual
INTEGER                                                 ::  i, k, lm_loc
INTEGER                                                 ::  here, there
INTEGER                                                 ::  ui, re, d, l



INTEGER                                                 ::  mk
INTEGER                                                 ::  INFO

REAL(KIND = idp), DIMENSION(1:4)                        :: timer

COMPLEX(idp), DIMENSION(1:NUM_R_NODES)                  ::  WORK_VEC
COMPLEX(idp), DIMENSION(1:NUM_R_NODES,1:NUM_R_NODES)    ::  WORK_MAT


COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: Resid_Vector
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: FVector
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: GVector
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: BVector
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: UVector
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: GVectorM
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: FVectorM
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: AMatrix
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: Work
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: Alpha

LOGICAL                                                 :: PR = .FALSE.

M = FP_Anderson_M
LWORK = 2*M

ALLOCATE( Resid_Vector(1:Beta_Prob_Dim) )
ALLOCATE( FVector(1:Beta_Prob_Dim,1:M) )
ALLOCATE( GVectorM(1:Beta_Prob_Dim) )
ALLOCATE( FVectorM(1:Beta_Prob_Dim) )
ALLOCATE( BVector(1:Beta_Prob_Dim) )
ALLOCATE( GVector(1:Beta_Prob_Dim,1:M) )
ALLOCATE( Work(1:LWORK) )
ALLOCATE( Alpha(1:M) )
ALLOCATE( AMatrix(1:Beta_Prob_Dim,1:M))
ALLOCATE( UVector(1:Beta_Prob_Dim) )





IF ( Verbose_Flag ) THEN
    PRINT*,"Begining Shift system. "
END IF



CALL Calc_XCFC_Beta_Source()

DO WHILE ( .NOT. CONVERGED  .AND. Cur_Iteration < Max_Iterations)



    Cur_Iteration = Cur_Iteration+1

    mk = MIN(Cur_Iteration, M)





    IF ( Cur_Iteration .NE. 1 ) THEN
        UVector = GVectorM
    ELSE
        UVector = FP_Coeff_Vector_Beta
    END IF



    !
    !   Solve Systems
    !
!    Call Solve_XCFC_System_Beta()



    GVector(:,mk) = FP_Coeff_Vector_Beta(:)
    FVector(:,mk) = GVector(:,mk) - UVector(:)



    IF ( mk == 1 ) THEN
        GVectorM = GVector(:,mk)
    ELSE
        BVector = -FVector(:,mk)
        AMatrix(:,1:mk-1) = FVector(:,1:mk-1) - SPREAD( FVector(:,mk), DIM=2, NCOPIES = mk-1)

        CALL ZGELS('N',Prob_Dim,mk-1,1,              &
                    AMatrix(:,1:mk-1), Prob_Dim,     &
                    BVector, Prob_Dim,                &
                    WORK, LWORK, INFO )

        Alpha(1:mk-1) = BVector(1:mk-1)
        Alpha(mk)     = 1.0_idp - SUM( Alpha(1:mk-1) )
        GVectorM = 0.0_idp
        DO i = 1,mk
            GVectorM = GVectorM + Alpha(i)*GVector(:,i)
        END DO
    END IF
    FVectorM = GVectorM - UVector(:)






    IF ( Verbose_Flag ) THEN
        PRINT*,"L_Inf(FVectorM) = ",MAXVAL(ABS(FVectorM))
    END IF

    IF ( ALL( ABS( FVectorM ) <= Convergence_Criteria ) ) THEN
        IF ( Verbose_Flag ) THEN
            PRINT*,"The Method has converged. The absolute update is less than the tolerance set. "
        END IF
        Converged = .TRUE.
        
    END IF




    IF ( mk == M .AND. .NOT. Converged ) THEN
        GVector = CSHIFT( GVector, SHIFT = +1, DIM = 2)
        FVector = CSHIFT( FVector, SHIFT = +1, DIM = 2)
    END IF



    IF ( Cur_Iteration == Max_Iterations ) THEN
        PRINT*,"FP_Accelerated has reached the maximum number of allowed iterations. "
        Converged = .TRUE.
    END IF




    !
    !   Calculate Source Vector with u_k
    !
    FP_Coeff_Vector_Beta = GVectorM
    CALL Calc_XCFC_Beta_Source()




    IF ( Verbose_Flag .EQV. .TRUE. ) THEN
        WRITE(*,'(A,1X,I3.3,/)') "End of Shift Iteration",Cur_Iteration
    END IF
!    STOP


END DO ! Converged Loop





DEALLOCATE( Resid_Vector )
DEALLOCATE( FVector )
DEALLOCATE( GVectorM)
DEALLOCATE( FVectorM )
DEALLOCATE( BVector )
DEALLOCATE( GVector )
DEALLOCATE( Work )
DEALLOCATE( Alpha )
DEALLOCATE( AMatrix )
DEALLOCATE( UVector )

END  SUBROUTINE Solve_Shift_System










!+204+##########################################################################!
!                                                                               !
!                       Solve_X_System                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE Define_Stared_Quantitites()

INTEGER                                                 :: re, te, pe
INTEGER                                                 :: rd, td, pd, tpd
INTEGER                                                 :: d, Here, There

COMPLEX(idp)                                            :: Tmp_Psi

! Calculate Psi^6




! Scale Sources
DO pe = 0,NUM_P_ELEMENTS-1
DO te = 0,NUM_T_ELEMENTS-1
DO re = 0,NUM_R_ELEMENTS-1


DO pd = 1,Num_P_Quad_Points
DO td = 1,Num_T_Quad_Points
DO rd = 1,Num_R_Quad_Points

DO d = 0,DEGREE
    Here = FP_FEM_Node_Map(re,d)
    tpd = FP_tpd_Map(td,pd)

    Tmp_Psi = Tmp_Psi                                   &
            + SUM( FP_Coeff_Vector( Here, :, 1 )        &
            * Ylm_Values( :, tpd, te, pe )       )      &
            * Lagrange_Poly_Table( d, rd, 0 )

END DO ! d


!E_Star(rd, td, pd, re, te, pe) = Block_Source_E(rd, td, pd, re, te, pe)*Tmp_Psi**6
!Si_Star(rd, td, pd, re, te, pe,:) = Block_Source_Si(rd, td, pd, re, te, pe,:)*Tmp_Psi**6



END DO ! pd Loop
END DO ! td Loop
END DO ! rd Loop

END DO ! pe Loop
END DO ! te Loop
END DO ! re Loop



END SUBROUTINE Define_Stared_Quantitites















END MODULE XCFC_Main_Module
