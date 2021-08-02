   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE XCFC_Lapse_Module                                                   !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!


USE Poseidon_Kinds_Module, &
        ONLY :  idp


USE Poseidon_Parameters, &
        ONLY :  DEGREE,                     &
                L_LIMIT,                    &
                Verbose_Flag,               &
                Convergence_Criteria,       &
                Max_Iterations


USE Variables_Derived, &
        ONLY :  LM_Length,                  &
                Num_R_Nodes,                &
                Prob_Dim,                   &
                Var_Dim


USE Variables_FP,  &
        ONLY :  FP_Coeff_Vector,            &
                FP_Source_Vector,           &
                FP_Update_Vector,           &
                FP_Anderson_M,              &
                MCF_Flag,                   &
                Factored_NNZ,               &
                Laplace_Factored_VAL,       &
                Laplace_Factored_ROW,       &
                Laplace_Factored_COL,       &
                CFA_Var_Map



USE FP_Beta_Banded, &
        ONLY :  Factorize_Beta_Banded,          &
                Jacobi_PC_MVL_Banded_Vector

USE FP_Functions_BC,  &
        ONLY :  DIRICHLET_BC_Beta_Banded

USE FP_Functions_Mapping, &
        ONLY :  FP_LM_Map

USE Poseidon_Cholesky_Module,   &
        ONLY :  CCS_Back_Substitution,          &
                CCS_Forward_Substitution,       &
                Cholesky_Factorization

USE FP_Functions_BC,  &
        ONLY :  Dirichlet_BC_CHOL,              &
                Neumann_BC_CCS

USE XCFC_Source_Vector_Module, &
        ONLY :  XCFC_Calc_Lapse_Source

USE IO_Print_Results, &
        ONLY :  Print_Results


IMPLICIT NONE

CONTAINS

!+101+##########################################################################!
!                                                                               !
!                       XCFC_Lapse_Solve                                      !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Lapse_Solve()



INTEGER                                                 ::  i
INTEGER                                                 ::  here, there
INTEGER                                                 ::  ui
INTEGER                                                 ::  lm_loc

INTEGER                                                 ::  M
INTEGER                                                 ::  LWORK
INTEGER                                                 ::  mk
INTEGER                                                 ::  INFO


COMPLEX(idp),DIMENSION(Var_Dim)                         :: BVector
COMPLEX(idp),DIMENSION(Var_Dim)                         :: UVector
COMPLEX(idp),DIMENSION(Var_Dim)                         :: GVectorM
COMPLEX(idp),DIMENSION(Var_Dim)                         :: FVectorM


COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: FVector
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: GVector
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: AMatrix
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: Work
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: Alpha

LOGICAL                                                 :: PR = .FALSE.
INTEGER                                                 ::  Cur_Iteration       = 0
LOGICAL                                                 ::  CONVERGED           = .FALSE.

ui = 2
M = FP_Anderson_M
LWORK = 2*M



ALLOCATE( FVector(1:Var_Dim,1:M) )
ALLOCATE( GVector(1:Var_Dim,1:M) )
ALLOCATE( Work(1:LWORK) )
ALLOCATE( Alpha(1:M) )
ALLOCATE( AMatrix(1:Var_Dim,1:M))



IF ( Verbose_Flag ) THEN
    PRINT*,"Begining Lapse Function system. "
END IF



CALL XCFC_Calc_Lapse_Source()

Cur_Iteration = 0

DO WHILE ( .NOT. CONVERGED  .AND. Cur_Iteration < Max_Iterations)



    Cur_Iteration = Cur_Iteration+1
    mk = MIN(Cur_Iteration, M)
!    PRINT*,"mk",mk,Cur_Iteration


    IF ( Cur_Iteration .NE. 1 ) THEN
        UVector = GVectorM
    ELSE
        DO lm_loc = 1,LM_Length
            here = (lm_loc-1)*Num_R_Nodes + 1
            there = lm_loc*Num_R_Nodes
            UVector(here:there) = FP_Coeff_Vector(:,lm_loc,ui)
        END DO
    END IF



    !
    !   Solve Systems
    !
    Call XCFC_Solve_Lapse_System()


    DO lm_loc = 1,LM_Length
        here = (lm_loc-1)*Num_R_Nodes + 1
        there = lm_loc*Num_R_Nodes
        GVector(here:there,mk) = FP_Coeff_Vector(:,lm_loc,ui)
    END DO
    FVector(:,mk) = GVector(:,mk) - UVector(:)


    IF ( mk == 1 ) THEN
        GVectorM = GVector(:,mk)
    ELSE
        BVector = -FVector(:,mk)
        AMatrix(:,1:mk-1) = FVector(:,1:mk-1) - SPREAD( FVector(:,mk), DIM=2, NCOPIES = mk-1)


        CALL ZGELS('N',Var_Dim,mk-1,1,              &
                    AMatrix(:,1:mk-1), Var_Dim,     &
                    BVector, Var_Dim,                &
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
        here = (lm_loc-1)*Num_R_Nodes + 1
        there = lm_loc*Num_R_Nodes
        FP_Coeff_Vector(:,lm_loc,ui) = GVectorM(here:There)
    END DO



    CALL XCFC_Calc_Lapse_Source()




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




DEALLOCATE( FVector )
DEALLOCATE( GVector )
DEALLOCATE( Work )
DEALLOCATE( Alpha )
DEALLOCATE( AMatrix )


END  SUBROUTINE XCFC_Lapse_Solve



!+301+###########################################################################!
!                                                                                !
!           Call Solve_FP_System                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Solve_Lapse_System()

COMPLEX(KIND = idp), DIMENSION(NUM_R_NODES)                     ::  WORK_VEC
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                  ::  WORK_ELEM_VAL

INTEGER                                                         ::  ui, l, m
INTEGER                                                         ::  lm_loc

!REAL(idp), DIMENSION(4)                                         ::  timer


IF ( Verbose_Flag ) THEN
    PRINT*,"--In XCFC Iteration, In XCFC_Solve_Lapse_System."
END IF


IF ( MCF_Flag == 0 ) THEN
    
    !
    !   This only needs to be done everytime the radial mesh is defined/redefined.
    !   This performs Cholesky factorization on the stiffness matrix and overwrites
    !   the stiffness matrix variables STF_ELEM_VAL, STF_ROW_IND, and STF_COL_PTR to
    !   represent the factorization matrix, L.  This matrix can then be reused to
    !   solve the linear system using forward and backward substitution.
    !

    CALL Cholesky_Factorization()
    MCF_Flag = 1

END IF


ALLOCATE( Work_Elem_Val(0:Factored_NNZ-1))



ui = 2
DO l = 0,L_LIMIT
DO m = -l,l


    lm_loc = FP_LM_Map(l,m)
    WORK_VEC = -FP_Source_Vector(:,lm_loc,ui)
    WORK_ELEM_VAL(:) = Laplace_Factored_VAL(:,l,CFA_Var_MAP(ui))



!                PRINT*,"Before Dirichelet_BC",ui
    CALL DIRICHLET_BC_CHOL( NUM_R_NODES,                &
                            Factored_NNZ,               &
                            l,                          &
                            m,                          &
                            Laplace_Factored_COL(:,l),  &
                            Laplace_Factored_ROW(:,l),  &
                            WORK_VEC,                   &
                            ui                          )



!                PRINT*,"Before Neumann_BC",ui
    CALL NEUMANN_BC_CCS(    NUM_R_NODES,                &
                            Factored_NNZ,               &
                            l,                          &
                            m,                          &
                            WORK_ELEM_VAL,              &
                            Laplace_Factored_COL(:,l),  &
                            Laplace_Factored_ROW(:,l),  &
                            WORK_VEC                    )


    CALL CCS_Forward_Substitution(  NUM_R_NODES,                    &
                                    Factored_NNZ,                   &
                                    WORK_ELEM_VAL,                  &
                                    Laplace_Factored_COL(:,l),      &
                                    Laplace_Factored_ROW(:,l),      &
                                    WORK_VEC                        )


    CALL CCS_Back_Substitution(     NUM_R_NODES,                    &
                                    Factored_NNZ,                   &
                                    WORK_ELEM_VAL,                  &
                                    Laplace_Factored_COL(:,l),      &
                                    Laplace_Factored_ROW(:,l),      &
                                    WORK_VEC                        )


    FP_Update_Vector(:,lm_loc,ui) = WORK_VEC(:)-FP_Coeff_Vector(:,lm_loc,ui)
    FP_Coeff_Vector( :,lm_loc,ui) = WORK_VEC(:)


    

END DO  ! m Loop
END DO  ! l Loop


DEALLOCATE( Work_Elem_Val )

END SUBROUTINE XCFC_Solve_Lapse_System











END MODULE XCFC_Lapse_Module


