   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE XCFC_Fixed_Point_Module                                               !##!
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

USE Poseidon_Message_Routines_Module, &
            ONLY :  Run_Message,        &
                    Warning_Message

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag,       &
                    Max_Iterations,     &
                    Convergence_Criteria

USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,                        &
                    iU_LF,                        &
                    iU_S1,                        &
                    iU_S2,                        &
                    iU_S3,                        &
                    iU_X1,                        &
                    iU_X2,                        &
                    iU_X3,                       &
                    iVB_S,                       &
                    iVB_X

USE Poseidon_IO_Parameters, &
            ONLY :  CFA_Var_Names

USE Variables_Mesh, &
            ONLY :  NUM_R_ELEMENTS,             &
                    NUM_T_ELEMENTS,             &
                    NUM_P_ELEMENTS

USE Variables_Derived, &
            ONLY :  Var_Dim

USE Variables_FP, &
            ONLY :  FP_Anderson_M

USE XCFC_Load_Vector_TypeA_Module, &
            ONLY :  XCFC_Calc_Load_Vector_TypeA

USE Linear_System_Solvers_TypeA_Module, &
            ONLY :  Solve_Linear_System_TypeA

USE Functions_Coeffs_Module, &
            ONLY :  Coeff_To_Vector_TypeA_SV,  &
                    Vector_To_Coeff_TypeA_SV

USE IO_Print_Results, &
            ONLY :  Print_Results
            
#ifdef POSEIDON_MEMORY_FLAG
USE Poseidon_Memory_Routines, &
            ONLY :  Poseidon_Mark_Memory
 

USE Memory_Variables_Module, &
            ONLY :  Memory_Method_Before_CF_LoadVector1, &
                    Memory_Method_Before_CF_LoadVector2, &
                    Memory_Method_Before_CF_LoadVector3, &
                    Memory_Method_Before_CF_LoadVector4, &
                    Memory_Method_After_CF_LoadVector1,  &
                    Memory_Method_Before_CF_Solve1,     &
                    Memory_Method_Before_CF_Solve2,     &
                    Memory_Method_Before_CF_Solve3,     &
                    Memory_Method_After_CF_Solve1,     &
                    Memory_Method_After_CF_Solve2,     &
                    Memory_Method_After_CF_Solve3,     &
                    Memory_Method_After_CF_FixedPoint,  &
                    Memory_Method_After_CF_DeallocWork, &
                    Memory_Method_Before_LF_LoadVector1, &
                    Memory_Method_Before_LF_LoadVector2, &
                    Memory_Method_Before_LF_LoadVector3, &
                    Memory_Method_Before_LF_LoadVector4, &
                    Memory_Method_Before_LF_Solve1,     &
                    Memory_Method_Before_LF_Solve2,     &
                    Memory_Method_Before_LF_Solve3,     &
                    Memory_Method_After_LF_Solve1,     &
                    Memory_Method_After_LF_Solve2,     &
                    Memory_Method_After_LF_Solve3,     &
                    Memory_Method_After_LF_LoadVector1,  &
                    Memory_Method_After_LF_FixedPoint,  &
                    Memory_Method_After_LF_DeallocWork, &
                    Memory_HWM
  
#endif

IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!                                                     				!
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Fixed_Point(iU)

INTEGER, INTENT(IN)                                     ::  iU

INTEGER, DIMENSION(3)                                   ::  iEU
INTEGER, DIMENSION(3)                                   ::  iEL

INTEGER                                                 ::  i

INTEGER                                                 ::  M
INTEGER                                                 ::  LWORK
INTEGER                                                 ::  mk
INTEGER                                                 ::  INFO

!REAL(KIND = idp), DIMENSION(1:4)                        :: timer

REAL(idp),  DIMENSION(Var_Dim)                          :: BVector
REAL(idp),  DIMENSION(Var_Dim)                          :: UVector
REAL(idp),  DIMENSION(Var_Dim)                          :: GVectorM
REAL(idp),  DIMENSION(Var_Dim)                          :: FVectorM


REAL(idp),  DIMENSION(Var_Dim,FP_Anderson_M)            :: FVector
REAL(idp),  DIMENSION(Var_Dim,FP_Anderson_M)            :: GVector
REAL(idp),  DIMENSION(Var_Dim,FP_Anderson_M)            :: AMatrix

REAL(idp),  DIMENSION(FP_Anderson_M)                    :: Alpha

REAL(idp),  DIMENSION(:),   ALLOCATABLE                 :: Work

INTEGER                                                 ::  Cur_Iteration
LOGICAL                                                 ::  CONVERGED

CHARACTER(LEN = 300)                                    ::  Message


IF ( Verbose_Flag ) THEN
    WRITE(Message,'(A,A,A)')'Beginning ',TRIM(CFA_Var_Names(iU)),' fixed point iterations.'
    CALL Run_Message(TRIM(Message))
END IF




INFO = 0
M = FP_Anderson_M
LWORK = 2*M

iEL = [0, 0, 0]
iEU = [Num_R_Elements-1,Num_T_Elements-1,Num_P_Elements-1]
ALLOCATE( Work(1:LWORK) )

#ifdef POSEIDON_MEMORY_FLAG
IF ( iU == iU_CF ) THEN
    CALL Poseidon_Mark_Memory(Memory_Method_Before_CF_LoadVector1,Memory_HWM)
    PRINT*,"Before First CF Load Vector         : ",Memory_Method_Before_CF_LoadVector1
ELSE IF ( iU == iU_LF ) THEN
    CALL Poseidon_Mark_Memory(Memory_Method_Before_LF_LoadVector1,Memory_HWM)
    PRINT*,"Before First LF Load Vector          : ",Memory_Method_Before_LF_LoadVector1
END IF
#endif

CALL XCFC_Calc_Load_Vector_TypeA( iU, iEU, iEL )

#ifdef POSEIDON_MEMORY_FLAG
IF ( iU == iU_CF ) THEN
    CALL Poseidon_Mark_Memory(Memory_Method_After_CF_LoadVector1,Memory_HWM)
    PRINT*,"After First CF Load Vector          : ",Memory_Method_After_CF_LoadVector1
ELSE IF ( iU == iU_LF ) THEN
    CALL Poseidon_Mark_Memory(Memory_Method_After_LF_LoadVector1,Memory_HWM)
    PRINT*,"After First LF Load Vector          : ",Memory_Method_After_LF_LoadVector1
END IF
#endif

Cur_Iteration = 0
CONVERGED     = .FALSE.
DO WHILE ( .NOT. CONVERGED  .AND. Cur_Iteration < Max_Iterations)
    Cur_Iteration = Cur_Iteration+1

    IF ( Verbose_Flag ) THEN
        WRITE(Message,'(A,A,A,I2.2,A)')'Starting ',TRIM(CFA_Var_Names(iU)),' Iteration ',Cur_Iteration,'.'
        CALL Run_Message(TRIM(Message))
    END IF

    mk = MIN(Cur_Iteration, M)

    IF ( Cur_Iteration .NE. 1 ) THEN
        UVector = GVectorM
    ELSE
        CALL Coeff_To_Vector_TypeA_SV( UVector, iU )
    END IF



#ifdef POSEIDON_MEMORY_FLAG
    IF ( Cur_Iteration == 1 ) THEN
        IF ( iU == iU_CF ) THEN
            CALL Poseidon_Mark_Memory(Memory_Method_Before_CF_Solve1,Memory_HWM)
            PRINT*,"Before First CF Liner Solve         : ",Memory_Method_Before_CF_Solve1
        ELSE IF ( iU == iU_LF ) THEN
            CALL Poseidon_Mark_Memory(Memory_Method_Before_LF_Solve1,Memory_HWM)
            PRINT*,"Before First LF Liner Solve         : ",Memory_Method_Before_LF_Solve1
        END IF
    ELSE IF ( Cur_Iteration == 2 ) THEN
        IF ( iU == iU_CF ) THEN
            CALL Poseidon_Mark_Memory(Memory_Method_Before_CF_Solve2,Memory_HWM)
            PRINT*,"Before Second CF Liner Solve        : ",Memory_Method_Before_CF_Solve2
        ELSE IF ( iU == iU_LF ) THEN
            CALL Poseidon_Mark_Memory(Memory_Method_Before_LF_Solve2,Memory_HWM)
            PRINT*,"Before Second LF Liner Solve        : ",Memory_Method_Before_LF_Solve2
        END IF
    ELSE IF ( Cur_Iteration == 3 ) THEN
        IF ( iU == iU_CF ) THEN
            CALL Poseidon_Mark_Memory(Memory_Method_Before_CF_Solve3,Memory_HWM)
            PRINT*,"Before Third CF Liner Solve         : ",Memory_Method_Before_CF_Solve3
        ELSE IF ( iU == iU_LF ) THEN
            CALL Poseidon_Mark_Memory(Memory_Method_Before_LF_Solve3,Memory_HWM)
            PRINT*,"Before Third LF Liner Solve         : ",Memory_Method_Before_LF_Solve3
        END IF
    END IF

#endif



    !
    !   Solve Systems
    !
    CALL Solve_Linear_System_TypeA(iU)
    
    
    
    
#ifdef POSEIDON_MEMORY_FLAG
    IF ( Cur_Iteration == 1 ) THEN
        IF ( iU == iU_CF ) THEN
            CALL Poseidon_Mark_Memory(Memory_Method_After_CF_Solve1,Memory_HWM)
            PRINT*,"After First CF Liner Solve          : ",Memory_Method_After_CF_Solve1
        ELSE IF ( iU == iU_LF ) THEN
            CALL Poseidon_Mark_Memory(Memory_Method_After_LF_Solve1,Memory_HWM)
            PRINT*,"After First LF Liner Solve          : ",Memory_Method_After_LF_Solve1
        END IF
    ELSE IF ( Cur_Iteration == 2 ) THEN
        IF ( iU == iU_CF ) THEN
            CALL Poseidon_Mark_Memory(Memory_Method_After_CF_Solve2,Memory_HWM)
            PRINT*,"After Second CF Liner Solve         : ",Memory_Method_After_CF_Solve2
        ELSE IF ( iU == iU_LF ) THEN
            CALL Poseidon_Mark_Memory(Memory_Method_After_LF_Solve2,Memory_HWM)
            PRINT*,"After Second LF Liner Solve         : ",Memory_Method_After_LF_Solve2
        END IF
    ELSE IF ( Cur_Iteration == 3 ) THEN
        IF ( iU == iU_CF ) THEN
            CALL Poseidon_Mark_Memory(Memory_Method_After_CF_Solve3,Memory_HWM)
            PRINT*,"After Third CF Liner Solve          : ",Memory_Method_After_CF_Solve3
        ELSE IF ( iU == iU_LF ) THEN
            CALL Poseidon_Mark_Memory(Memory_Method_After_LF_Solve3,Memory_HWM)
            PRINT*,"After Third LF Liner Solve          : ",Memory_Method_After_LF_Solve3
        END IF
    END IF
#endif
    
    
    
    CALL Coeff_To_Vector_TypeA_SV( GVector(:,mk), iU )
    
    FVector(:,mk) = GVector(:,mk) - UVector(:)


    IF ( mk == 1 ) THEN
        GVectorM = GVector(:,mk)
    ELSE

        BVector = -FVector(:,mk)
        AMatrix(:,1:mk-1) = FVector(:,1:mk-1) - SPREAD( FVector(:,mk), DIM=2, NCOPIES = mk-1)

        CALL DGELS('N',Var_Dim,mk-1,1,              &
                    AMatrix(:,1:mk-1), Var_Dim,     &
                    BVector, Var_Dim,               &
                    WORK, LWORK, INFO )

        IF ( INFO .NE. 0 ) THEN
            WRITE(Message,'(A,I1.1,A,I1.1)')'In XCFC_Fixed_Point, iU = ',iU,' : ZGELS failed with INFO = ',INFO
            CALL WARNING_MESSAGE(Message)
        END IF
        Alpha(1:mk-1) = BVector(1:mk-1)
        Alpha(mk)     = 1.0_idp - SUM( Alpha(1:mk-1) )

        GVectorM = 0.0_idp
        DO i = 1,mk
            GVectorM = GVectorM + Alpha(i)*GVector(:,i)
        END DO
    END IF

    FVectorM = GVectorM - UVector(:)


    
    ! Check for Convergence
    CALL Convergence_Check( FVectorM, Cur_Iteration, Converged )



    IF ( mk == M .AND. .NOT. Converged ) THEN
        GVector = CSHIFT( GVector, SHIFT = +1, DIM = 2)
        FVector = CSHIFT( FVector, SHIFT = +1, DIM = 2)
    END IF


    !   Update Coefficient Vector
    CALL Vector_To_Coeff_TypeA_SV( GVectorM, iU )

#ifdef POSEIDON_MEMORY_FLAG
    IF ( Cur_Iteration == 1 ) THEN
        IF ( iU == iU_CF ) THEN
            CALL Poseidon_Mark_Memory(Memory_Method_Before_CF_LoadVector2,Memory_HWM)
            PRINT*,"After Second CF Load Vector          : ",Memory_Method_Before_CF_LoadVector2
        ELSE IF ( iU == iU_LF ) THEN
            CALL Poseidon_Mark_Memory(Memory_Method_Before_LF_LoadVector2,Memory_HWM)
            PRINT*,"After Second LF Load Vector          : ",Memory_Method_Before_LF_LoadVector2
        END IF
    ELSE IF ( Cur_Iteration == 2 ) THEN
        IF ( iU == iU_CF ) THEN
            CALL Poseidon_Mark_Memory(Memory_Method_Before_CF_LoadVector3,Memory_HWM)
            PRINT*,"After Third CF Load Vector           : ",Memory_Method_Before_CF_LoadVector3
        ELSE IF ( iU == iU_LF ) THEN
            CALL Poseidon_Mark_Memory(Memory_Method_Before_LF_LoadVector3,Memory_HWM)
            PRINT*,"After Third LF Load Vector           : ",Memory_Method_Before_LF_LoadVector3
        END IF
    ELSE IF ( Cur_Iteration == 3 ) THEN
        IF ( iU == iU_CF ) THEN
            CALL Poseidon_Mark_Memory(Memory_Method_Before_CF_LoadVector4,Memory_HWM)
            PRINT*,"After Forth CF Load Vector           : ",Memory_Method_Before_CF_LoadVector4
        ELSE IF ( iU == iU_LF ) THEN
            CALL Poseidon_Mark_Memory(Memory_Method_Before_LF_LoadVector4,Memory_HWM)
            PRINT*,"After Forth LF Load Vector           : ",Memory_Method_Before_LF_LoadVector4
        END IF
    END IF
#endif


    !   Calculate Source Vector with updated solution
    CALL XCFC_Calc_Load_Vector_TypeA( iU, iEU, iEL )
    


    IF ( Verbose_Flag ) THEN
        WRITE(Message,'(A,A,A,I2.2,A)')'Ending ',TRIM(CFA_Var_Names(iU)),' Iteration ',Cur_Iteration,'.'
        CALL Run_Message(TRIM(Message))
        WRITE(*,'()')
    END IF


END DO ! Converged Loop


#ifdef POSEIDON_MEMORY_FLAG
IF ( iU == iU_CF ) THEN
    CALL Poseidon_Mark_Memory(Memory_Method_After_CF_FixedPoint,Memory_HWM)
    PRINT*,"After First Fixed Point Load Vector : ",Memory_Method_After_CF_FixedPoint
ELSE IF ( iU == iU_LF ) THEN
    CALL Poseidon_Mark_Memory(Memory_Method_After_LF_FixedPoint,Memory_HWM)
    PRINT*,"After First Fixed Point Load Vector : ",Memory_Method_After_LF_FixedPoint
END IF
#endif

DEALLOCATE( Work )


#ifdef POSEIDON_MEMORY_FLAG
IF ( iU == iU_CF ) THEN
    CALL Poseidon_Mark_Memory(Memory_Method_After_CF_DeallocWork,Memory_HWM)
    PRINT*,"After First Fixed Point Load Vector : ",Memory_Method_After_CF_DeallocWork
ELSE IF ( iU == iU_LF ) THEN
    CALL Poseidon_Mark_Memory(Memory_Method_After_LF_DeallocWork,Memory_HWM)
    PRINT*,"After First Fixed Point Load Vector : ",Memory_Method_After_LF_DeallocWork
END IF
#endif

END SUBROUTINE XCFC_Fixed_Point



















 !+201+########################################################!
!                                                               !
!           Convergence_Check                                   !
!                                                               !
 !#############################################################!
SUBROUTINE Convergence_Check( Update, Iter, Flag )


REAL(idp),DIMENSION(Var_Dim),   INTENT(IN)          :: Update
INTEGER,                        INTENT(IN)          :: Iter
LOGICAL,                        INTENT(INOUT)       :: Flag

CHARACTER(LEN = 300)                                ::  Message


IF ( Verbose_Flag ) THEN
    WRITE(Message,'(A,ES22.15,A)')'Maximum change in the coefficient vector : ',MAXVAL(ABS(Update)),'.'
    CALL Run_Message(TRIM(Message))
END IF



IF ( ALL( ABS( Update ) <= Convergence_Criteria ) ) THEN
    IF ( Verbose_Flag ) THEN
        CALL Run_Message("The Method has converged.")
        CALL Run_Message("The absolute update is less than the tolerance set. ")
    END IF
    Flag = .TRUE.
    
END IF


IF ( Iter == Max_Iterations ) THEN
    CALL Warning_Message('FP_Accelerated has reached the maximum number of allowed iterations.')
    Flag = .TRUE.
END IF




END SUBROUTINE Convergence_Check






END MODULE XCFC_Fixed_Point_Module
