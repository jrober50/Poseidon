   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE FP_Functions_Laplace                                                         !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
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
                    L_LIMIT,                    &
                    Verbose_Flag

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,             &
                    rlocs

USE Variables_Derived, &
            ONLY :  Num_R_Nodes,                &
                    ULM_Length,                 &
                    LM_Length,                  &
                    Var_Dim

USE Variables_Tables, &
            ONLY :  LPT_LPT

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature

USE Variables_FP, &
            ONLY :  Matrix_Format,              &
                    Laplace_Matrix_Full,        &
                    Laplace_Matrix_VAL,         &
                    Laplace_Matrix_ROW,         &
                    Laplace_Matrix_COL,         &
                    Laplace_Factored_VAL,       &
                    Laplace_Factored_ROW,       &
                    Laplace_Factored_COL,       &
                    CFA_EQ_Flags,               &
                    CFA_EQ_Map,                 &
                    CFA_Mat_Map,                &
                    MCF_Flag


USE Poseidon_Cholesky_Module,   &
            ONLY :  Cholesky_Factorization

USE FP_Functions_Laplace_Beta, &
            ONLY :  Initialize_Laplace_Matrices_Beta

USE FP_Beta_Banded, &
            ONLY :  Initialize_Beta_MVL_Banded

USE IO_FP_Linear_System, &
            ONLY :  Output_Laplace

USE Poseidon_IO_Module, &
            ONLY :  Clock_In

USE FP_Functions_Mapping, &
            ONLY :  FP_FEM_Node_Map

USE MPI

IMPLICIT NONE

CONTAINS


!+101+###########################################################################!
!                                                                                !
!           Initialize_Laplace_Matrices                                          !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_Laplace_Matrices()

LOGICAL                                 ::  Success_Flag
REAL(idp), DIMENSION(1:3)               ::  Timer


Success_Flag = .FALSE.
IF ( Matrix_Format == 'Full' ) THEN
    
    
    IF ( ANY(CFA_EQ_Flags(1:2) == 1 ) ) THEN
        timer(1) = MPI_WTime()
        CALL Initialize_Laplace_Matrices_Full()
        timer(2) = MPI_Wtime()
        Call Clock_In(timer(2)-timer(1),1)
    END IF
    IF ( ANY(CFA_EQ_Flags(3:5) == 1) ) THEN
        timer(1) = MPI_WTime()
        CALL Initialize_Laplace_Matrices_Beta()
        timer(2) = MPI_Wtime()
        Call Clock_In(timer(2)-timer(1),2)
        
    END IF

    Success_Flag = .TRUE.



ELSEIF ( Matrix_Format == 'CCS' ) THEN


    IF ( ANY(CFA_EQ_Flags(1:2) == 1 ) ) THEN
        timer(1) = MPI_WTime()
        CALL Initialize_Laplace_Matrices_CCS()
        timer(2) = MPI_Wtime()
        Call Clock_In(timer(2)-timer(1),1)
    
    END IF
    IF ( ANY(CFA_EQ_Flags(3:5) == 1) ) THEN
        timer(1) = MPI_WTime()
        CALL Initialize_Beta_MVL_Banded()
        timer(2) = MPI_Wtime()
        Call Clock_In(timer(2)-timer(1),2)
    END IF

    Success_Flag = .TRUE.
END IF





IF ( Success_Flag .EQV. .FALSE.) THEN
    PRINT*,"WARNING: Poseidon failed to initalize the Laplace matrices."
    PRINT*,"         Matrix_Format = ",Matrix_Format
    STOP
END IF


END SUBROUTINE Initialize_Laplace_Matrices




!+201+###########################################################################!
!                                                                                !
!           Initialize_Laplace_Matrices_Full                                     !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_Laplace_Matrices_Full()

INTEGER                                                 ::  l, re, rd, d, dp
INTEGER                                                 ::  i, j


REAL(KIND = idp)                                        ::  deltar, TODR
REAL(KIND = idp)                                        ::  L_Lp1


INTEGER                                                 ::  Int_Degree
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  CUR_R_LOCS
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  R_SQUARE
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  Int_Locs
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  Int_Weights



IF ( Verbose_Flag ) THEN
    PRINT*,"-Initializing Laplace Matrix.  Format: Full. "
END IF


Int_Degree = NUM_R_QUAD_POINTS

Laplace_Matrix_Full = 0.0_idp

ALLOCATE( CUR_R_LOCS(1:Int_Degree)  )
ALLOCATE( R_SQUARE(1:Int_Degree)    )
ALLOCATE( Int_Locs(1:Int_Degree)    )
ALLOCATE( Int_Weights(1:Int_Degree) )


CALL Initialize_LG_Quadrature(Int_Degree, Int_Locs, Int_Weights )


DO l = 0,L_LIMIT
    L_Lp1 = REAL( l*(l+1), idp )
    DO re = 0,NUM_R_ELEMENTS-1

        deltar = rlocs(re+1) - rlocs(re)
        
        TODR = 2.0_idp/deltar
        CUR_R_LOCS(:) = (deltar/2.0_idp) * (Int_Locs(:)+1.0_idp) + rlocs(re)
        R_SQUARE = CUR_R_LOCS**2

        DO dp = 0,DEGREE
            DO d = 0,DEGREE
                i = FP_FEM_Node_Map(re,d)
                j = FP_FEM_Node_Map(re,dp)

!                Laplace_Matrix_Full(i, j, l) = Laplace_Matrix_Full(i, j, l)           &
!                                             + SUM( R_SQUARE(:) * LPT_LPT(:,d,dp,1,1) &
!                                                    * TODR * Int_Weights(:)           &
!                                                   - L_Lp1 * LPT_LPT(:,d,dp,0,0)      &
!                                                    * Int_Weights(:)                )
                DO rd = 1,Int_Degree


                    
                    Laplace_Matrix_Full(i, j, l) = Laplace_Matrix_Full(i, j, l) &
                                        - R_SQUARE(rd) * LPT_LPT(rd,d,dp,1,1)       &
                                            * TODR * Int_Weights(rd)                &
                                        + L_Lp1 * LPT_LPT(rd,d,dp,0,0)              &
                                            * Int_Weights(rd)
                    
                END DO  ! rd Loop

            END DO  ! dp Loop
        END DO  ! d Loop
    END DO  ! re Loop
END DO  ! l Loop


!DO l = 0,L_LIMIT
!    PRINT*,Laplace_Matrix_Full(:, :, l)
!    PRINT*," "
!END DO

!Call Output_Laplace(Laplace_Matrix_Full(:,:,0), Num_R_Nodes, Num_R_Nodes, "W")


END SUBROUTINE Initialize_Laplace_Matrices_Full





!+202+###########################################################################!
!                                                                                !
!           Initialize_Laplace_Matrices_CCS                                      !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_Laplace_Matrices_CCS()


INTEGER                                                 ::  l, re, rd, d, dp
INTEGER                                                 ::  Here


REAL(KIND = idp)                                        ::  deltar, TODR
REAL(KIND = idp)                                        ::  L_Lp1

INTEGER                                                 ::  Int_Degree
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  CUR_R_LOCS
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  R_SQUARE
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  Int_Locs
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  Int_Weights


IF ( Verbose_Flag ) THEN
    PRINT*,"-Initializing Laplace Matrix.  Format: CCS. "
END IF


Int_Degree = NUM_R_QUAD_POINTS

ALLOCATE( CUR_R_LOCS(1:Int_Degree)  )
ALLOCATE( R_SQUARE(1:Int_Degree)    )
ALLOCATE( Int_Locs(1:Int_Degree)    )
ALLOCATE( Int_Weights(1:Int_Degree) )

CALL Initialize_LG_Quadrature(Int_Degree, Int_Locs, Int_Weights )


Laplace_Matrix_COL(0,:) = 0
Laplace_Matrix_COL(1,:) = Laplace_Matrix_COL(0,:) + (DEGREE+1)
HERE = 2


DO re = 1,NUM_R_ELEMENTS-1
    DO d = 1,DEGREE - 1

        Laplace_Matrix_COL(Here,:) = Laplace_Matrix_COL(Here - 1,:) +  (DEGREE + 1)
        Here = Here + 1
    END DO

    Laplace_Matrix_COL(Here,:) = Laplace_Matrix_COL(Here - 1,:) + (2*DEGREE + 1)
    Here = Here + 1

END DO

DO d = 1, DEGREE

    Laplace_Matrix_COL(Here,:) = Laplace_Matrix_COL(Here - 1,:) + (DEGREE + 1)
    Here = Here + 1

END DO



  !                             !
 !!                             !!
!!!    ROW_IND INITIALIZATION   !!!
 !!                             !!
  !                             !
Here = 0
DO re = 0, NUM_R_ELEMENTS - 1
    DO d = 0,DEGREE - 1
    DO dp = 0, DEGREE

        Laplace_Matrix_ROW(Here,:) = re*DEGREE + dp
        Here = Here + 1

    END DO ! dp Loop
    END DO ! d Loop

    DO d = 0,DEGREE - 1

        Laplace_Matrix_ROW(Here,:) = re*DEGREE + d
        Here = Here + 1

    END DO ! d Loop
END DO ! re Loop
Laplace_Matrix_ROW(Here,:) = DEGREE * NUM_R_ELEMENTS



Laplace_Matrix_VAL = 0.0_idp

DO l = 0,L_LIMIT
    L_Lp1 = REAL( l*(l+1), idp )
    Here = 0
    DO re = 0,NUM_R_ELEMENTS-1

        deltar = rlocs(re+1) - rlocs(re)
        TODR = 2.0_idp/deltar
        CUR_R_LOCS(:) = (deltar/2.0_idp) * (Int_Locs(:)+1.0_idp) + rlocs(re)
        R_SQUARE = CUR_R_LOCS**2

        DO dp = 0,DEGREE
            DO d = 0,DEGREE


!                Laplace_Matrix_VAL(Here,l) = Laplace_Matrix_VAL(Here,l)   &
!                                + SUM( R_SQUARE(:) * LPT_LPT(:,d,dp,1,1)  &
!                                        * TODR * Int_Weights(:)           &
!                                       - L_Lp1 * LPT_LPT(:,d,dp,0,0)      &
!                                        * Int_Weights(:)                )
                DO rd = 1,Int_Degree


                    
                    Laplace_Matrix_VAL(Here,l,:) = Laplace_Matrix_VAL(Here,l,:) &
                                        - R_SQUARE(rd) * LPT_LPT(rd,d,dp,1,1)   &
                                            * TODR * Int_Weights(rd)            &
                                        + L_Lp1 * LPT_LPT(rd,d,dp,0,0)          &
                                            * Int_Weights(rd)
                                            
                END DO  ! rd Loop
                Here = Here + 1

            END DO  ! dp Loop
        END DO  ! d Loop
        Here = Here - 1 !!! Take one step back due to overlap of first value
                        !!! in next element with last in current element
    END DO  ! re Loop
END DO  ! l Loop


!DO l = 0,L_LIMIT
!    PRINT*,Laplace_Matrix_VAL(:,l,1)
!    PRINT*," "
!END DO

MCF_Flag = 0
Laplace_Factored_VAL = -Laplace_Matrix_VAL
Laplace_Factored_ROW = Laplace_Matrix_ROW
Laplace_Factored_COL = Laplace_Matrix_COL


!CALL Cholesky_Factorization()
!MCF_Flag = 1



!PRINT*,"Row"
!PRINT*,-Laplace_Matrix_VAL
!PRINT*,"col"
!PRINT*,Laplace_Matrix_ROW
!PRiNT*,"Val"
!PRINT*,Laplace_Matrix_COL


END SUBROUTINE Initialize_Laplace_Matrices_CCS










END MODULE FP_Functions_Laplace
