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
                    L_LIMIT

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
                    CFA_Mat_Map


USE FP_Functions_Laplace_Beta, &
            ONLY :  Initialize_Laplace_Matrices_Beta, &
                    Output_Laplace


IMPLICIT NONE

CONTAINS


!+101+###########################################################################!
!                                                                                !
!           Initialize_Laplace_Matrices                                          !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_Laplace_Matrices()

LOGICAL                                 :: Success_Flag

Success_Flag = .FALSE.
IF ( Matrix_Format == 'Full' ) THEN
    
    CALL Initialize_Laplace_Matrices_Full()
    CALL Initialize_Laplace_Matrices_Beta()
    PRINT*,"Poseidon Initialized the Laplace Matrices. Format : ",Matrix_Format
    Success_Flag = .TRUE.
ELSEIF ( Matrix_Format == 'CCS' ) THEN

    CALL Initialize_Laplace_Matrices_CCS()
    PRINT*,"Poseidon Initialized the Laplace Matrices. Format : ",Matrix_Format
    Success_Flag = .TRUE.
END IF





IF ( Success_Flag .EQV. .FALSE.) THEN
    PRINT*,"WARNING: Poseidon did not initalize the Laplace matrices."
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

INTEGER                                                 ::  Mat_Loc

INTEGER                                                 ::  Int_Degree
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  CUR_R_LOCS
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  R_SQUARE
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  Int_Locs
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  Int_Weights

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
                i = re*DEGREE+d+1
                j = re*DEGREE+dp+1

!                Laplace_Matrix_Full(i, j, l) = Laplace_Matrix_Full(i, j, l)           &
!                                             + SUM( R_SQUARE(:) * LPT_LPT(:,d,dp,1,1) &
!                                                    * TODR * Int_Weights(:)           &
!                                                   - L_Lp1 * LPT_LPT(:,d,dp,0,0)      &
!                                                    * Int_Weights(:)                )
                DO rd = 1,Int_Degree


                    
                    Laplace_Matrix_Full(i, j, l,:) = Laplace_Matrix_Full(i, j, l,:) &
                                        - R_SQUARE(rd) * LPT_LPT(rd,d,dp,1,1)       &
                                            * TODR * Int_Weights(rd)                &
                                        + L_Lp1 * LPT_LPT(rd,d,dp,0,0)              &
                                            * Int_Weights(rd)
                    
                END DO  ! rd Loop

            END DO  ! dp Loop
        END DO  ! d Loop
    END DO  ! re Loop
END DO  ! l Loop


Call Output_Laplace(Laplace_Matrix_Full(:,:, 0,CFA_Mat_Map(3)), Num_R_Nodes, Num_R_Nodes, "W")


IF ( CFA_EQ_Flags(3) == 1 ) THEN

    Mat_Loc = CFA_Mat_Map(3)
    DO l = 0,L_LIMIT
        DO re = 0,NUM_R_ELEMENTS-1

            DO dp = 0,DEGREE
                DO d = 0,DEGREE
                    i = re*DEGREE+d+1
                    j = re*DEGREE+dp+1

                    DO rd = 1,Int_Degree


                        
                        Laplace_Matrix_Full(i, j, l,Mat_Loc)             &
                                            = Laplace_Matrix_Full(i, j, l,Mat_Loc)    &
                                            - 2 * LPT_LPT(rd,d,dp,0,0)                  &
                                                / TODR * Int_Weights(rd)

!                        PRINT*,LPT_LPT(rd,d,dp,0,0)                  &
!                        / TODR * Int_Weights(rd),   re, rd, d, dp
                        
                    END DO  ! rd Loop

                END DO  ! dp Loop
            END DO  ! d Loop
        END DO  ! re Loop
    END DO  ! l Loop

Call Output_Laplace(Laplace_Matrix_Full(:,:, 0,CFA_Mat_Map(3)), Num_R_Nodes, Num_R_Nodes, "X")
   
END IF



END SUBROUTINE Initialize_Laplace_Matrices_Full





!+202+###########################################################################!
!                                                                                !
!           Initialize_Laplace_Matrices_CCS                                      !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_Laplace_Matrices_CCS()


INTEGER                                                 ::  l, re, rd, d, dp
INTEGER                                                 ::  i, j, Here


REAL(KIND = idp)                                        ::  deltar, TODR
REAL(KIND = idp)                                        ::  L_Lp1

INTEGER                                                 ::  Int_Degree
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  CUR_R_LOCS
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  R_SQUARE
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  Int_Locs
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  Int_Weights

Int_Degree = 10

Laplace_Matrix_Full = 0.0_idp

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

            Laplace_Matrix_ROW(Here,:) = i*DEGREE + dp
            Here = Here + 1

        END DO ! dp Loop
    END DO ! d Loop

    DO d = 0,DEGREE - 1

        Laplace_Matrix_ROW(Here,:) = i*DEGREE + d
        Here = Here + 1

    END DO ! d Loop
END DO ! re Loop
Laplace_Matrix_ROW(Here,:) = DEGREE * NUM_R_ELEMENTS



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
                                        + R_SQUARE(rd) * LPT_LPT(rd,d,dp,1,1)   &
                                            * TODR * Int_Weights(rd)            &
                                        - L_Lp1 * LPT_LPT(rd,d,dp,0,0)          &
                                            * Int_Weights(rd)
                                            
                END DO  ! rd Loop
                Here = Here + 1

            END DO  ! dp Loop
        END DO  ! d Loop
        Here = Here - 1 !!! Take one step back due to overlap of first value
                        !!! in next element with last in current element
    END DO  ! re Loop
END DO  ! l Loop


Laplace_Factored_VAL = Laplace_Matrix_VAL
Laplace_Factored_ROW = Laplace_Matrix_ROW
Laplace_Factored_COL = Laplace_Matrix_COL


END SUBROUTINE Initialize_Laplace_Matrices_CCS










END MODULE FP_Functions_Laplace
