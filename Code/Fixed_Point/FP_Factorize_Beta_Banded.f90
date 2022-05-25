   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE FP_Factorize_Beta_Banded                                      !##!
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

USE Poseidon_Parameters, &
            ONLY :  DEGREE,                     &
                    L_LIMIT,                    &
                    Verbose_Flag

USE Poseidon_IO_Parameters, &
            ONLY :  Poseidon_LinSys_Dir

    
USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,          &
                    Num_T_Quad_Points,          &
                    Num_P_Quad_Points,          &
                    Num_TP_Quad_Points,         &
                    Int_R_Locations,            &
                    Int_T_Locations,            &
                    Int_P_Locations,            &
                    Int_R_Weights,              &
                    Int_T_Weights,              &
                    Int_P_Weights


USE Variables_Mesh, &
            ONLY :  Num_R_Elements,             &
                    Num_T_Elements,             &
                    Num_P_Elements,             &
                    rlocs,                      &
                    tlocs,                      &
                    plocs

USE Variables_Derived, &
            ONLY :  Num_R_Nodes,                &
                    ULM_Length,                 &
                    LM_Length,                  &
                    Var_Dim,                    &
                    Beta_Prob_Dim


USE Variables_IO, &
            ONLY :  File_Suffix

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature

USE Variables_FP, &
            ONLY :  Matrix_Format,              &
                    Beta_Diagonals,             &
                    Beta_Bandwidth,             &
                    Beta_IPIV,                  &
                    Beta_MVL_Banded,            &
                    Beta_MVL_Diagonal,          &
                    First_Column_Beta_Storage,  &
                    Last_Column_Beta_Storage,   &
                    Beta_Factorized_Flag

USE Poseidon_IO_Module, &
            ONLY :  CLOCK_IN

USE IO_FP_Linear_System, &
            ONLY :  Output_Laplace,             &
                    Output_Laplace_Beta,        &
                    Output_Source_Beta

USE Variables_BC, &
            ONLY :  INNER_CFA_BC_VALUES,        &
                    OUTER_CFA_BC_VALUES,        &
                    INNER_CFA_BC_TYPE,          &
                    OUTER_CFA_BC_TYPE

USE Maps_Fixed_Point, &
            ONLY :  FP_Beta_Array_Map

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                 &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_XCFC_Banded_Factorization


USE MPI

IMPLICIT NONE


CONTAINS



!+501+###########################################################################!
!                                                                                !
!                   Factorize_Beta_Banded                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Factorize_Beta_Banded()



INTEGER                                                 :: INFO
COMPLEX(idp),   DIMENSION(2*Beta_Prob_Dim)              :: WORK
REAL(idp),      DIMENSION(Beta_Prob_Dim)                :: RWORK
REAL(idp)                                               :: RCOND_One
REAL(idp)                                               :: NORM

INTEGER                                                 :: i


IF ( Verbose_Flag ) THEN
    PRINT*,"--In Factorize_Beta_Banded."
END IF
CALL TimerStart( Timer_XCFC_Banded_Factorization )

!   Dirichlet BCs modify the stiffness matrix so we modify it now.
!   But to apply the BCs we will need values from the original matrix,
!   so those values are stored before we modify the matrix.
!

!PRINT*,"Beta_MVL_Banded in Factorize_Beta_Banded"
!PRINT*,Beta_MVL_Banded


CALL DIRICHLET_BC_Beta_Banded_Mat()


CALL Jacobi_PC_MVL_Banded()




CALL ZGBTRF( Beta_Prob_Dim,             &
             Beta_Prob_Dim,             &
             Beta_Diagonals,            &
             Beta_Diagonals,            &
             Beta_MVL_Banded,           &
             3*Beta_Diagonals+1,        &
             Beta_IPIV,                 &
             INFO                       )

IF (INFO .NE. 0) THEN
    print*,"ZGBTRF has failed with INFO = ",INFO
ELSE

    Beta_Factorized_Flag = .TRUE.

END IF



!PRINT*,"Beta_MVL_Banded in Factorize_Beta_Banded"
!PRINT*,REAL(Beta_MVL_Banded, kind = idp)






IF ( Verbose_Flag ) THEN

    NORM = 0.0_idp
    DO i = 1,Beta_Prob_Dim

        NORM = MAX( NORM, ABS(SUM(Beta_MVL_Banded(:,i) ) ) )

    END DO

    CALL ZGBCON( '1',                   &
                 Beta_Prob_Dim,         &
                 Beta_Diagonals,        &
                 Beta_Diagonals,        &
                 Beta_MVL_Banded,       &
                 3*Beta_Diagonals+1,    &
                 Beta_IPIV,             &
                 NORM,                  &
                 RCOND_One,             &
                 WORK,                  &
                 RWORK,                 &
                 INFO                   )
    IF (INFO .NE. 0) THEN
        print*,"ZGBCON has failed with INFO = ",INFO
    ELSE
            PRINT*,"RCOND = ",RCOND_One," Norm = ",Norm
    END IF
END IF


CALL TimerStop( Timer_XCFC_Banded_Factorization )



!PRINT*,"STOPing in Factorize_Beta_Banded"
!STOP

END SUBROUTINE Factorize_Beta_Banded







!+501+###########################################################################!
!                                                                                !
!                   Jacobi_PC_MVL_Banded                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Jacobi_PC_MVL_Banded()

INTEGER                                         ::  i
INTEGER                                         ::  Row, Col
INTEGER                                         ::  RE
INTEGER                                         ::  ui, d, l
INTEGER                                         ::  uj, dp, lp

! Store inverse diagonal
DO i = 1,Beta_Prob_Dim

    Row = Beta_Bandwidth + i
    Beta_MVL_Diagonal(i) = 1.0_idp/Beta_MVL_Banded(Beta_Bandwidth,i)

END DO






! Row : RE = 0, D = 0
DO ui = 1,3
DO l  = 1,LM_Length

    Row = Beta_Bandwidth + FP_Beta_Array_Map(0,0,ui,l)


    ! Column : RE = 0, D = 0
    DO uj = 1,3
    DO lp = 1,LM_Length
    
        Col = FP_Beta_Array_Map(0,0,uj,lp)

        Beta_MVL_Banded(Row-Col,Col) = Beta_MVL_Banded(Row-Col,Col)*Beta_MVL_Diagonal(Row-Beta_Bandwidth)

!        PRINT*,Row-Beta_Bandwidth,Col,Beta_MVL_Banded(Row-Col,Col)

    END DO  ! lp
    END DO  ! uj


END DO ! ui
END DO ! l







DO RE = 0,Num_R_Elements-1

    ! Row : RE , D = 0
    DO ui = 1,3
    DO l  = 1,LM_Length

        Row = Beta_Bandwidth + FP_Beta_Array_Map(re, 0, ui, l)


        ! Column : RE , D = 1,DEGREE
        DO dp = 1,DEGREE
        DO uj = 1,3
        DO lp = 1,LM_Length

            Col = FP_Beta_Array_Map(RE,dp,uj,lp)

            Beta_MVL_Banded(Row-Col,Col) = Beta_MVL_Banded(Row-Col,Col)*Beta_MVL_Diagonal(Row-Beta_Bandwidth)
!            PRINT*,Row-Beta_Bandwidth,Col,Beta_MVL_Banded(Row-Col,Col)

        END DO  ! lp
        END DO  ! uj
        END DO  ! dp

    END DO ! l
    END DO ! ui




    ! Row RE D = 1,Degree
    DO d  = 1,DEGREE
    DO ui = 1,3
    DO l  = 1,LM_Length

        Row = Beta_Bandwidth + FP_Beta_Array_Map(RE,d,ui,l)


        ! Column RE, D = 0
        DO uj = 1,3
        DO lp = 1,LM_Length

            Col = FP_Beta_Array_Map(RE,0,uj,lp)

            Beta_MVL_Banded(Row-Col,Col) = Beta_MVL_Banded(Row-Col,Col)*Beta_MVL_Diagonal(Row-Beta_Bandwidth)
!            PRINT*,Row-Beta_Bandwidth,Col,Beta_MVL_Banded(Row-Col,Col)


        END DO  ! lp
        END DO  ! uj





        ! Column RE, D = 1,DEGREE
        DO dp = 1,DEGREE
        DO uj = 1,3
        DO lp = 1,LM_Length

            Col = FP_Beta_Array_Map(RE,dp,uj,lp)

            Beta_MVL_Banded(Row-Col,Col) = Beta_MVL_Banded(Row-Col,Col)*Beta_MVL_Diagonal(Row-Beta_Bandwidth)
!            PRINT*,Row-Beta_Bandwidth,Col,Beta_MVL_Banded(Row-Col,Col)


        END DO  ! lp
        END DO  ! uj
        END DO  ! dp



    END DO  ! l
    END DO  ! ui
    END DO  ! d


END DO ! RE





END SUBROUTINE Jacobi_PC_MVL_Banded









!+501+###########################################################################!
!                                                                                !
!                   Jacobi_PC_MVL_Banded                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Jacobi_PC_MVL_Banded_Vector( Work_Vec )

COMPLEX(idp),   DIMENSION(1:Beta_Prob_Dim),     INTENT(INOUT)       ::  Work_Vec



! Multiply diagonal and work vec

Work_Vec(:)  = Work_Vec(:)*Beta_MVL_Diagonal(:)


END SUBROUTINE Jacobi_PC_MVL_Banded_Vector








!+501+###########################################################################!
!                                                                                !
!                   DIRICHLET_BC_Beta_Banded_Mat                                       !
!                                                                                !
!################################################################################!
SUBROUTINE DIRICHLET_BC_Beta_Banded_Mat()


INTEGER                                                 :: Col, Row
INTEGER                                                 :: d, ui, lm
INTEGER                                                 :: dp, uj, lp



DO ui = 1,3

    !############################################!
    !#                                          #!
    !#          Inner Boundary Condition        #!
    !#                                          #!
    !############################################!
    IF (INNER_CFA_BC_TYPE(ui) == "D") THEN

        Col = FP_Beta_Array_Map(0,0,ui,0)
    

        DO d = 0,Degree
        DO lm = 1,LM_Length

            Row = Beta_Bandwidth                    &
                + FP_Beta_Array_Map(0,d,ui,lm)

            First_Column_Beta_Storage(lm,d,ui) = Beta_MVL_Banded(Row-Col,Col)
   
        END DO ! l Loop
        END DO ! d Loop



       DO d = 0,Degree
       DO lm = 1,LM_Length

           Row = Beta_Bandwidth                    &
               + FP_Beta_Array_Map(0,d,ui,lm)


           Beta_MVL_Banded(Row-Col,Col) = 0.0_idp

       END DO ! l Loop
       END DO ! d Loop




       DO d = 0,Degree
       DO lm = 1,LM_Length

           Row = Beta_Bandwidth + FP_Beta_Array_Map(0,0,ui,0)
           Col = FP_Beta_Array_Map(0,d,ui,lm)

           Beta_MVL_Banded(Row-Col,Col) = 0.0_idp

       END DO ! l Loop
       END DO ! d Loop

       

       DO lm = 1,LM_Length

           Row = FP_Beta_Array_Map(0,0,ui,lm)
           Col = FP_Beta_Array_Map(0,0,ui,lm)
 
           Beta_MVL_Banded(Row-Col,Col) = 1.0_idp

       END DO ! l Loop
        
    END IF






    !############################################!
    !#                                          #!
    !#          Outer Boundary Condition        #!
    !#                                          #!
    !############################################!
    IF (OUTER_CFA_BC_TYPE(ui) == "D") THEN

        !
        ! Save the needed values first
        !


        Col = FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui,1)
        DO d = 0,Degree
        DO lm = 1,LM_Length

            Row = FP_Beta_Array_Map(Num_R_Elements-1,d,ui,lm)+Beta_Bandwidth
            
            Last_Column_Beta_Storage(lm,d,ui) = Beta_MVL_Banded(Row-Col,Col)
    
        END DO ! l Loop
        END DO ! d Loop



        ! Clear the Rows
        DO lm = 1,LM_Length

            Row = FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui,lm) + Beta_Bandwidth

            DO uj = 1,3
            DO dp = 0,Degree
            DO lp = 1,LM_Length

                Col = FP_Beta_Array_Map(Num_R_Elements-1,dp,uj,lp)

                Beta_MVL_Banded(Row-Col,Col) = 0.0_idp
            

            END DO ! lp
            END DO ! dp
            END DO ! uj
        END DO ! l Loop





        ! Clear the Columns
        DO lm = 1,LM_Length

            Col = FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui,lm)

            Beta_MVL_Banded(:,Col) = 0.0_idp

        END DO ! l Loop




        
        DO lm = 1,LM_Length

            Row = Beta_Bandwidth + FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui,lm)
            Col = FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui,lm)

            Beta_MVL_Banded(Row-Col,Col) = 1.0_idp

        END DO ! lm Loop


    END IF

END DO


END SUBROUTINE DIRICHLET_BC_Beta_Banded_Mat



END MODULE FP_Factorize_Beta_Banded
