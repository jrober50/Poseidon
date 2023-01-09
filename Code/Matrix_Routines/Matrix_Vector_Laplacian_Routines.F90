   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Matrix_Vector_Laplacian_Routines                                      !##!
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
            ONLY :  Init_Message,           &
                    Run_Message,            &
                    Warning_Message

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
                    iVB_Prob_Dim


USE Variables_IO, &
            ONLY :  File_Suffix

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature

USE Variables_Matrices, &
            ONLY :  Matrix_Format,              &
                    iMB_Diagonals,              &
                    iMB_Bandwidth,              &
                    iMB_IPIV,                   &
                    zMB_Matrix_Banded,          &
                    zMB_Matrix_Diagonal,        &
                    zMB_First_Col_Storage,      &
                    zMB_Last_Col_Storage

USE IO_Condition_Number_Output_Module, &
            ONLY :  IO_Output_Condition_Number

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
            ONLY :  Timer_Banded_Factorization


USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Matrices_Flags,    &
                    iPF_Init_Matrices_Type_B_LU


#ifdef POSEIDON_MEMORY_FLAG
USE Poseidon_Memory_Routines, &
            ONLY :  Poseidon_Mark_Memory

USE Memory_Variables_Module, &
            ONLY :  Memory_Method_Before_Bnd_Factorize,     &
                    Memory_Method_After_Bnd_Factorize
#endif
USE MPI

IMPLICIT NONE


CONTAINS



!+501+###########################################################################!
!                                                                                !
!           Factorize_Vector_Laplacian                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Factorize_Vector_Laplacian()



INTEGER                                                 :: INFO

CHARACTER(LEN = 300)                                    ::  Message

IF ( Verbose_Flag ) CALL Init_Message('Factorizing vector Laplace matrix using LAPAK LU factorization.')



#ifdef POSEIDON_MEMORY_FLAG
CALL Poseidon_Mark_Memory(Memory_Method_Before_Bnd_Factorize)
PRINT*,"Before Banded Factorization  : ",Memory_Method_Before_Bnd_Factorize
#endif




CALL TimerStart( Timer_Banded_Factorization )

!   Dirichlet BCs modify the stiffness matrix so we modify it now.
!   But to apply the BCs we will need values from the original matrix,
!   so those values are stored before we modify the matrix.
!
CALL Dirichlet_BC_Beta_Banded_Mat()




CALL Jacobi_PC_MVL_Banded()




CALL ZGBTRF( iVB_Prob_Dim,              &
             iVB_Prob_Dim,              &
             iMB_Diagonals,             &
             iMB_Diagonals,             &
             zMB_Matrix_Banded,           &
             3*iMB_Diagonals+1,         &
             iMB_IPIV,                  &
             INFO                       )

IF (INFO .NE. 0) THEN
    WRITE(Message,'(A,I1.1)')"ZGBTRF has failed with INFO = ",INFO
    CALL Warning_Message(TRIM(Message))
END IF





CALL TimerStop( Timer_Banded_Factorization )




CALL IO_Output_Condition_Number()


lPF_Init_Matrices_Flags(iPF_Init_Matrices_Type_B_LU) = .TRUE.


#ifdef POSEIDON_MEMORY_FLAG
CALL Poseidon_Mark_Memory(Memory_Method_After_Bnd_Factorize)
PRINT*,"After Banded Factorization   : ",Memory_Method_After_Bnd_Factorize
#endif


END SUBROUTINE Factorize_Vector_Laplacian







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



IF ( Verbose_Flag ) CALL Init_Message('Preconditioning the Modified Vector Laplacian matrix using Jacobi preconditioner.')

! Store inverse diagonal
DO i = 1,iVB_Prob_Dim

    Row = iMB_Bandwidth + i
    zMB_Matrix_Diagonal(i) = 1.0_idp/zMB_Matrix_Banded(iMB_Bandwidth,i)

END DO






! Row : RE = 0, D = 0
DO ui = 1,3
DO l  = 1,LM_Length

    Row = iMB_Bandwidth + FP_Beta_Array_Map(0,0,ui,l)


    ! Column : RE = 0, D = 0
    DO uj = 1,3
    DO lp = 1,LM_Length
    
        Col = FP_Beta_Array_Map(0,0,uj,lp)

        zMB_Matrix_Banded(Row-Col,Col) = zMB_Matrix_Banded(Row-Col,Col)*zMB_Matrix_Diagonal(Row-iMB_Bandwidth)

!        PRINT*,Row-iMB_Bandwidth,Col,zMB_Matrix_Banded(Row-Col,Col)

    END DO  ! lp
    END DO  ! uj


END DO ! ui
END DO ! l







DO RE = 0,Num_R_Elements-1

    ! Row : RE , D = 0
    DO ui = 1,3
    DO l  = 1,LM_Length

        Row = iMB_Bandwidth + FP_Beta_Array_Map(re, 0, ui, l)


        ! Column : RE , D = 1,DEGREE
        DO dp = 1,DEGREE
        DO uj = 1,3
        DO lp = 1,LM_Length

            Col = FP_Beta_Array_Map(RE,dp,uj,lp)

            zMB_Matrix_Banded(Row-Col,Col) = zMB_Matrix_Banded(Row-Col,Col)*zMB_Matrix_Diagonal(Row-iMB_Bandwidth)
!            PRINT*,Row-iMB_Bandwidth,Col,zMB_Matrix_Banded(Row-Col,Col)

        END DO  ! lp
        END DO  ! uj
        END DO  ! dp

    END DO ! l
    END DO ! ui




    ! Row RE D = 1,Degree
    DO d  = 1,DEGREE
    DO ui = 1,3
    DO l  = 1,LM_Length

        Row = iMB_Bandwidth + FP_Beta_Array_Map(RE,d,ui,l)


        ! Column RE, D = 0
        DO uj = 1,3
        DO lp = 1,LM_Length

            Col = FP_Beta_Array_Map(RE,0,uj,lp)

            zMB_Matrix_Banded(Row-Col,Col) = zMB_Matrix_Banded(Row-Col,Col)*zMB_Matrix_Diagonal(Row-iMB_Bandwidth)
!            PRINT*,Row-iMB_Bandwidth,Col,zMB_Matrix_Banded(Row-Col,Col)


        END DO  ! lp
        END DO  ! uj





        ! Column RE, D = 1,DEGREE
        DO dp = 1,DEGREE
        DO uj = 1,3
        DO lp = 1,LM_Length

            Col = FP_Beta_Array_Map(RE,dp,uj,lp)

            zMB_Matrix_Banded(Row-Col,Col) = zMB_Matrix_Banded(Row-Col,Col)*zMB_Matrix_Diagonal(Row-iMB_Bandwidth)
!            PRINT*,Row-iMB_Bandwidth,Col,zMB_Matrix_Banded(Row-Col,Col)


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

COMPLEX(idp),   DIMENSION(1:iVB_Prob_Dim),     INTENT(INOUT)       ::  Work_Vec

IF ( Verbose_Flag ) CALL Init_Message('Preconditioning the RHS vector using Jacobi preconditioner.')

! Multiply diagonal and work vec

Work_Vec(:)  = Work_Vec(:)*zMB_Matrix_Diagonal(:)


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

            Row = iMB_Bandwidth                    &
                + FP_Beta_Array_Map(0,d,ui,lm)

            zMB_First_Col_Storage(lm,d,ui) = zMB_Matrix_Banded(Row-Col,Col)
   
        END DO ! l Loop
        END DO ! d Loop



       DO d = 0,Degree
       DO lm = 1,LM_Length

           Row = iMB_Bandwidth                    &
               + FP_Beta_Array_Map(0,d,ui,lm)


           zMB_Matrix_Banded(Row-Col,Col) = 0.0_idp

       END DO ! l Loop
       END DO ! d Loop




       DO d = 0,Degree
       DO lm = 1,LM_Length

           Row = iMB_Bandwidth + FP_Beta_Array_Map(0,0,ui,0)
           Col = FP_Beta_Array_Map(0,d,ui,lm)

           zMB_Matrix_Banded(Row-Col,Col) = 0.0_idp

       END DO ! l Loop
       END DO ! d Loop

       

       DO lm = 1,LM_Length

           Row = FP_Beta_Array_Map(0,0,ui,lm)
           Col = FP_Beta_Array_Map(0,0,ui,lm)
 
           zMB_Matrix_Banded(Row-Col,Col) = 1.0_idp

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

            Row = FP_Beta_Array_Map(Num_R_Elements-1,d,ui,lm)+iMB_Bandwidth
            
            zMB_Last_Col_Storage(lm,d,ui) = zMB_Matrix_Banded(Row-Col,Col)
    
        END DO ! l Loop
        END DO ! d Loop



        ! Clear the Rows
        DO lm = 1,LM_Length

            Row = FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui,lm) + iMB_Bandwidth

            DO uj = 1,3
            DO dp = 0,Degree
            DO lp = 1,LM_Length

                Col = FP_Beta_Array_Map(Num_R_Elements-1,dp,uj,lp)

                zMB_Matrix_Banded(Row-Col,Col) = 0.0_idp
            

            END DO ! lp
            END DO ! dp
            END DO ! uj
        END DO ! l Loop





        ! Clear the Columns
        DO lm = 1,LM_Length

            Col = FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui,lm)

            zMB_Matrix_Banded(:,Col) = 0.0_idp

        END DO ! l Loop




        
        DO lm = 1,LM_Length

            Row = iMB_Bandwidth + FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui,lm)
            Col = FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui,lm)

            zMB_Matrix_Banded(Row-Col,Col) = 1.0_idp

        END DO ! lm Loop


    END IF

END DO


END SUBROUTINE DIRICHLET_BC_Beta_Banded_Mat



END MODULE Matrix_Vector_Laplacian_Routines
