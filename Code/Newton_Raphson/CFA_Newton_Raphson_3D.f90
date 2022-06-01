   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE CFA_Newton_Raphson_3D_Module                                                 !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!        Contains the functions used to calculate the components of the          !##!
!##!    extended Jacobian matrix as well as the derivative coefficients. These      !##!
!##!    are used to construct the stiffness matrix.                                 !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   CFA_Newton_Raphson                                                  !##!
!##!                                                                                !##!
!##!    +201+   CFA_Solve                                                           !##!
!##!                                                                                !##!
!##!    +301+   CFA_Coefficient_Update_All                                          !##!
!##!    +302+   CFA_Coefficient_Update                                              !##!
!##!                                                                                !##!
!##!    +401+   CFA_Coefficient_Share                                               !##!
!##!    +402+   CFA_Update_Share                                                    !##!
!##!                                                                                !##!
!##!    +501+   Output_Iteration_Report                                             !##!
!##!                                                                                !##!
!##!    +601+   Convergence_Check                                                   !##!
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

USE Poseidon_Numbers_Module, &
                    ONLY : pi

USE Poseidon_Units_Module, &
                    ONLY :  C_Square,                   &
                            Centimeter,                 &
                            Meter,                      &
                            GravPot_Units,              &
                            Shift_Units



USE Variables_Functions, &
                    ONLY :  Potential_Solution,         &
                            Calc_3D_Values_At_Location

USE Variables_IO, &
                    ONLY :  Report_IDs,        &
                            Iter_Report_Num_Samples,    &
                            Frame_Update_Table,         &
                            Iter_Time_Table,            &
                            Total_Run_Iters,            &
                            Iteration_Histogram,        &
                            iRF_Frame,                  &
                            iRF_Iter

USE Variables_Mesh, &
                    ONLY :  Num_R_Elements,             &
                            R_Inner,                    &
                            R_Outer

USE Poseidon_Parameters, &
                    ONLY :  DOMAIN_DIM,                 &
                            DEGREE,                     &
                            CUR_ITERATION,              &
                            MAX_ITERATIONS,             &
                            CONVERGENCE_CRITERIA,       &
                            CONVERGENCE_FLAG,           &
                            Verbose_Flag

USE Variables_NR,   &
                    ONLY :  NR_Coeff_Vector,         &
                            Block_RHS_Vector,           &
                            NR_Update_Vector,              &
                            BLOCK_ELEM_STF_MATVEC

USE Variables_MPI, &
                    ONLY :  myID_Poseidon,              &
                            myID_PETSC,                 &
                            myID_Shell,                 &
                            Local_Length,               &
                            Num_SubShells,              &
                            Num_R_Elems_Per_SubShell,   &
                            POSEIDON_COMM_PETSC,        &
                            POSEIDON_COMM_SHELL
                            

USE Variables_Derived, &
                    ONLY :  PROB_DIM,                   &
                            SUBSHELL_PROB_DIM,          &
                            ULM_LENGTH,                 &
                            LM_LENGTH
                            


USE CFA_1D_Master_Build_Module, &
                    ONLY :  CFA_1D_Master_Build

USE CFA_3D_Master_Build_Module, &
                    ONLY :  CFA_3D_Master_Build

USE IO_Write_Final_Results, &
            ONLY :  Write_Final_Results


USE IO_Linear_System, &
                    ONLY :  OUTPUT_JACOBIAN_MATRIX,             &
                            OUTPUT_RHS_VECTOR,                  &
                            OUTPUT_UPDATE_VECTOR,               &
                            READ_COEFFICIENT_VECTOR,            &
                            OUTPUT_COEFFICIENT_VECTOR_MATLAB,   &
                            OUTPUT_COEFFICIENT_VECTOR_FORTRAN

USE IO_Print_Results, &
                    ONLY :  Print_Results

USE Poseidon_Petsc_Solver, &
                    ONLY : PETSC_Distributed_Solve


USE Functions_Mapping, &
                    ONLY :  CFA_ALL_Matrix_Map


USE MPI

IMPLICIT NONE



CONTAINS







!+101+###########################################################################!
!                                                                                !
!                  CFA_Newton_Raphson                                            !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_Newton_Raphson_3D()


LOGICAL                                         :: CONVERGED

REAL(KIND = idp)                                :: timea, timeb, timec

LOGICAL                                                 :: PR = .FALSE.

timea = 0.0_idp
timeb = 0.0_idp
timec = 0.0_idp



IF ( lPF_IO_Flags(iPF_IO_Print_Results) ) THEN
    PR = .TRUE.
    WRITE(*,'(A)')"Initial Guess Values"
    CALL Print_Results()
    PRINT*," "
END IF


Cur_Iteration = 1
CONVERGED = .FALSE.
DO WHILE ( CONVERGED .EQV. .FALSE. )
!    PRINT*,"Begining of Iter ",Cur_Iteration


    IF (myID_Poseidon == 0 ) THEN
!        CALL OPEN_ITER_REPORT_FILE(Cur_Iteration, myID_Poseidon)
    END IF

    timea = MPI_Wtime()

    !*!
    !*! Clean the System
    !*!
    Block_RHS_Vector = 0.0_idp
    BLOCK_ELEM_STF_MATVEC = 0.0_idp


    !*!
    !*! Create the System
    !*!
    IF ( DOMAIN_DIM == 1 ) THEN
        PRINT*,"Before CFA_1D_Master_Build"
        CALL CFA_1D_Master_Build()
        PRINT*,"After CFA_1D_Master_Build"
    ELSEIF ( DOMAIN_DIM .NE. 1 ) THEN
        CALL CFA_3D_Master_Build()
    END IF



    !*!
    !*! Solve the System
    !*!
    CALL CFA_Solve()
    

    !*!
    !*! Update Coefficient_Vector
    !*!
    CALL CFA_Update_Share

    !*!
    !*!  Share Coefficient Vector
    !*!
    CALL CFA_Coefficient_Update_All


    !*!
    !*! Check for convergence
    !*!
    CALL CONVERGENCE_CHECK(CONVERGED, Cur_Iteration)



    !*!
    !*! Output Iteration Report
    !*!
    IF ( myID_Poseidon == 0 ) THEN
!        CALL OUTPUT_ITERATION_REPORT(Cur_Iteration, myID_Poseidon)
!        CALL CLOSE_ITER_REPORT_FILE()
        !CALL OUTPUT_COEFFICIENT_VECTOR_FORTRAN()
    END IF

    IF ( PR ) THEN
        WRITE(*,'(A,I3.3,A)')'Iteration ',Cur_Iteration,' Results'
        CALL Print_Results()
        PRINT*," "
    END IF





    IF ( Verbose_Flag .EQV. .TRUE. ) THEN
        WRITE(*,'(A,1X,I3.3,/)') "End of Iteration",Cur_Iteration
    END IF

   



    Cur_Iteration = Cur_Iteration + 1
    Total_Run_Iters = Total_Run_Iters + 1

END DO
Iteration_Histogram(Cur_Iteration-1) = Iteration_Histogram(Cur_Iteration-1) + 1







CALL Write_Final_Results()






END SUBROUTINE CFA_Newton_Raphson_3D













!+201+###########################################################################!
!                                                                                !
!                  CFA_Solve                                                     !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_Solve()


IF ( DOMAIN_DIM == 1 ) THEN


ELSE IF ( ( DOMAIN_DIM == 2 ) .OR. ( DOMAIN_DIM == 3 ) ) THEN


    CALL PETSC_Distributed_Solve( Block_Elem_STF_MatVec, Block_RHS_Vector, NR_Update_Vector)


END IF


END SUBROUTINE CFA_Solve










!+301+##########################################################################!
!                                                                               !
!               CFA_Coefficient_Update_All                                      !
!                                                                               !
!###############################################################################!
SUBROUTINE CFA_Coefficient_Update_All( )


!CALL CFA_Update_Modifier()

NR_Coeff_Vector = NR_Coeff_Vector + NR_Update_Vector

!PRINT*,"Update Vector is being output in CFA_Coefficient_Update_All,z.CFA_Newton_Raphson_3D.F90"
!DO i = 0,PROB_DIM-1
!    IF ( Coefficient_Vector(i) .NE. 0.0_idp ) THEN
!        PRINT*, NR_Update_Vector(i), NR_Update_Vector(i)/Coefficient_Vector(i)
!    ELSE
!        PRINT*,NR_Update_Vector(i)
!    END IF
!END DO


END SUBROUTINE CFA_Coefficient_Update_All


!+301a+##########################################################################!
!                                                                               !
!               CFA_Update_Modifier                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE CFA_Update_Modifier( )

INTEGER                                             :: re, d, lm_loc, ui, Here

IF ( DOMAIN_DIM .NE. 1 ) THEN
    PRINT*,"The Update Vector is Modified, in CFA_Update_Modifier, z.CFA_Newton_Raphon.F90"
    DO re = 0,NUM_R_ELEMENTS - 1
        DO d = 0,DEGREE
            DO ui = 4,5
                DO lm_loc = 0,LM_LENGTH-1

                    Here = CFA_ALL_Matrix_Map(ui, lm_loc, re, d)

                    NR_Update_Vector(Here) = 0.0_idp

                END DO
            END DO
        END DO
    END DO
END IF

END SUBROUTINE CFA_Update_Modifier



!+302+###########################################################################!
!                                                                                !
!               CFA_Coefficient_Update                                    !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_Coefficient_Update( )

INTEGER                                                 ::  Start_Here
INTEGER                                                 ::  End_Here


IF ( POSEIDON_COMM_PETSC .NE. MPI_COMM_NULL ) THEN

    Start_Here = myID_PETSc*NUM_R_ELEMS_PER_SUBSHELL*DEGREE*ULM_LENGTH
    End_Here = Start_Here + Local_Length - 1

    NR_Coeff_Vector(Start_Here:End_Here) = NR_Coeff_Vector(Start_Here:End_Here)   &
                                            + NR_Update_Vector(Start_Here:End_Here)

END IF

END SUBROUTINE CFA_Coefficient_Update







!+401+##########################################################################!
!                                                                               !
!               CFA_Coefficient_Share                                           !
!                                                                               !
!###############################################################################!
SUBROUTINE CFA_Coefficient_Share( )


INTEGER                                                 ::  i, ierr
INTEGER                                                 ::  Start_Here,     &
                                                            End_Here

INTEGER                                                 ::  Common_Length

INTEGER, DIMENSION(0:NUM_SUBSHELLS-1)                   ::  recvcounts
INTEGER, DIMENSION(0:NUM_SUBSHELLS-1)                   ::  displs
COMPLEX(kind = idp), DIMENSION(0:Local_Length-1)        ::  Send_Buffer_Block
COMPLEX(kind = idp), DIMENSION(0:PROB_DIM-1)            ::  Recieve_Buffer


!
!  Gather Parts of Coefficient_Vector from PETSc Processes
!
IF ( POSEIDON_COMM_PETSC .NE. MPI_COMM_NULL ) THEN


    Common_Length = SUBSHELL_PROB_DIM - ULM_LENGTH



    !
    !   Create recieve counts array
    !
    recvcounts(:) = Common_Length
    recvcounts(NUM_SUBSHELLS-1) = SUBSHELL_PROB_DIM


    !
    ! Create displacement array
    !
    displs(0) = 0
    DO i = 1,NUM_SUBSHELLS-1
        displs(i) = displs(i-1)+Common_Length
    END DO



    !
    !   Load Send Buffer
    !
    Start_Here = myID_PETSc*NUM_R_ELEMS_PER_SUBSHELL*DEGREE*ULM_LENGTH
    End_Here = Start_Here + Local_Length - 1

    Send_Buffer_Block = NR_Coeff_Vector(Start_Here:End_Here)



    !
    !   Share Data
    !
    CALL MPI_Allgatherv(Send_Buffer_Block,   &
                        Local_Length,        &
                        MPI_DOUBLE_COMPLEX,  &
                        Recieve_Buffer,      &
                        recvcounts,          &
                        displs,              &
                        MPI_DOUBLE_COMPLEX,  &
                        POSEIDON_COMM_PETSC, &
                        ierr                 )




    NR_Coeff_Vector = Recieve_Buffer

END IF  ! Poseidon_Comm_PETSC .NE. MPI_COMM_NULL








    !
    !   Share Coefficient Vector with blocks in Shell
    !
IF ( myID_Shell == 0 ) THEN

    CALL MPI_BCAST( NR_Coeff_Vector,     &
                    PROB_DIM,               &
                    MPI_DOUBLE_COMPLEX,     &
                    0,                      &
                    POSEIDON_COMM_SHELL,    &
                    ierr                    )

ELSE

    CALL MPI_BCAST( Recieve_Buffer,         &
                    PROB_DIM,               &
                    MPI_DOUBLE_COMPLEX,     &
                    0,                      &
                    POSEIDON_COMM_SHELL,    &
                    ierr                    )
   
   NR_Coeff_Vector = Recieve_Buffer

END IF




END SUBROUTINE CFA_Coefficient_Share













!+402+###########################################################################!
!                                                                                !
!               CFA_Update_Share                                    !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_Update_Share( )


INTEGER                                                 ::  i, ierr
INTEGER                                                 ::  Start_Here,     &
                                                            End_Here

INTEGER                                                 ::  Common_Length

INTEGER, DIMENSION(0:NUM_SUBSHELLS-1)                   ::  recvcounts
INTEGER, DIMENSION(0:NUM_SUBSHELLS-1)                   ::  displs
COMPLEX(kind = idp), DIMENSION(0:Local_Length-1)        ::  Send_Buffer_Block
COMPLEX(kind = idp), DIMENSION(0:PROB_DIM-1)            ::  Recieve_Buffer


!
!  Gather Parts of Coefficient_Vector from PETSc Processes
!
IF ( POSEIDON_COMM_PETSC .NE. MPI_COMM_NULL ) THEN


    Common_Length = SUBSHELL_PROB_DIM - ULM_LENGTH



    !
    !   Create recieve counts array
    !
    recvcounts(:) = Common_Length
    recvcounts(NUM_SUBSHELLS-1) = SUBSHELL_PROB_DIM


    !
    ! Create displacement array
    !
    displs(0) = 0
    DO i = 1,NUM_SUBSHELLS-1
        displs(i) = displs(i-1)+Common_Length
    END DO



    !
    !   Load Send Buffer
    !
    Start_Here = myID_PETSc*NUM_R_ELEMS_PER_SUBSHELL*DEGREE*ULM_LENGTH
    End_Here = Start_Here + Local_Length - 1

    Send_Buffer_Block = NR_Update_Vector(Start_Here:End_Here)



    !
    !   Share Data
    !
    CALL MPI_Allgatherv(Send_Buffer_Block,   &
                        Local_Length,        &
                        MPI_DOUBLE_COMPLEX,  &
                        Recieve_Buffer,      &
                        recvcounts,          &
                        displs,              &
                        MPI_DOUBLE_COMPLEX,  &
                        POSEIDON_COMM_PETSC, &
                        ierr                 )




    NR_Update_Vector = Recieve_Buffer

END IF  ! Poseidon_Comm_PETSC .NE. MPI_COMM_NULL








    !
    !   Share Coefficient Vector with blocks in Shell
    !
IF ( myID_Shell == 0 ) THEN

    CALL MPI_BCAST( NR_Update_Vector,     &
                    PROB_DIM,               &
                    MPI_DOUBLE_COMPLEX,     &
                    0,                      &
                    POSEIDON_COMM_SHELL,    &
                    ierr                    )

ELSE

    CALL MPI_BCAST( Recieve_Buffer,         &
                    PROB_DIM,               &
                    MPI_DOUBLE_COMPLEX,     &
                    0,                      &
                    POSEIDON_COMM_SHELL,    &
                    ierr                    )

   NR_Update_Vector = Recieve_Buffer

END IF


END SUBROUTINE CFA_Update_Share

















!+501+##########################################################################!
!                                                                               !
!                   OUTPUT_ITERATION_REPORT                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE OUTPUT_ITERATION_REPORT(Iter, Rank)

INTEGER, INTENT(IN)                 :: Iter, Rank

INTEGER, DIMENSION(0:1)                         ::  FILE_ID
INTEGER                                         ::  i, j
REAL(KIND = idp)                                ::  r, theta, phi, deltar
REAL(KIND = idp)                                ::  Analytic_Val
REAL(KIND = idp)                                ::  Return_Psi, Return_AlphaPsi
REAL(KIND = idp)                                ::  Return_Beta1, Return_Beta2, Return_Beta3
REAL(KIND = idp)                                ::  PsiPot_Val, AlphaPsiPot_Val


120 FORMAT (A61)
121 FORMAT (A1)
123 FORMAT (A38,ES22.15)

110 FORMAT (11X,A1,16X,A18,9X,A13,10X,A18,10X,A11,14X,A11,14X,A11)
111 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)





FILE_ID = -1
IF ( lPF_IO_Flags(iPF_IO_Write_Frame_Report) ) THEN

    FILE_ID(0) = Report_IDs(iRF_Frame)

END IF

IF ( lPF_IO_Flags(iPF_IO_Write_Iter_Report) ) THEN
    FILE_ID(1) = Report_IDs(iRF_Iter)
END IF



! Write Title to File
DO i = 0,1

    IF ( FILE_ID(i) .NE. -1 ) THEN

        WRITE(FILE_ID(i),'(A)')"                                     Timing Results"
        WRITE(FILE_ID(i),'(A)')"            ============================================================="
        WRITE(FILE_ID(i),'(A)')" "
        WRITE(FILE_ID(i),123)"                    Initialize Time : ",ITER_TIME_TABLE(1)
        WRITE(FILE_ID(i),123)" Input/Communicate Source Data Time : ",ITER_TIME_TABLE(2)
        WRITE(FILE_ID(i),123)"     Input Boundary Conditions Time : ",ITER_TIME_TABLE(3)
        WRITE(FILE_ID(i),123)"        CFA_3D_Apply_BCs_Part1 Time : ",ITER_TIME_TABLE(4)
        WRITE(FILE_ID(i),120)"-------------------------------------------------------------"
        WRITE(FILE_ID(i),123)" ||     Calc_3D_Current_Values Time : ",ITER_TIME_TABLE(5)
        WRITE(FILE_ID(i),123)" ||    CREATE_3D_SubJcbn_Terms Time : ",ITER_TIME_TABLE(6)
        WRITE(FILE_ID(i),123)" ||       CREATE_3D_RHS_VECTOR Time : ",ITER_TIME_TABLE(7)
        WRITE(FILE_ID(i),123)"\  /     CREATE_3D_JCBN_MATRIX Time : ",ITER_TIME_TABLE(8)
        WRITE(FILE_ID(i),120)"-\/ ---------------------------------------------------------"
        WRITE(FILE_ID(i),123)"CREATE_3D_NONLAPLACIAN_STF_MAT Time : ",ITER_TIME_TABLE(9)
        WRITE(FILE_ID(i),123)"REDUCE_3D_NONLAPLACIAN_STF_MAT Time : ",ITER_TIME_TABLE(10)
        WRITE(FILE_ID(i),123)"FINISH_3D_NONLAPLACIAN_STF_MAT Time : ",ITER_TIME_TABLE(11)
        WRITE(FILE_ID(i),123)"          FINISH_3D_RHS_VECTOR Time : ",ITER_TIME_TABLE(12)
        WRITE(FILE_ID(i),123)"        CFA_3D_Apply_BCs_Part2 Time : ",ITER_TIME_TABLE(13)
        WRITE(FILE_ID(i),123)"                    CFA_Solver Time : ",ITER_TIME_TABLE(14)
        WRITE(FILE_ID(i),123)"        CFA_Coefficient_Update Time : ",ITER_TIME_TABLE(15)
        WRITE(FILE_ID(i),123)"   CFA_Coefficient_Share_PETSc Time : ",ITER_TIME_TABLE(16)
        WRITE(FILE_ID(i),123)"         CFA_Convergence_Check Time : ",ITER_TIME_TABLE(17)
        WRITE(FILE_ID(i),123)"               Total Iteration Time : ",ITER_TIME_TABLE(18)
        WRITE(FILE_ID(i),123)"             Poseidon_Dist_Sol Time : ",ITER_TIME_TABLE(19)
        WRITE(FILE_ID(i),120)"============================================================="
        WRITE(FILE_ID(i),121)" "
        WRITE(FILE_ID(i),121)" "
        WRITE(FILE_ID(i),121)" "
        WRITE(FILE_ID(i),121)" "

        WRITE(FILE_ID(i),'(A,I2.2,A)')"                        Iteration ",Iter," Results"
        WRITE(FILE_ID(i),'(A)')"            ============================================================="
        WRITE(FILE_ID(i),'(A)')" "
        WRITE(FILE_ID(i),110)"r","Analytic Potential","Psi Potential","AlphaPsi Potential","Beta Value1","Beta Value2","Beta Value3"

    END IF
END DO






IF ( ASSOCIATED(Potential_Solution) ) THEN

    ! Write Results Table Header to Screen
    IF ( lPF_IO_Flags(iPF_IO_Print_Frame_Report) ) THEN
        IF ( Rank == 0 ) THEN
            PRINT*,"+++++++++++++++++++ myID,",Rank," Iteration",Iter,"++++++++++++++++++++++"
            WRITE(*,110)"r","Analytic Potential","Psi Potential","AlphaPsi Potential","Beta Value1","Beta Value2","Beta Value3"
        END IF
    END IF

    deltar = ( R_OUTER - R_INNER )/ REAL(ITER_REPORT_NUM_SAMPLES, KIND = idp)
    DO i = 0,ITER_REPORT_NUM_SAMPLES

        r = i*deltar + R_INNER
        theta = pi/2.0_idp
        phi = pi/2.0_idp


        CALL Calc_3D_Values_At_Location( r, theta, phi,                              &
                                        Return_Psi, Return_AlphaPsi,                &
                                        Return_Beta1, Return_Beta2, Return_Beta3    )

        ! Determine the Newtonian Potential at the location r, theta, phi
        Analytic_Val = Potential_Solution(r,theta,phi)/GravPot_Units


        ! AlphaPsi_to_Pot   =   2*C_Square*(AlphaPsi - 1)
        ! Psi_to_Pot        =   2*C_Square*(1 - Psi)

        ! Calculate Conformal Factor value from Newtonian Potential
        PsiPot_Val = 2.0_idp*C_Square*(1.0_idp - Return_Psi)/GravPot_Units

        ! Calculate the product of the Conformal Factor and Lapse Function from Newtonian Potential
        AlphaPsiPot_Val = 2.0_idp*C_Square*(Return_AlphaPsi - 1.0_idp)/GravPot_Units


        ! Write Results to Screen
        IF ( lPF_IO_Flags(iPF_IO_Print_Frame_Report) ) THEN
            IF ( Rank == 0 ) THEN
                WRITE(*,111) r/Centimeter,              &
                             Analytic_Val,              &
                             PsiPot_Val,                &
                             AlphaPsiPot_Val,           &
                             Return_Beta1/Shift_Units,  &
                             Return_Beta2,              &
                             Return_Beta3
            END IF
        END IF

        ! Write Results to File
        DO j = 0,1
            IF ( FILE_ID(j) .NE. -1 ) THEN
                WRITE(FILE_ID(j),111) r/Centimeter,              &
                                      Analytic_Val,              &
                                      PsiPot_Val,                &
                                      AlphaPsiPot_Val,           &
                                      Return_Beta1/Shift_Units,  &
                                      Return_Beta2,              &
                                      Return_Beta3
            END IF
        END DO ! j Loop

    END DO  ! i Loop

ELSE

    ! Write Results Table Header to Screen
    IF ( lPF_IO_Flags(iPF_IO_Print_Frame_Report) ) THEN
        IF ( Rank == 0 ) THEN
            PRINT*,"+++++++++++++++++++ myID,",Rank," Iteration",Iter,"++++++++++++++++++++++"
            WRITE(*,110)"r","Psi Potential","AlphaPsi Potential","Beta Value1","Beta Value2","Beta Value3"
        END IF
    END IF

    deltar = ( R_OUTER - R_INNER )/ REAL(ITER_REPORT_NUM_SAMPLES, KIND = idp)
    DO i = 0,ITER_REPORT_NUM_SAMPLES

        r = i*deltar + R_INNER
        theta = pi/2.0_idp
        phi = pi/2.0_idp


        CALL Calc_3D_Values_At_Location( r, theta, phi,                              &
                                        Return_Psi, Return_AlphaPsi,                &
                                        Return_Beta1, Return_Beta2, Return_Beta3    )

        ! AlphaPsi_to_Pot   =   2*C_Square*(AlphaPsi - 1)
        ! Psi_to_Pot        =   2*C_Square*(1 - Psi)

        ! Calculate Conformal Factor value from Newtonian Potential
        PsiPot_Val = 2.0_idp*C_Square*(1.0_idp - Return_Psi)/GravPot_Units

        ! Calculate the product of the Conformal Factor and Lapse Function from Newtonian Potential
        AlphaPsiPot_Val = 2.0_idp*C_Square*(Return_AlphaPsi - 1.0_idp)/GravPot_Units


        ! Write Results to Screen
        IF ( lPF_IO_Flags(iPF_IO_Print_Frame_Report) ) THEN
            IF ( Rank == 0 ) THEN
                WRITE(*,111) r/Centimeter,              &
                             PsiPot_Val,                &
                             AlphaPsiPot_Val,           &
                             Return_Beta1/Shift_Units,  &
                             Return_Beta2,              &
                             Return_Beta3
            END IF
        END IF

        ! Write Results to File
        DO j = 0,1
            IF ( FILE_ID(j) .NE. -1 ) THEN
                WRITE(FILE_ID(j),111) r/Centimeter,              &
                                      PsiPot_Val,                &
                                      AlphaPsiPot_Val,           &
                                      Return_Beta1/Shift_Units,  &
                                      Return_Beta2,              &
                                      Return_Beta3
            END IF
        END DO ! j Loop

    END DO  ! i Loop




END IF



WRITE( Report_IDs(iRF_Iter), '(4/)')

END SUBROUTINE OUTPUT_ITERATION_REPORT









!+601+##########################################################################!
!                                                                               !
!                       CONVERGENCE_CHECK                                       !
!                                                                               !
!###############################################################################!
SUBROUTINE CONVERGENCE_CHECK(CONVERGED, Iteration)

LOGICAL, INTENT(INOUT)              :: CONVERGED
INTEGER, INTENT(IN)                 :: Iteration

INTEGER, DIMENSION(0:1)                                 ::  FILE_ID
INTEGER                                                 ::  i

COMPLEX(KIND = idp)                                     ::  RMS_VALUE
COMPLEX(KIND = idp)                                     ::  Euclidean

REAL(KIND = idp)                                        ::  MaxA
INTEGER                                                 ::  MaxA_Loc
REAL(KIND = idp)                                        ::  Ratio_Real, Ratio_Imag

FILE_ID = -1
IF ( lPF_IO_Flags(iPF_IO_Write_Frame_Report) ) THEN
    FILE_ID(0) = Report_IDs(iRF_Frame)
END IF

IF ( lPF_IO_Flags(iPF_IO_Write_Iter_Report) ) THEN
    FILE_ID(1) = Report_IDs(iRF_Iter)
END IF


MaxA = 0.0_idp
!PRINT*,"In CONVERGENCE_CHECK"
DO i = 0,PROB_DIM-1

    if ( NR_Coeff_Vector(i) .NE. 0.0_idp ) THEN
    Ratio_Real = REAL( NR_Update_Vector(i)/NR_Coeff_Vector(i) , KIND = idp)
    Ratio_Imag = REAL(IMAG(  NR_Update_Vector(i)/NR_Coeff_Vector(i) ), KIND = idp)


        IF (  abs(Ratio_Real) > MaxA  .OR. abs(Ratio_Imag) > MaxA   ) THEN
            MaxA = abs(NR_Update_Vector(i)/NR_Coeff_Vector(i))
            MaxA_Loc = i
        END IF
    END IF
END DO

!PRINT*,"Max Ratio",MaxA, MaxA_Loc





IF (myID_Poseidon == 0) THEN

    RMS_VALUE = SQRT( SUM(NR_Update_Vector**2)/size(NR_Update_Vector) )
    Euclidean = SQRT( SUM(ABS(NR_Update_Vector)**2) )


    IF ( MAXVAL(ABS(NR_Update_Vector)) .LE. Convergence_Criteria ) THEN

        CONVERGED = .TRUE.
        CONVERGENCE_FLAG = 1

    ELSE IF ( Iteration .GE. Max_Iterations ) THEN

        CONVERGED = .TRUE.
        CONVERGENCE_FLAG = 2

    END IF


!    PRINT*,"RMS ",RMS_VALUE," Euclidean ",Euclidean," MAXVAL ",MAXVAL(ABS(NR_Update_Vector))

!   Add the current iteration's convergence check value to the frame table
    Frame_Update_Table(Iteration,1) =    MAXVAL(ABS(NR_Update_Vector))

!   Output data to files
    DO i = 0,1
        IF ( FILE_ID(i) .NE. -1 ) THEN
            WRITE(FILE_ID(i),'(A)')"                                 Convergence Results"
            WRITE(FILE_ID(i),'(A)')"            ============================================================="
            WRITE(FILE_ID(i),'(A)')""
            WRITE(FILE_ID(i),'(A,ES22.15,A,ES22.15)')"Convergence Check: MAXVAL(ABS(NR_Update_Vector)) ",      &
                                                     MAXVAL(ABS(NR_Update_Vector)),                            &
                                                     "  Criteria ",Convergence_Criteria
            IF ( CONVERGENCE_FLAG == 1 ) THEN
                WRITE(FILE_ID(i),'(A,I2.2,A,ES22.15,A,ES22.15)')"CONVERGED IN ", Iteration," Iterations"
            ELSE IF ( CONVERGENCE_FLAG == 2 ) THEN
                WRITE(FILE_ID(i),'(A)')"Iteration Exceeded Max Iteration Limit"
            END IF
            WRITE(FILE_ID(i),'(3/)')



        END IF
    END DO


!   Output data to screen
    IF ( lPF_IO_Flags(iPF_IO_Print_Frame_Report) ) THEN

        WRITE(*,'(A,ES22.15,A,ES22.15)')"Convergence Check: MAXVAL(ABS(NR_Update_Vector)) ",MAXVAL(ABS(NR_Update_Vector)),   &
                                              "  Criteria ",Convergence_Criteria
    END IF



END IF






END SUBROUTINE CONVERGENCE_CHECK






END MODULE CFA_Newton_Raphson_3D_Module
