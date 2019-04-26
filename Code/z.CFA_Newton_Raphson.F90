   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE CFA_Newton_Raphson_Module                                                    !##!
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
!##!    +301+   CFA_Coefficient_Update                                              !##!
!##!    +302+   CFA_Coefficient_Share                                               !##!
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
                        ONLY : idp, pi, speed_of_light

USE CHIMERA_Parameters, &
                        ONLY :  Analytic_Solution

USE Poseidon_Parameters, &
                        ONLY :  DOMAIN_DIM,                 &
                                DEGREE,                     &
                                L_LIMIT,                    &
                                NUM_SHELLS,                 &
                                NUM_SUBSHELLS,              &
                                nPROCS_POSEIDON,            &
                                STF_MAPPING_FLAG,           &
                                NUM_R_ELEMS_PER_BLOCK,      &
                                NUM_R_ELEMS_PER_SUBSHELL,   &
                                MAX_ITERATIONS,             &
                                CONVERGENCE_CRITERIA

USE Global_Variables_And_Parameters, &
                        ONLY :  NUM_R_ELEMENTS,             &
                                NUM_T_ELEMENTS,             &
                                NUM_P_ELEMENTS,             &
                                rlocs, tlocs, plocs,        &
                                NUM_R_NODES,                &
                                INT_R_LOCATIONS,            &
                                INT_T_LOCATIONS,            &
                                INT_P_LOCATIONS,            &
                                INT_R_WEIGHTS,              &
                                INT_T_WEIGHTS,              &
                                INT_P_WEIGHTS,              &
                                VAR_DIM,                    &
                                Lagrange_Poly_Table,        &
                                Coefficient_Vector,         &
                                RHS_Vector,                 &
                                Block_RHS_Vector,           &
                                myID_Poseidon,              &
                                myID_PETSC,                 &
                                myID_Shell,                 &
                                myShell,                    &
                                BLOCK_ELEM_STF_MATVEC,      &
                                POSEIDON_COMM_WORLD,        &
                                POSEIDON_COMM_PETSC,        &
                                POSEIDON_COMM_SHELL,        &
                                PROB_DIM,                   &
                                Block_Prob_Dim,             &
                                SUBSHELL_PROB_DIM,          &
                                Local_Length,               &
                                ULM_LENGTH,                 &
                                LM_LENGTH,                  &
                                Update_Vector,              &
                                NUM_OFF_DIAGONALS,          &
                                R_INNER, R_OUTER



USE CFA_3D_Master_Build_Module, &
                        ONLY :  CFA_3D_Master_Build,      &
                                CFA_3D_Apply_BCs_Part1,   &
                                CFA_3D_Apply_BCs_Part2,   &
                                Calc_3D_Values_At_Location


USE Jacobian_Internal_Functions_Module,  &
                        ONLY :  Initialize_Guess_Values,        &
                                Initialize_Special_Guess_Values

USE IO_Functions_Module, &
                        ONLY : Clock_In

USE Poseidon_Petsc_Solver, &
                        ONLY : PETSC_Distributed_Solve


USE Additional_Functions_Module, &
                                ONLY :  Lagrange_Poly,              &
                                        CFA_ALL_Matrix_Map

USE MPI

IMPLICIT NONE



CONTAINS







!+101+###########################################################################!
!                                                                                !
!                  CFA_Newton_Raphson                                            !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_Newton_Raphson()


LOGICAL                                         :: CONVERGED = .FALSE.

INTEGER                                         :: Iteration

INTEGER, DIMENSION(1:2+DOMAIN_DIM)              :: VARIABLE_RESTRICT

INTEGER                                         :: i
INTEGER                                         :: ierr


INTEGER                                         :: j, NUM_SAMPLES
REAL(KIND = idp)                                :: r, theta,phi, deltar,     &
                                                   Analytic_Val, Solver_Val,  Error_Val
REAL(KIND = idp)                                :: Return_Psi, Return_AlphaPsi, Return_Beta1,  &
                                                   Return_Beta2, Return_Beta3

REAL(KIND = idp)                                :: PsiPot_Val, AlphaPsiPot_Val

REAL(KIND = idp)       :: csqr
REAL(KIND = idp)       :: timea, timeb, timec



110 FORMAT (11X,A1,16X,A18,9X,A13,10X,A18,10X,A10)
111 FORMAT (11X,A1,16X,A18,9X,A13,10X,A18,10X,A11,10X,A11)
113 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)
114 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)


timea = 0.0_idp
timeb = 0.0_idp
timec = 0.0_idp


CALL Initialize_Special_Guess_Values()
!CALL Initialize_Guess_Values()


csqr = speed_of_light*speed_of_light




NUM_SAMPLES = 20
!IF ( .TRUE. ) THEN
!IF ( .FALSE. ) THEN
IF ( myID_Poseidon == 0 ) THEN
    PRINT*,"+++++++++++++++++++ myID,",myID_Poseidon," Initial Guess ++++++++++++++++++++++"
    deltar = ( R_OUTER - R_INNER )/ REAL(NUM_SAMPLES, KIND = idp)


!    PRINT*,"    r            Analytic       Solver        Error"
    WRITE(*,110)"r","Analytic Potential","Psi Potential","AlphaPsi Potential","Beta Value"

    DO i = 0,NUM_SAMPLES

        r = i*deltar + R_INNER
        theta = pi/2.0_idp
        phi = pi/2.0_idp
           





        CALL Calc_3D_Values_At_Location( r, theta, phi,                              &
                                         Return_Psi, Return_AlphaPsi,                &
                                         Return_Beta1, Return_Beta2, Return_Beta3    )

        Analytic_Val = Analytic_Solution(r,theta,phi)
        Solver_Val = 2.0_idp*csqr*(1.0_idp - Return_Psi)
        Error_Val = ABS((Analytic_Val - Solver_Val)/Analytic_Val)

        PsiPot_Val = 2.0_idp*csqr*(1.0_idp - Return_Psi)
        AlphaPsiPot_Val = 2.0_idp*csqr*(Return_AlphaPsi - 1.0_idp)

    ! AlphaPsi_to_Pot = 2*(AlphaPsi - 1)
    ! Psi_to_Pot = 2*(1 - Psi) or 0.5(1 - Psi^4)

!       print*,r,Analytic_Val, Solver_Val,Error_Val
        WRITE(*,113) r,Analytic_Val, PsiPot_Val, AlphaPsiPot_Val, Return_Beta1

   END DO
END IF
!CALL MPI_BARRIER(POSEIDON_COMM_WORLD, ierr)
!CALL SLEEP(1)












Iteration = 1
DO WHILE ( CONVERGED .EQV. .FALSE. )

    timea = MPI_Wtime()

    Block_RHS_Vector = 0.0_idp
    BLOCK_ELEM_STF_MATVEC = 0.0_idp

    IF ( DOMAIN_DIM == 1 ) THEN




    ELSE IF ( DOMAIN_DIM == 2 ) THEN

        timeb = MPI_Wtime()
        CALL CFA_3D_Apply_BCs_Part1()
        timec = MPI_Wtime()
        CALL Clock_In(timec-timeb, 4)







        CALL CFA_3D_Master_Build()





        timeb = MPI_Wtime()
        CALL CFA_3D_Apply_BCs_Part2()
        timec = MPI_Wtime()
        CALL Clock_In(timec-timeb, 13)
                                          




    ELSE IF ( DOMAIN_DIM == 3 ) THEN
        

        timeb = MPI_Wtime()
        CALL CFA_3D_Apply_BCs_Part1()
        timec = MPI_Wtime()
        CALL Clock_In(timec-timeb, 4)






        CALL CFA_3D_Master_Build()




        timeb = MPI_Wtime()
        CALL CFA_3D_Apply_BCs_Part2()
        timec = MPI_Wtime()
        CALL Clock_In(timec-timeb, 13)




    END IF



    IF (myID_PETSC == -130) THEN
        PRINT*,"BEFORE SOLVE Matrix",myID_Poseidon
        DO i = 0,NUM_SUBSHELLS-1
            DO j = 0,NUM_R_ELEMS_PER_BLOCK-1
                PRINT*,"RE = ",j
                PRINT*,REAL(BLOCK_ELEM_STF_MATVEC(:,j),KIND = idp)
                PRINT*," "
            END DO
        END DO
    END IF
    PRINT*," "






    !*!
    !*! Solve the System
    !*!


    timeb = MPI_Wtime()
    CALL CFA_Solve()
    timec = MPI_Wtime()
    CALL Clock_In(timec-timeb, 14)






     !*!
     !*! Update Coefficient_Vector
     !*!


     timeb = MPI_Wtime()

     CALL CFA_Update_Share
     timec = MPI_Wtime()
     CALL Clock_In(timec-timeb, 16)




     !*!
     !*!  Share Coefficient Vector
     !*! 


     timeb = MPI_Wtime()
     CALL CFA_Coefficient_Update_All
     timec = MPI_Wtime()
     CALL Clock_In(timec-timeb, 15)






    !*!
    !*! Check for convergence
    !*!

    timeb = MPI_Wtime()

    IF (myID_Poseidon == 0) THEN

        PRINT*,"Convergence Check: MAXVAL(ABS(Update_Vector))",MAXVAL(ABS(Update_Vector)),"Criteria ",Convergence_Criteria 

    END IF
    IF ( MAXVAL(ABS(Update_Vector)) .LE. Convergence_Criteria ) THEN
        IF (myID_Poseidon == 0) THEN
            PRINT*,"CONVERGED IN",Iteration,"Iterations",MAXVAL(ABS(Update_Vector)),"<=",Convergence_Criteria
            PRINT*," "
        END IF
        CONVERGED = .TRUE.

    ELSE IF ( Iteration .GE. Max_Iterations ) THEN
        IF (myID_Poseidon == 0) THEN
            PRINT*,"Iteration Exceeded Max Iteration Limit"
            PRINT*," "
        END IF
        CONVERGED = .TRUE.

    END IF





    timec = MPI_Wtime()
    IF (myID_Poseidon == 0) THEN
        Print*,"         CFA_Convergence_Check :",timec-timeb
    END IF
    CALL Clock_In(timec-timeb, 17)


    Iteration = Iteration + 1

    timeb = MPI_Wtime()
!    Print*,"          Total Iteration Time :",timeb-timea

    CALL Clock_In(timeb-timea,18)


    CALL SLEEP(1)
    DO j = 0,0 !nPROCS_POSEIDON-1
      IF ( myID_Poseidon == 0 ) THEN
        PRINT*,"+++++++++++++++++++ myID,",myID_Poseidon," Iteration",Iteration-1,"++++++++++++++++++++++"
        deltar = ( R_OUTER - R_INNER )/ REAL(NUM_SAMPLES, KIND = idp)


!        WRITE(*,110)"r","Analytic Potential","Psi Potential","AlphaPsi Potential","Beta Value"
        WRITE(*,111)"r","Analytic Potential","Psi Potential","AlphaPsi Potential","Beta Value1","Beta Value2"

        DO i = 0,NUM_SAMPLES

           r = i*deltar + R_INNER
           theta = pi/2.0_idp
           phi = pi/2.0_idp






           CALL Calc_3D_Values_At_Location( r, theta, phi,                              &
                                            Return_Psi, Return_AlphaPsi,                &
                                            Return_Beta1, Return_Beta2, Return_Beta3    )


          Analytic_Val = Analytic_Solution(r,theta,phi)
          PsiPot_Val = 2.0_idp*csqr*(1.0_idp - Return_Psi)
          AlphaPsiPot_Val = 2.0_idp*csqr*(Return_AlphaPsi - 1.0_idp)
          Error_Val = ABS((Analytic_Val - Solver_Val)/Analytic_Val)

          ! AlphaPsi_to_Pot = 2*(AlphaPsi - 1)
          ! Psi_to_Pot = 2*(1 - Psi) or 0.5(1 - Psi^4)

!          WRITE(*,113) r,Analytic_Val, PsiPot_Val, AlphaPsiPot_Val, Return_Beta1
          WRITE(*,114) r,Analytic_Val, PsiPot_Val, AlphaPsiPot_Val, Return_Beta1,Return_Beta2

       END DO
      END IF
!     CALL MPI_BARRIER(POSEIDON_COMM_WORLD, ierr)
!     CALL SLEEP(1)
   END DO



END DO


END SUBROUTINE CFA_Newton_Raphson
















!+201+###########################################################################!
!                                                                                !
!                  CFA_Solve                                                     !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_Solve()


CHARACTER(LEN = 1)                                  :: TRANS

INTEGER                                             :: M, N, NRHS, LDA, LDB, INFO

INTEGER, DIMENSION(:), ALLOCATABLE                  :: IPIV





IF ( DOMAIN_DIM == 1 ) THEN




ELSE IF ( ( DOMAIN_DIM == 2 ) .OR. ( DOMAIN_DIM == 3 ) ) THEN


    CALL PETSC_Distributed_Solve()


END IF







END SUBROUTINE CFA_Solve





!+301+###########################################################################!
!                                                                                !
!               CFA_Coefficient_Update                                    !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_Coefficient_Update_All( )

INTEGER         :: re, d, lm_loc, CUR_I_LOC



!DO re = 0,NUM_R_ELEMS_PER_BLOCK -1
!  DO d = 0,DEGREE
!    DO lm_loc = 0,LM_LENGTH-1

!        Cur_i_Loc = CFA_ALL_Matrix_Map(1, lm_loc, re, d)

!        PRINT*,"***********   Update_Vector(cur_i_loc + 0:0) = 0.0    ************"

!        Update_Vector(Cur_i_Loc+0:Cur_i_Loc+0) = 0.0_idp

!    END DO  ! l Loop
!  END DO  ! d Loop
!END DO




Coefficient_Vector = Coefficient_Vector + Update_Vector






END SUBROUTINE CFA_Coefficient_Update_All






!+301+###########################################################################!
!                                                                                !
!               CFA_Coefficient_Update                                    !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_Coefficient_Update( )


INTEGER                                                 ::  i
INTEGER                                                 ::  Start_Here,     &
                                                            End_Here

INTEGER   :: ierr

!Start_Here = myID_PETSc*NUM_R_ELEMS_PER_BLOCK*DEGREE*ULM_LENGTH
!End_Here = Start_Here + Local_Length
!Print*,"Update",myID,"Orig Start_Here ",Start_Here," End_Here ",End_Here

!DO  i = 0,nPROCS-1
!    IF (myID == i) THEN

!        PRINT*,"COEFFS_BEFORE update"
!        PRINT*,Coefficient_Vector
!        PRINT*," "

!    END IF                 
!    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!END DO



IF ( POSEIDON_COMM_PETSC .NE. MPI_COMM_NULL ) THEN


    Start_Here = myID_PETSc*NUM_R_ELEMS_PER_SUBSHELL*DEGREE*ULM_LENGTH
    End_Here = Start_Here + Local_Length - 1
!Print*,"Update",myID,"New  Start_Here ",Start_Here," End_Here ",End_Here


    Coefficient_Vector(Start_Here:End_Here) = Coefficient_Vector(Start_Here:End_Here)   &
                                            + Update_Vector(Start_Here:End_Here)

END IF

!DO  i = 0,nPROCS-1
!    IF (myID == i) THEN
    
!        PRINT*,"COEFFS_AFTER update"
!        PRINT*,Coefficient_Vector
!        PRINT*," "

!    END IF                 
!    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!END DO



END SUBROUTINE CFA_Coefficient_Update







!+301+###########################################################################!
!                                                                                !
!               CFA_Coefficient_Share                                    !
!                                                                                !
!################################################################################!
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

    Send_Buffer_Block = Coefficient_Vector(Start_Here:End_Here)



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




    Coefficient_Vector = Recieve_Buffer

END IF  ! Poseidon_Comm_PETSC .NE. MPI_COMM_NULL








    !
    !   Share Coefficient Vector with blocks in Shell
    !
IF ( myID_Shell == 0 ) THEN

    CALL MPI_BCAST( Coefficient_Vector,     &
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
   
   Coefficient_Vector = Recieve_Buffer

END IF

!DO  i = 0,nPROCS-1
!    IF (myID == i) THEN
    
!        PRINT*,"COEFFS After Share"
!        PRINT*,Coefficient_Vector
!        PRINT*," "

!    END IF                 
!    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!END DO




END SUBROUTINE CFA_Coefficient_Share













!+301+###########################################################################!
!                                                                                !
!               CFA_Coefficient_Share                                    !
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

    Send_Buffer_Block = Update_Vector(Start_Here:End_Here)



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




    Update_Vector = Recieve_Buffer

END IF  ! Poseidon_Comm_PETSC .NE. MPI_COMM_NULL








    !
    !   Share Coefficient Vector with blocks in Shell
    !
IF ( myID_Shell == 0 ) THEN

    CALL MPI_BCAST( Update_Vector,     &
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

   Update_Vector = Recieve_Buffer

END IF

!DO  i = 0,nPROCS-1
!    IF (myID == i) THEN

!        PRINT*,"COEFFS After Share"
!        PRINT*,Coefficient_Vector
!        PRINT*," "

!    END IF                 
!    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!END DO




END SUBROUTINE CFA_Update_Share
























END MODULE CFA_Newton_Raphson_Module
