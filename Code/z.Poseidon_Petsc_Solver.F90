   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Petsc_Solver                                                        !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!================================================================================!##!
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


USE OMP_LIB


USE Poseidon_Constants_Module, &
                                ONLY :  idp

USE Poseidon_Parameters, &
                                ONLY :  DOMAIN_DIM,             &
                                        DEGREE,                 &
                                        L_LIMIT,                &
                                        nPROCS_POSEIDON,        &
                                        NUM_CFA_VARS,            &
                                        NUM_SHELLS,             &
                                        NUM_SUBSHELLS,          &
                                        NUM_SUBSHELLS_PER_SHELL,&
                                        NUM_R_ELEMS_PER_BLOCK,  &
                                        NUM_T_ELEMS_PER_BLOCK,  &
                                        NUM_P_ELEMS_PER_BLOCK,  &
                                        NUM_R_ELEMS_PER_SUBSHELL,&
                                        NUM_BLOCK_THETA_ROWS,   &
                                        NUM_BLOCK_PHI_COLUMNS

USE Global_Variables_And_Parameters, &
                                ONLY :  NUM_R_ELEMENTS,                 &
                                        NUM_T_ELEMENTS,                 &
                                        NUM_P_ELEMENTS,                 &
                                        NUM_R_NODES,                    &
                                        VAR_DIM,                    &
                                        Block_Source_E,             &
                                        Block_Source_S,             &
                                        Block_Source_Si,            &
                                        Coefficient_Vector,         &
                                        RHS_Vector,                 &
                                        Block_RHS_Vector,           &
                                        nPROCS_PETSC,               &
                                        myID_Poseidon,              &
                                        myID_Shell,                 &
                                        myID_Petsc,                 &
                                        myID_SUBSHELL,              &
                                        ierr,                       &
                                        PROB_DIM,                   &
                                        Coefficient_Vector,         &
                                        LM_LENGTH,                  &
                                        ULM_LENGTH,                 &
                                        POSEIDON_COMM_WORLD,        &
                                        POSEIDON_COMM_PETSC,        &
                                        Local_Length,               &
                                        ELEM_PROB_DIM,              &
                                        ELEM_PROB_DIM_SQR,          &
                                        BLOCK_PROB_DIM,             &
                                        SUBSHELL_PROB_DIM,          &
                                        BLOCK_ELEM_STF_MATVEC,      &
                                        myShell,                    &
                                        Update_Vector
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscpc.h>

use petscksp
use petscsys
use petscmat
use petscvec
use petscpc


IMPLICIT NONE




CONTAINS



!+201+###########################################################################!
!                                                                                !
!                  CFA_Solve                                                     !
!                                                                                !
!################################################################################!
SUBROUTINE PETSC_Distributed_Solve()


double precision info(MAT_INFO_SIZE)

PetscInt                        ::   my_m, my_n, my_xa, my_ya,my_xb,my_yb
PetscInt,ALLOCATABLE            ::   Index_Array_RHS(:)
PetscInt                        ::   Index_Array_MAT(0:ELEM_PROB_DIM-1)
PetscInt,ALLOCATABLE            ::   Index_Array_Sol(:)
PetscInt                        ::   Start_here, Mid_Here, End_here
PetscInt                        ::   i, j, re
PetscInt                        ::   i_here, j_here
PetscInt                        ::   length_val_a, length_val_b
PetscInt,ALLOCATABLE            ::   on_diag_nnz_a(:)
PetscInt,ALLOCATABLE            ::   off_diag_nnz_a(:)


PetscInt                        ::   Iter_Count, maxits
double precision                ::   rtol, abstol, dtol
              

PetscScalar                     ::   Petsc_Matvec(0:NUM_R_ELEMS_PER_SUBSHELL*ELEM_PROB_DIM_SQR-1)
PetscScalar                     ::   Petsc_MatValues(0:ELEM_PROB_DIM-1)

PetscScalar                     ::   Petsc_SourceValues(0:Local_Length-1)

PetscScalar                     ::   Sol_Vec(0:Local_Length-1)
PetscInt                        ::   ncols

Vec                             ::   Source, Solution
Mat                             ::   Stiffness
KSP                             ::   ksp
PC                              ::   pc

PetscErrorCode                  ::   jerr

INTEGER                         ::   Block_Re, Global_re
INTEGER                         ::   Gm, Gn, Lm, Ln, SRC_G, SRC_L, SOL_G, SOL_L

REAL(KIND = idp)                ::   timea, timeb, timec

INTEGER                         ::   Start_Here_b, End_Here_b


INTEGER                         ::   Local_Length_RHS


IF ( POSEIDON_COMM_PETSC .NE. MPI_COMM_NULL ) THEN


    CALL MPI_COMM_DUP(POSEIDON_COMM_PETSC, PETSC_COMM_WORLD, jerr)


    CALL PetscInitialize(PETSC_NULL_CHARACTER, jerr)
    IF (jerr .ne. 0 ) THEN
        PRINT*,"Unable to Initialize PETSc"
    END IF



!
! Calculate dimension of submatrices
!
    Length_Val_a = ELEM_PROB_DIM - ULM_LENGTH
    Length_Val_b = ULM_LENGTH

    


    ALLOCATE( Index_Array_Sol(0:Local_Length-1) )
    ALLOCATE( on_diag_nnz_a(0:Local_Length-1), off_diag_nnz_a(0:Local_Length-1) )



   



    !
    !   Create Petsc Vectors
    !
    CALL VecCreateMPI(POSEIDON_COMM_PETSC, Local_Length, PROB_DIM, Source, jerr)
    CALL VecSetFromOptions( Source, jerr)
    CALL VecDuplicate(Source, Solution, jerr)



    Local_Length_RHS = SUBSHELL_PROB_DIM
    ALLOCATE( Index_Array_RHS(0:Local_Length_RHS-1) )



    start_here=myID_PETSc*NUM_R_ELEMS_PER_SUBSHELL*DEGREE*ULM_LENGTH
    end_here=start_here + Local_Length_RHS - 1
    Index_Array_RHS =  (/(i,i=start_here,end_here,1)/)

    start_here=myID_SHELL*NUM_R_ELEMS_PER_SUBSHELL*DEGREE*ULM_LENGTH
    end_here=start_here + Local_Length_RHS - 1









    CALL VecSetValues(Source, Local_Length_RHS, Index_Array_RHS,              &
                      Block_RHS_VECTOR(START_HERE:END_HERE), ADD_VALUES, jerr)




!    PRINT*,"----------------------  Assemble Vector  ---------------------------",myID_Poseidon

    CALL VecAssemblyBegin(Source, jerr)
    CALL VecAssemblyEnd(Source, jerr)

  
!    PRINT*,"----------------------  Vector Assembled  ---------------------------",myID_Poseidon
!    CALL SLEEP(1)








!
! Create Matrix 
!
    CALL MatCreate(POSEIDON_COMM_PETSC, Stiffness, jerr)
    CALL MatSetType(Stiffness, MATAIJ, jerr)

    CALL MatSetSizes(Stiffness, LOCAL_LENGTH,LOCAL_LENGTH , PROB_DIM, PROB_DIM, jerr)

    CALL MatSetFromOptions(Stiffness, jerr)




!
!   Set Preallocation Arrays
!
    ! Handle First Element of Block
    IF ( NUM_R_ELEMS_PER_SUBSHELL .NE. 1 ) THEN
        IF (myID_PETSc == 0 ) THEN

            on_diag_nnz_a(0:Length_Val_a-1) = ELEM_PROB_DIM
            off_diag_nnz_a(0:Length_Val_a-1) = 0

        ELSE

            on_diag_nnz_A(0:Length_Val_a-1) = ELEM_PROB_DIM

            off_diag_nnz_a(0:Length_Val_b-1) = Length_Val_a
            off_diag_nnz_a(Length_Val_b:Length_Val_a-1) = 0

        END IF


        ! Handle Middle Elements of Block
        DO re = 1,NUM_R_ELEMS_PER_SUBSHELL-2

            Start_Here = re*Length_Val_a
            Mid_Here = Start_Here + Length_Val_b
            End_Here = (re+1)*Length_Val_a

            on_diag_nnz_a(Start_Here:Mid_Here-1) = ELEM_PROB_DIM + Length_Val_a
            on_diag_nnz_a(Mid_Here:End_Here-1) = ELEM_PROB_DIM

            off_diag_nnz_a(Start_Here:End_Here-1) = 0

        END DO


        ! Handle Last Element of Block
        IF ( myID_PETSc == NUM_SUBSHELLS-1 ) THEN

            Start_Here = (NUM_R_ELEMS_PER_SUBSHELL-1)*Length_Val_a
            Mid_Here = Start_Here + Length_Val_b
            End_Here = Start_Here + ELEM_PROB_DIM



            on_diag_nnz_a(Start_Here:Mid_Here-1) = ELEM_PROB_DIM + Length_Val_a
            on_diag_nnz_a(Mid_Here:End_Here-1) = ELEM_PROB_DIM

            off_diag_nnz_a(Start_Here:End_Here-1) = 0



        ELSE

            Start_Here = (NUM_R_ELEMS_PER_SUBSHELL-1)*Length_Val_a
            Mid_Here = Start_Here + Length_Val_b
            End_Here = NUM_R_ELEMS_PER_SUBSHELL*Length_Val_a

            on_diag_nnz_a(Start_Here:Mid_Here-1) = 2*Length_Val_a
            on_diag_nnz_a(Mid_Here:End_Here-1) = Length_Val_a

            off_diag_nnz_a(Start_Here:End_Here-1) = Length_Val_b

        END IF

     ELSE 


        IF ( myID_PETSC == NUM_SUBSHELLS-1) THEN
            on_diag_nnz_a(0:ELEM_PROB_DIM-1) = ELEM_PROB_DIM
            off_diag_nnz_a = 0
            off_diag_nnz_a(0:ULM_LENGTH-1) = ULM_LENGTH
        ELSE IF ( myID_PETSC == 0 ) THEN
       
            on_diag_nnz_a(0:Length_Val_a-1) = Length_Val_a
       
            off_diag_nnz_a(0:Length_Val_a-1) = ULM_Length
       ELSE 
            on_diag_nnz_a(0:Length_Val_a-1) = Length_Val_a
            off_diag_nnz_a(0:Length_Val_a-1) = 2*ULM_LENGTH

       END IF

     END IF





    IF (nPROCS_PETSC == 1 ) THEN

        CALL MatSeqAIJSetPreallocation(Stiffness, 0, on_diag_nnz_a(0:Local_Length-1), jerr)

    ELSE 
       
        CALL MatMPIAIJSetPreallocation(Stiffness, 0, on_diag_nnz_a(0:Local_Length-1),       &
                                                  0, off_diag_nnz_a(0:Local_Length-1), jerr)
    END IF 









    !
    !   Set Matrix Values
    !
    timea = MPI_WTime()
    DO re = 0,NUM_R_ELEMS_PER_SUBSHELL-1

!
!       Calc Shell_re from subshell_re
!

        Block_Re = myID_Shell*NUM_R_ELEMS_PER_SUBSHELL + re
        Global_Re = myID_PETSc*NUM_R_ELEMS_PER_SUBSHELL + re

        start_here=(myID_PETSc*NUM_R_ELEMS_PER_SUBSHELL + re)*DEGREE*ULM_LENGTH
        end_here=start_here+Elem_Prob_Dim-1
        Index_Array_MAT(0:ELEM_PROB_DIM-1) =  (/(i,i=start_here,end_here,1)/)






        DO i = 0,ELEM_PROB_DIM-1

            i_here = start_here+i


            CALL MatSetValues(  Stiffness, 1, i_here, ELEM_PROB_DIM, Index_Array_MAT,                     &
                                BLOCK_ELEM_STF_MATVEC(i*ELEM_PROB_DIM:(i+1)*ELEM_PROB_DIM-1, Block_re),   &
                                ADD_VALUES, jerr)

 
        END DO

    END DO
    timeb = MPI_WTime()
    IF (myID_Poseidon == 0 ) THEN
        PRINT*,"Set Matrix Values Time ",timeb-timea
    END IF


















!    PRINT*,"***********************  Assemble matrix  ***************************",myID_Poseidon



    !
    !  Tell Petsc to assemble matrix
    ! 
    CALL MatAssemblyBegin(Stiffness, MAT_FINAL_ASSEMBLY,jerr)
    CALL MatAssemblyEND(Stiffness, MAT_FINAL_ASSEMBLY,jerr)

!    CALL SLEEP(3)
!    PRINT*,"++++++++++++++++++++++  Matrix Assembled +++++++++++++++++++++++++++",myID_Poseidon
!    CALL MPI_BARRIER(POSEIDON_COMM_PETSC, jerr)











!     start_here=myID_PETSc*NUM_R_ELEMS_PER_SUBSHELL*DEGREE*ULM_LENGTH
!     end_here=start_here + Local_Length - 1
!     Index_Array_RHS =  (/(i,i=start_here,end_here,1)/)
!     CALL VecGetValues(Source, Local_Length, Index_Array_RHS, Petsc_SourceValues,jerr)
!    DO i = 0,NUM_SUBSHELLS-1
!       IF ( myID_PETSC == 0) THEN
!           PRINT*,"Get RHS_Vector (PETSC)",i,Start_Here,End_Here
!           PRINT*,Petsc_SourceValues
!           PRINT*," "
!       END IF
!       CALL SLEEP(1)
!       CALL MPI_BARRIER(POSEIDON_COMM_PETSC, jerr)
!    END DO
!    PRINT*," "



!    DO re = 0,NUM_R_ELEMS_PER_SUBSHELL-1
!
!!
!!       Calc Shell_re from subshell_re
!!
!
!        Block_Re = myID_Shell*NUM_R_ELEMS_PER_SUBSHELL + re
!        Global_Re = myID_PETSc*NUM_R_ELEMS_PER_SUBSHELL + re
!
!        start_here=(myID_PETSc*NUM_R_ELEMS_PER_SUBSHELL + re)*DEGREE*ULM_LENGTH
!        end_here=start_here+Elem_Prob_Dim-1
!        Index_Array_MAT(0:ELEM_PROB_DIM-1) =  (/(i,i=start_here,end_here,1)/)
!
!        DO i = 0,ELEM_PROB_DIM-1
!            i_here = start_here+i
!
!            CALL MatGetValues( Stiffness, 1, i_here, ELEM_PROB_DIM, Index_Array_Mat, Petsc_MatValues, jerr )
!
!            DO j = 0,NUM_SUBSHELLS-1
!              IF ( myID_PETSC == j) THEN
!                 PRINT*,"MatGetValues myID",j,"Re",Block_re," row",i,i_here
!                 PRINT*,Petsc_MatValues
!                 PRINT*,"========================================================"
!              END IF
!              CALL MPI_BARRIER(POSEIDON_COMM_PETSC, jerr)
!            END DO
!            CALL SLEEP(1)
!            PRINT*," "
!
!        END DO
!
!
!    END DO
































    !   Create Solver
    !
    CALL KSPCreate(POSEIDON_COMM_PETSC,ksp, jerr)
    CALL KSPSetOperators(ksp, Stiffness, Stiffness, jerr)
 
    CALL KSPGetPC(ksp, pc, jerr)

!    CALL PCSetType(pc, PCCHOLESKY, jerr)


    CALL KSPSetType(ksp, KSPBCGS, jerr)
!    CALL KSPSetType(ksp, KSPPREONLY, jerr)
   
    CALL KSPSetFromOptions(ksp, jerr)


!    PRINT*,"BEFORE KSPSolve",myID_Poseidon
!    CALL MPI_BARRIER(POSEIDON_COMM_PETSC, jerr)


    !
    !  Solve the System
    !
 
    timea = MPI_WTime()
    CALL KSPSolve(ksp, Source, Solution, jerr)
    timeb = MPI_WTime()

    CALL KSPGetIterationNumber(ksp, Iter_Count, jerr)
    CALL KSPGetTolerances(ksp, rtol, abstol, dtol, maxits,jerr)
    IF (myID_Poseidon == 0) THEN
!        PRINT*,"Solver Type : KSPPREONLY, PCCHOLESKY"
        PRINT*,"Solver Type : KSPBCGS, PCCHOLESKY"
        PRINT*,"Solve Time ",timeb-timea
        PRINT*,"Solve Iterations :",Iter_Count
        PRINT*,"Solve Tolerances :",rtol, abstol, dtol, maxits
    END IF
!    CALL SLEEP(3)


    !
    !   Retrieve Solution
    !
    Start_Here = myID_PETSc*Length_Val_a*NUM_R_ELEMS_PER_SUBSHELL
    End_Here = Start_Here + Local_Length-1
    Index_Array_Sol = (/ (i, i=Start_Here,End_Here,1)/)


!    PRINT*,"BEFORE VecGetValues",myID_Poseidon
    CALL VecGetValues(Solution, Local_Length, Index_Array_Sol, Sol_Vec, jerr )

    Update_Vector = 0.0_idp
    Update_Vector(Start_Here:End_Here) = Sol_Vec


!    PRINT*,"In Poseidon_Petsc_Solver"
!    PRINT*,"Update_Vector"
!    PRiNT*,REAL(Update_Vector, KIND = idp)


    !
    !  Destroy Vectors, Matrix, and KSP
    ! 
    CALL VecDestroy(Source, jerr)
    CALL VecDestroy(Solution, jerr)
    CALL MatDestroy(Stiffness, jerr)
    CALL KSPDestroy(ksp, jerr)



    !
    ! Close PETSC
    !
    CALL PetscFinalize(jerr)

END IF




END SUBROUTINE PETSC_Distributed_Solve

























END MODULE Poseidon_Petsc_Solver
