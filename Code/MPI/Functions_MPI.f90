   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Functions_MPI                                                                !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the subroutine responsible for creating the various MPI groups,    !##!
!##!    and communicators to facilitate data (re)distribution and solver use.       !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!


USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Variables_MPI, &
            ONLY :  NUM_SHELLS,                 &
                    NUM_SUBSHELLS,              &
                    NUM_SUBSHELLS_PER_SHELL,    &
                    NUM_BLOCKS_PER_SHELL,       &
                    nPROCS_POSEIDON,            &
                    ierr,                       &
                    myID_Poseidon,              &
                    myID_Shell,                 &
                    myID_Dist,                  &
                    myID_PETSc,                 &
                    myID_SUBSHELL,              &
                    nPROCS_PETSc,               &
                    myShell,                    &
                    nPROCS_SHELL,               &
                    Local_Length,               &
                    POSEIDON_COMM_WORLD,        &
                    POSEIDON_COMM_DIST,         &
                    POSEIDON_COMM_SHELL,        &
                    POSEIDON_COMM_PETSC

USE Variables_Derived, &
            ONLY :  ULM_Length,                 &
                    SubShell_Prob_Dim


USE mpi


IMPLICIT NONE

CONTAINS


! !+101+############################################################################!
!!                                                                                   !
!!       CREATE_POSEIDON_COMMUNICATORS                                               !
!!                                                                                   !
! !#################################################################################!
!SUBROUTINE CREATE_POSEIDON_COMMUNICATORS()
!
!
!INTEGER, ALLOCATABLE, DIMENSION(:)                  ::  Block_Array
!INTEGER, ALLOCATABLE, DIMENSION(:)                  ::  Shell_Array
!INTEGER, ALLOCATABLE, DIMENSION(:)                  ::  Workers_Array
!
!
!INTEGER                                             ::  nPROCS
!INTEGER                                             ::  i, j
!
!INTEGER                                             ::  key, color
!
!INTEGER                                             ::  MPI_GROUP_WORLD,            &
!                                                        POSEIDON_GROUP_WORLD,       &
!                                                        POSEIDON_GROUP_INTERSHELL
!INTEGER                                             ::  POSEIDON_GROUP_PETSC
!
!INTEGER                                             ::  here
!
!
!
!ALLOCATE(   Block_Array(0:NUM_BLOCKS_PER_SHELL-1)   )
!ALLOCATE(   Shell_Array(0:NUM_SUBSHELLS-1)             )
!ALLOCATE(   Workers_Array(0:nPROCS_POSEIDON-1)      )
!
!
!myID_Poseidon = -1
!CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nPROCS,ierr)
!
!
!!
!!   Create POSEIDON_COMM_WORLD
!!
!IF ( nPROCS_POSEIDON < nPROCS ) THEN
!
!
!
!!
!!   Create POSEIDON_COMM_WORLD
!!       This COMMUNICATOR contains all processes assigned to the running
!!           of the Poseidon code.
!!
!    CALL MPI_Comm_group(MPI_COMM_WORLD, MPI_GROUP_WORLD, ierr)
!
!
!    Workers_Array = (/(i, i=0,nPROCS_POSEIDON-1, 1)/)
!
!
!    CALL MPI_Group_incl(    MPI_GROUP_WORLD,            &
!                            nPROCS_POSEIDON,            &
!                            Workers_Array,              &
!                            POSEIDON_GROUP_WORLD,       &
!                            ierr                        )
!
!
!    CALL MPI_Comm_create_group( MPI_COMM_WORLD,             &
!                                POSEIDON_GROUP_WORLD,       &
!                                0,                          &
!                                POSEIDON_COMM_WORLD,        &
!                                ierr                        )
!
!
!    IF ( MPI_COMM_NULL .NE. POSEIDON_COMM_WORLD ) THEN
!
!        CALL MPI_COMM_RANK(POSEIDON_COMM_WORLD, myID_Poseidon, ierr)
!
!    END IF
!
!
!ELSE IF ( nPROCS_POSEIDON > nPROCS ) THEN
!
!    PRINT*," Poseidon expects more processes than in MPI_COMM_WORLD "
!    PRINT*," Expected : ",nPROCS_POSEIDON," Found : ", nPROCS
!
!ELSE
!    POSEIDON_COMM_WORLD = MPI_COMM_WORLD
!    CALL MPI_Comm_group(POSEIDON_COMM_WORLD, POSEIDON_GROUP_WORLD, ierr)
!    CALL MPI_COMM_RANK(POSEIDON_COMM_WORLD, myID_Poseidon, ierr)
!
!END IF
!
!
!
!
!
!
!!
!!   Create POSEIDON_COMM_DIST
!!       This communicator divides MPI_COMM_WORLD (as evenly as possible) so that
!!       each group contains one (and only one) member of POSEIDON_COMM_PETSC(see below).
!!       This will allow distribution of the solution to all processes.
!!
!
!color = MOD(myID_Poseidon, NUM_SHELLS)
!key = myID_Poseidon/NUM_SHELLS
!
!
!
!CALL MPI_Comm_Split(MPI_COMM_WORLD, color, key, POSEIDON_COMM_DIST, ierr)
!CALL MPI_COMM_RANK(POSEIDON_COMM_DIST, myID_DIST, ierr)
!
!
!
!
!
!
!
!
!
!
!
!IF ( POSEIDON_COMM_WORLD .NE. MPI_COMM_NULL ) THEN
!
!        myShell = myID_Poseidon/NUM_BLOCKS_PER_SHELL
!
!
!        !   Create Communicator for each Shell
!        !       Created such that NUM_BLOCKS_PER_SHELL consecutive processes are
!        !       assigned to each shell.
!        !
!        Block_Array = (/ (i+NUM_BLOCKS_PER_SHELL*myShell, i=0,NUM_BLOCKS_PER_SHELL-1, 1)/)
!
!
!
!        CALL MPI_Group_incl(    POSEIDON_GROUP_WORLD,       &
!                                NUM_BLOCKS_PER_SHELL,       &
!                                Block_Array,                &
!                                POSEIDON_GROUP_INTERSHELL,  &
!                                ierr                        )
!
!        CALL MPI_Comm_create_group( POSEIDON_COMM_WORLD,            &
!                                    POSEIDON_GROUP_INTERSHELL,      &
!                                    0,                              &
!                                    POSEIDON_COMM_SHELL,            &
!                                    ierr                            )
!
!
!        CALL MPI_COMM_RANK(POSEIDON_COMM_SHELL, myID_Shell, ierr)
!        CALL MPI_COMM_SIZE(POSEIDON_COMM_SHELL, nPROCS_SHELL, ierr)
!
!
!
!
!        !
!        !   Create Communicator for Distributed Solver
!        !    Currently done so that Rank 0 of each shell is included.
!        !    Inner most shell is rank 0, second inner most is rank 1, etc.
!        !
!        DO i = 0,NUM_SHELLS-1
!
!            here = i*NUM_SUBSHELLS_PER_SHELL
!            Shell_Array(here:here+NUM_SUBSHELLS_PER_SHELL-1) = i*NUM_BLOCKS_PER_SHELL          &
!                                                   + (/ (j,j=0,NUM_SUBSHELLS_PER_SHELL-1,1) /)
!        END DO
!
!        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!        CALL MPI_Group_incl(    POSEIDON_GROUP_WORLD,       &
!                                NUM_SUBSHELLS,              &
!                                Shell_Array,                &
!                                POSEIDON_GROUP_PETSC,       &
!                                ierr                        )
!
!        CALL MPI_Comm_create_group( POSEIDON_COMM_WORLD,            &
!                                    POSEIDON_GROUP_PETSC,           &
!                                    0,                              &
!                                    POSEIDON_COMM_PETSC,            &
!                                    ierr                            )
!
!        IF ( POSEIDON_COMM_PETSC .NE. MPI_COMM_NULL ) THEN
!
!            CALL MPI_COMM_RANK(POSEIDON_COMM_PETSC, myID_PETSc, ierr)
!            CALL MPI_COMM_SIZE(POSEIDON_COMM_PETSC, nPROCS_PETSc, ierr)
!
!            myID_SubShell = MOD(myID_PETSC, NUM_SUBSHELLS_PER_SHELL)
!
!        END IF
!
!
!        IF ( myID_PETSc == NUM_SUBSHELLS-1) THEN
!
!          Local_Length = SUBSHELL_PROB_DIM
!
!        ELSE
!
!          Local_Length = SUBSHELL_PROB_DIM - ULM_LENGTH
!
!        END IF
!
!!        PRINT*,"myID",myID,"myID_PETSC",myID_PETSC,"myShell",myShell,"myID_SubShell",myID_SubShell
!
!!       PRINT*,"myID",myID,"myID_Poseidon",myID_Poseidon,"myShell",myShell,"myID_Shell",myID_Shell,"myID_Petsc",myID_PETSc
!
!
!END IF
!
!
!
!
!!CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!
!
!
!END SUBROUTINE CREATE_POSEIDON_COMMUNICATORS




 !+101+############################################################################!
!                                                                                   !
!       CREATE_POSEIDON_COMMUNICATORS                                               !
!                                                                                   !
 !#################################################################################!
SUBROUTINE CREATE_POSEIDON_COMMUNICATORS()


INTEGER, ALLOCATABLE, DIMENSION(:)                  ::  Block_Array
INTEGER, ALLOCATABLE, DIMENSION(:)                  ::  Shell_Array
INTEGER, ALLOCATABLE, DIMENSION(:)                  ::  Workers_Array


INTEGER                                             ::  nPROCS
INTEGER                                             ::  i, j

INTEGER                                             ::  key, color

INTEGER                                             ::  MPI_GROUP_WORLD,            &
                                                        POSEIDON_GROUP_WORLD,       &
                                                        POSEIDON_GROUP_INTERSHELL
INTEGER                                             ::  POSEIDON_GROUP_PETSC

INTEGER                                             ::  here



ALLOCATE(   Block_Array(0:NUM_BLOCKS_PER_SHELL-1)   )
ALLOCATE(   Shell_Array(0:NUM_SUBSHELLS-1)             )
ALLOCATE(   Workers_Array(0:nPROCS_POSEIDON-1)      )


myID_Poseidon = -1
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nPROCS,ierr)

!
!   Create POSEIDON_COMM_WORLD
!
IF ( nPROCS_POSEIDON < nPROCS ) THEN

!
!   Create POSEIDON_COMM_WORLD
!       This COMMUNICATOR contains all processes assigned to the running
!           of the Poseidon code.
!
    CALL MPI_Comm_group(MPI_COMM_WORLD, MPI_GROUP_WORLD, ierr)


    Workers_Array = (/(i, i=0,nPROCS_POSEIDON-1, 1)/)


    CALL MPI_Group_incl(    MPI_GROUP_WORLD,            &
                            nPROCS_POSEIDON,            &
                            Workers_Array,              &
                            POSEIDON_GROUP_WORLD,       &
                            ierr                        )


    CALL MPI_Comm_create_group( MPI_COMM_WORLD,             &
                                POSEIDON_GROUP_WORLD,       &
                                0,                          &
                                POSEIDON_COMM_WORLD,        &
                                ierr                        )


    IF ( MPI_COMM_NULL .NE. POSEIDON_COMM_WORLD ) THEN

        CALL MPI_COMM_RANK(POSEIDON_COMM_WORLD, myID_Poseidon, ierr)

    END IF


ELSE IF ( nPROCS_POSEIDON > nPROCS ) THEN

    PRINT*," Poseidon expects more processes than in MPI_COMM_WORLD "
    PRINT*," Expected : ",nPROCS_POSEIDON," Found : ", nPROCS

ELSE
    POSEIDON_COMM_WORLD = MPI_COMM_WORLD
    CALL MPI_Comm_group(POSEIDON_COMM_WORLD, POSEIDON_GROUP_WORLD, ierr)
    CALL MPI_COMM_RANK(POSEIDON_COMM_WORLD, myID_Poseidon, ierr)

END IF




!
!   Create POSEIDON_COMM_DIST
!       This communicator divides MPI_COMM_WORLD (as evenly as possible) so that
!       each group contains one (and only one) member of POSEIDON_COMM_PETSC(see below).
!       This will allow distribution of the solution to all processes.
!

color = MOD(myID_Poseidon, NUM_SHELLS)
key = myID_Poseidon/NUM_SHELLS



CALL MPI_Comm_Split(MPI_COMM_WORLD, color, key, POSEIDON_COMM_DIST, ierr)
CALL MPI_COMM_RANK(POSEIDON_COMM_DIST, myID_DIST, ierr)











IF ( POSEIDON_COMM_WORLD .NE. MPI_COMM_NULL ) THEN

        myShell = myID_Poseidon/NUM_BLOCKS_PER_SHELL


        !   Create Communicator for each Shell
        !       Created such that NUM_BLOCKS_PER_SHELL consecutive processes are
        !       assigned to each shell.
        !
        Block_Array = (/ (i+NUM_BLOCKS_PER_SHELL*myShell, i=0,NUM_BLOCKS_PER_SHELL-1, 1)/)



        CALL MPI_Group_incl(    POSEIDON_GROUP_WORLD,       &
                                NUM_BLOCKS_PER_SHELL,       &
                                Block_Array,                &
                                POSEIDON_GROUP_INTERSHELL,  &
                                ierr                        )

        CALL MPI_Comm_create_group( POSEIDON_COMM_WORLD,            &
                                    POSEIDON_GROUP_INTERSHELL,      &
                                    0,                              &
                                    POSEIDON_COMM_SHELL,            &
                                    ierr                            )


        CALL MPI_COMM_RANK(POSEIDON_COMM_SHELL, myID_Shell, ierr)
        CALL MPI_COMM_SIZE(POSEIDON_COMM_SHELL, nPROCS_SHELL, ierr)


   

        !
        !   Create Communicator for Distributed Solver
        !    Currently done so that Rank 0 of each shell is included.
        !    Inner most shell is rank 0, second inner most is rank 1, etc.
        !
        DO i = 0,NUM_SHELLS-1

            here = i*NUM_SUBSHELLS_PER_SHELL
            Shell_Array(here:here+NUM_SUBSHELLS_PER_SHELL-1) = i*NUM_BLOCKS_PER_SHELL          &
                                                   + (/ (j,j=0,NUM_SUBSHELLS_PER_SHELL-1,1) /)
        END DO

        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

        CALL MPI_Group_incl(    POSEIDON_GROUP_WORLD,       &
                                NUM_SUBSHELLS,              &
                                Shell_Array,                &
                                POSEIDON_GROUP_PETSC,       &
                                ierr                        )

        CALL MPI_Comm_create_group( POSEIDON_COMM_WORLD,            &
                                    POSEIDON_GROUP_PETSC,           &
                                    0,                              &
                                    POSEIDON_COMM_PETSC,            &
                                    ierr                            )

        IF ( POSEIDON_COMM_PETSC .NE. MPI_COMM_NULL ) THEN

            CALL MPI_COMM_RANK(POSEIDON_COMM_PETSC, myID_PETSc, ierr)
            CALL MPI_COMM_SIZE(POSEIDON_COMM_PETSC, nPROCS_PETSc, ierr)

            myID_SubShell = MOD(myID_PETSC, NUM_SUBSHELLS_PER_SHELL)

        END IF


        IF ( myID_PETSc == NUM_SUBSHELLS-1) THEN

          Local_Length = SUBSHELL_PROB_DIM
          
        ELSE

          Local_Length = SUBSHELL_PROB_DIM - ULM_LENGTH

        END IF

!        PRINT*,"myID",myID,"myID_PETSC",myID_PETSC,"myShell",myShell,"myID_SubShell",myID_SubShell

!       PRINT*,"myID",myID,"myID_Poseidon",myID_Poseidon,"myShell",myShell,"myID_Shell",myID_Shell,"myID_Petsc",myID_PETSc


END IF




!CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)




END SUBROUTINE CREATE_POSEIDON_COMMUNICATORS
















END MODULE Functions_MPI
