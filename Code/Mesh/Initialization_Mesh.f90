   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Initialization_Mesh                                                          !##!
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

USE Poseidon_Parameters, &
                    ONLY :  Verbose_Flag

USE Variables_Mesh, &
                    ONLY :  Num_R_Elements,         &
                            Num_T_Elements,         &
                            Num_P_Elements,         &
                            Num_Loc_R_Elements,     &
                            Num_Loc_T_Elements,     &
                            Num_Loc_P_Elements,     &
                            R_Inner,                &
                            R_Outer,                &
                            R_Coarsen_Factor,       &
                            T_Coarsen_Factor,       &
                            P_Coarsen_Factor,       &
                            rlocs,                  &
                            tlocs,                  &
                            plocs,                  &
                            drlocs,                 &
                            dtlocs,                 &
                            dplocs,                 &
                            Radial_Mesh_Set_Flag,   &
                            Theta_Mesh_Set_Flag,    &
                            Phi_Mesh_Set_Flag,      &
                            locs_Set,               &
                            dlocs_Set

   
USE Functions_Mesh, &
                    ONLY :  Generate_Defined_Coarse_Mesh


IMPLICIT NONE

CONTAINS

 !+101+################################################################################!
!                                                                                       !
!       Initialize_Mesh                                                                 !
!                                                                                       !
 !#####################################################################################!
SUBROUTINE Initialize_Mesh()


IF ( Verbose_Flag ) THEN
    PRINT*,"-Initializing Mesh. "
END IF

IF ( ( locs_Set(1) .EQV. .FALSE. ) .AND. ( dlocs_Set(1) .EQV. .FALSE. ) ) THEN
    PRINT*,"!*  Fatal Error in Poseidon   *!"
    PRINT*,"    No information passed about mesh. "
    PRINT*,"!*  POSEIDON STOPPING CODE  *!"
    STOP
END IF
IF ( ( locs_Set(2) .EQV. .FALSE. ) .AND. ( dlocs_Set(2) .EQV. .FALSE. ) ) THEN
    PRINT*,"!*  Fatal Error in Poseidon   *!"
    PRINT*,"    No information passed about mesh. "
    PRINT*,"!*  POSEIDON STOPPING CODE  *!"
    STOP
END IF
IF ( ( locs_Set(3) .EQV. .FALSE. ) .AND. ( dlocs_Set(3) .EQV. .FALSE. ) ) THEN
    PRINT*,"!*  Fatal Error in Poseidon   *!"
    PRINT*,"    No information passed about mesh. "
    PRINT*,"!*  POSEIDON STOPPING CODE  *!"
    STOP
END IF

CALL Generate_Defined_Coarse_Mesh(NUM_R_ELEMENTS, NUM_R_ELEMENTS, R_Coarsen_Factor,     &
                                  R_INNER, 1, rlocs, drlocs)

RADIAL_MESH_SET_FLAG = .TRUE.






CALL Generate_Defined_Coarse_Mesh(NUM_T_ELEMENTS, NUM_T_ELEMENTS, T_Coarsen_Factor,     &
                                  0.0_idp, 2, tlocs, dtlocs)
THETA_MESH_SET_FLAG = .TRUE.


CALL Generate_Defined_Coarse_Mesh(NUM_P_ELEMENTS, NUM_P_ELEMENTS, P_Coarsen_Factor,     &
                                  0.0_idp, 3, plocs, dplocs)
PHI_MESH_SET_FLAG = .TRUE.


END SUBROUTINE Initialize_Mesh





 !+501+############################################################################!
!                                                                                   !
!                     OPEN_NEW_FILE                                                 !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OPEN_NEW_FILE(File_Name, File_Number, Suggested_Number)



CHARACTER(LEN = *), INTENT(IN)                          ::  File_Name
INTEGER,            INTENT(INOUT)                       ::  File_Number
INTEGER, OPTIONAL,  INTENT(IN)                          ::  Suggested_Number

INTEGER                                                 ::  Temp_Number
INTEGER                                                 ::  istat = 0
LOGICAL                                                 ::  FLAG, OP
LOGICAL                                                 ::  UNIT_FLAG, NAME_FLAG


UNIT_FLAG = .FALSE.
NAME_FLAG = .FALSE.


!  Assigned an unused number, and assign it to new file
FLAG = .TRUE.
IF ( Present(Suggested_Number) ) THEN
    Temp_Number = Suggested_Number
ELSE
    Temp_Number = 1000
END IF


DO WHILE (FLAG)
    INQUIRE( UNIT = Temp_Number, OPENED = OP )

    IF ( OP ) THEN
        Temp_Number = Temp_Number + 1
    ELSE
        File_Number = Temp_Number
        FLAG = .FALSE.
        UNIT_FLAG = .TRUE.
    END IF
END DO




! Check if file already exists !
!INQUIRE( FILE = File_Name, EXIST = EX )
!IF ( EX ) THEN
!    PRINT*,"File ",File_Name," is already opened"
!
!ELSE
!    PRINT*,File_Name," is not already opened."
!    NAME_FLAG = .TRUE.
!END IF



! Open New File
IF ( UNIT_FLAG  ) THEN

    OPEN( Unit = File_Number, File = File_Name, IOSTAT = istat )
    IF ( istat .NE. 0 ) THEN

        PRINT*,"WARNING: Could not open file at ", File_Name, istat

    END IF
END IF


END SUBROUTINE OPEN_NEW_FILE



END MODULE Initialization_Mesh

