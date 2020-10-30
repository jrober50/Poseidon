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





END MODULE Initialization_Mesh

