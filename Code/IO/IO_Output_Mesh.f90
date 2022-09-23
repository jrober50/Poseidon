   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE IO_Output_Mesh_Module                                                 !##!
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
            ONLY : idp

USE Poseidon_Numbers_Module, &
            ONLY : pi

USE Poseidon_Parameters, &
            ONLY :  DEGREE

USE Variables_Derived, &
            ONLY :  Num_R_Nodes


USE Variables_Mesh, &
            ONLY :  NUM_R_ELEMENTS,         &
                    rlocs,                  &
                    tlocs,                  &
                    R_Inner,                &
                    R_Outer

USE Variables_IO, &
            ONLY :  File_Suffix

USE Variables_FEM_Module, &
            ONLY :  FEM_Node_xlocs

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature_Locations,     &
                    Initialize_LGL_Quadrature_Locations

USE Poseidon_File_Routines_Module, &
            ONLY :  Open_New_File


IMPLICIT NONE


CONTAINS


!+202+###########################################################################!
!                                                                                !
!                   OUTPUT_JACOBIAN_MATRIX                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Output_Mesh( Mesh, row_in, flag )


INTEGER,                        INTENT(IN)                  ::  row_in
REAL(idp), DIMENSION(1:Row_In), INTENT(IN)               ::  Mesh


CHARACTER(LEN = 1 ), INTENT(IN), OPTIONAL                   ::  flag

CHARACTER(LEN = 300)                                        ::  FILE_NAME
CHARACTER(LEN = 300)                                        ::  FILE_NAMEb
CHARACTER(LEN = 40)                                         ::  fmt


INTEGER                                                     ::  rows

INTEGER                                                     ::  FILE_ID
INTEGEr                                                     ::  i

CHARACTER(LEN = 29)        :: Poseidon_Mesh_Dir      = "Poseidon_Output/Objects/Mesh/"

100 FORMAT (A,A,A,A,A,A)
101 FORMAT (A,A,A,A)

fmt = '(ES24.16E3,SP,ES24.16E3,"i")'


rows = size( Mesh, 1 )

IF ( present(flag) ) THEN
    WRITE(FILE_NAMEb,100) Poseidon_Mesh_Dir,"Mesh_Radial_Dim_",trim(File_Suffix),"_",flag,".out"
ELSE
    WRITE(FILE_NAMEb,101) Poseidon_Mesh_Dir,"Mesh_Radial_Dim_",trim(File_Suffix),".out"
END IF


CALL OPEN_NEW_FILE( trim(FILE_NAMEb), FILE_ID, 300 )
WRITE(FILE_ID,*) rows
CLOSE(FILE_ID)



IF ( present(flag) ) THEN
    WRITE(FILE_NAME,100) Poseidon_Mesh_Dir,"Mesh_Radial_Loc_",trim(File_Suffix),"_",flag,".out"
ELSE
    WRITE(FILE_NAME,100) Poseidon_Mesh_Dir,"Mesh_Radial_Loc_",trim(File_Suffix),".out"
END IF

CALL OPEN_NEW_FILE( trim(FILE_NAME), FILE_ID, 300 )
DO i = 1,rows
    WRITE(FILE_ID,fmt) Mesh(i)
END DO
CLOSE( FILE_ID )



END SUBROUTINE Output_Mesh




!+202+###########################################################################!
!                                                                                !
!                   Output_Nodal_Mesh                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Output_Nodal_Mesh( Mesh, row_in, flag )


INTEGER,                        INTENT(IN)                  ::  row_in
REAL(idp), DIMENSION(1:Row_In), INTENT(IN)               ::  Mesh


CHARACTER(LEN = 1 ), INTENT(IN), OPTIONAL                   ::  flag

CHARACTER(LEN = 300)                                        ::  FILE_NAME
CHARACTER(LEN = 300)                                        ::  FILE_NAMEb
CHARACTER(LEN = 40)                                         ::  fmt

REAL(KIND = idp)                                            ::  DROT
REAL(KIND = idp), DIMENSION(0:Degree)                       ::  Cur_R_Locs

INTEGER                                                     ::  rows

INTEGER                                                     ::  FILE_ID
INTEGER                                                     ::  re, d

CHARACTER(LEN = 29)        :: Poseidon_Mesh_Dir      = "Poseidon_Output/Objects/Mesh/"

100 FORMAT (A,A,A,A,A,A)
101 FORMAT (A,A,A,A)

fmt = '(ES24.16E3,SP,ES24.16E3,"i")'


rows = size( Mesh, 1 )

IF ( present(flag) ) THEN
    WRITE(FILE_NAMEb,100) Poseidon_Mesh_Dir,"Mesh_Nodal_Dim_",trim(File_Suffix),"_",flag,".out"
ELSE
    WRITE(FILE_NAMEb,101) Poseidon_Mesh_Dir,"Mesh_Nodal_Dim_",trim(File_Suffix),".out"
END IF


CALL OPEN_NEW_FILE( trim(FILE_NAMEb), FILE_ID, 300 )
WRITE(FILE_ID,*) Num_R_Nodes
CLOSE(FILE_ID)



IF ( present(flag) ) THEN
    WRITE(FILE_NAME,100) Poseidon_Mesh_Dir,"Mesh_Nodal_Loc_",trim(File_Suffix),"_",flag,".out"
ELSE
    WRITE(FILE_NAME,100) Poseidon_Mesh_Dir,"Mesh_Nodal_Loc_",trim(File_Suffix),".out"
END IF

CALL OPEN_NEW_FILE( trim(FILE_NAME), FILE_ID, 300 )
WRITE(File_ID,fmt) Mesh(1)
DO re = 1,NUM_R_ELEMENTS

    DROT = 0.50_idp*(Mesh(re+1) - Mesh(re))
    CUR_R_LOCS(:) = DROT * (FEM_Node_xlocs(:)+1.0_idp) + Mesh(re)

    DO d = 1,Degree
        WRITE(FILE_ID,fmt) Cur_R_Locs(d)
    END DO
END DO
CLOSE( FILE_ID )



END SUBROUTINE Output_Nodal_Mesh



END MODULE IO_Output_Mesh_Module
