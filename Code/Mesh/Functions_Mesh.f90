   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Functions_Mesh                                                               !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!    This module contains a few basic functions for creating various grid        !##!
!##!    meshes.                                                                     !##!
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

USE Poseidon_Numbers_Module, &
                ONLY :  pi

USE Variables_Mesh, &
                ONLY :  locs_Set,   &
                        dlocs_Set

IMPLICIT NONE

CONTAINS


!+101+##############################################################################!
!                                                                                   !
!           Create_3D_Mesh                                                          !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Calls appropriote mesh generating subroutine based on mesh_type input.          !
!                                                                                   !
!###################################################################################!
SUBROUTINE Create_3D_Mesh(  Mesh_Type,                  &
                            Inner_Radius, Outer_Radius, &
                            nx, ny, nz,                 &
                            x_e, x_c, dx_c,             &
                            y_e, y_c, dy_c,             &
                            z_e, z_c, dz_c,             &
                            Zoom, nc, Core_Radius   )

INTEGER,                            INTENT(IN)              ::  Mesh_Type

REAL(KIND = idp),                   INTENT(IN)              ::  Inner_Radius, Outer_Radius

INTEGER,                            INTENT(IN)              ::  nx, ny, nz

REAL(KIND = idp), DIMENSION(0:nx),  INTENT(OUT)             ::  x_e
REAL(KIND = idp), DIMENSION(1:nx),  INTENT(OUT)             ::  x_c
REAL(KIND = idp), DIMENSION(1:nx),  INTENT(OUT)             ::  dx_c

REAL(KIND = idp), DIMENSION(0:ny),  INTENT(OUT)             ::  y_e
REAL(KIND = idp), DIMENSION(1:ny),  INTENT(OUT)             ::  y_c
REAL(KIND = idp), DIMENSION(1:ny),  INTENT(OUT)             ::  dy_c

REAL(KIND = idp), DIMENSION(0:nz),  INTENT(OUT)             ::  z_e
REAL(KIND = idp), DIMENSION(1:nz),  INTENT(OUT)             ::  z_c
REAL(KIND = idp), DIMENSION(1:nz),  INTENT(OUT)             ::  dz_c

REAL(KIND = idp),                   INTENT(IN), OPTIONAL    ::  Zoom
INTEGER         ,                   INTENT(IN), OPTIONAL    ::  nc
REAL(KIND = idp),                   INTENT(IN), OPTIONAL    ::  Core_Radius





IF ( Mesh_Type == 1 ) THEN

    CALL Create_Uniform_3D_Mesh( Inner_Radius, Outer_Radius,       &
                                 nx, ny, nz,                       &
                                 x_e, x_c, dx_c,                   &
                                 y_e, y_c, dy_c,                   &
                                 z_e, z_c, dz_c                    )



ELSE IF ( Mesh_Type == 2) THEN

    CALL Create_Logarithmic_3D_Mesh( Inner_Radius, Outer_Radius,       &
                                     nx, ny, nz,                       &
                                     x_e, x_c, dx_c,                   &
                                     y_e, y_c, dy_c,                   &
                                     z_e, z_c, dz_c                    )



ELSE IF ( Mesh_Type == 3) THEN

    IF ( PRESENT(nc) .AND. PRESENT(Core_Radius) ) THEN

        CALL Create_Split_3D_Mesh( Inner_Radius, Core_Radius, Outer_Radius,      &
                                   nx, nc, ny, nz,                               &
                                   x_e, x_c, dx_c,                               &
                                   y_e, y_c, dy_c,                               &
                                   z_e, z_c, dz_c                                )

    ELSE

        PRINT*,"!*  Fatal Error in Poseidon *!"
        PRINT*," Error in Create_3D_Mesh() "
        PRINT*,"    Split Mesh (Mesh_Type == 3) requires the specification of : "
        PRINT*," - Elements in the mesh type (nc)           "
        IF ( PRESENT(nc) ) THEN
            PRINT*,"   - SPECIFIED "
        ELSE
            PRINT*,"   - NOT SPECIFIED "
        END IF
        PRINT*," - Radius of the inner mesh (Core_Radius)   "
        IF ( PRESENT(Core_Radius) ) THEN
            PRINT*,"   - SPECIFIED "
        ELSE
            PRINT*,"   - NOT SPECIFIED "
        END IF
        PRINT*,"!*  POSEIDON STOPPING CODE  *!"
        STOP

    END IF

ELSE IF ( Mesh_Type == 4) THEN
    

    IF ( PRESENT(Zoom) ) THEN

        CALL Create_Geometric_3D_Mesh( Inner_Radius, Outer_Radius,         &
                                         nx, ny, nz,                       &
                                         x_e, x_c, dx_c,                   &
                                         y_e, y_c, dy_c,                   &
                                         z_e, z_c, dz_c,                   &
                                         Zoom                               )

    ELSE

        PRINT*,"!*  Fatal Error in Poseidon *!"
        PRINT*," Error in Create_3D_Mesh() "
        PRINT*,"    Geometric Mesh (Mesh_Type == 4) requires the specification of : "
        PRINT*," - Zoom Factor (Zoom)"
        PRINT*,"!*  POSEIDON STOPPING CODE  *!"
        STOP

    END IF


END IF


END SUBROUTINE Create_3D_Mesh




!+201+##############################################################################!
!                                                                                   !
!           Create_Uniform_3D_Mesh                                                  !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Creates a uniform 3D mesh.  Assumes y_range = (0:pi) and z_range = (0:2pi).     !
!                                                                                   !
!###################################################################################!
SUBROUTINE Create_Uniform_3D_Mesh(  Inner_Radius, Outer_Radius,                     &
                                    nx, ny, nz,                                     &
                                    x_e, x_c, dx_c,                                 &
                                    y_e, y_c, dy_c,                                 &
                                    z_e, z_c, dz_c                                  )

REAL(KIND = idp),   INTENT(IN)                          ::  Inner_Radius,           &
                                                            Outer_Radius

INTEGER,            INTENT(IN)                          ::  nx, ny, nz

REAL(KIND = idp), DIMENSION(0:nx),  INTENT(OUT)         ::  x_e
REAL(KIND = idp), DIMENSION(1:nx),  INTENT(OUT)         ::  x_c
REAL(KIND = idp), DIMENSION(1:nx),  INTENT(OUT)         ::  dx_c

REAL(KIND = idp), DIMENSION(0:ny),  INTENT(OUT)         ::  y_e
REAL(KIND = idp), DIMENSION(1:ny),  INTENT(OUT)         ::  y_c
REAL(KIND = idp), DIMENSION(1:ny),  INTENT(OUT)         ::  dy_c

REAL(KIND = idp), DIMENSION(0:nz),  INTENT(OUT)         ::  z_e
REAL(KIND = idp), DIMENSION(1:nz),  INTENT(OUT)         ::  z_c
REAL(KIND = idp), DIMENSION(1:nz),  INTENT(OUT)         ::  dz_c

INTEGER                                                 ::  i


! Create Uniform Radial Mesh !
CALL Create_Uniform_1D_Mesh( Inner_Radius, Outer_Radius, nx, x_e, x_c, dx_c )

! Create Uniform Theta Mesh !
CALL Create_Uniform_1D_Mesh( 0.0_idp, pi, ny, y_e, y_c, dy_c )

! Create Uniform Phi Mesh !
CALL Create_Uniform_1D_Mesh( 0.0_idp, 2.0_idp*pi, nz, z_e, z_c, dz_c )




END SUBROUTINE Create_Uniform_3D_Mesh












!+202+##############################################################################!
!                                                                                   !
!           Create_Logarithmic_Mesh                                                 !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Creates a uniform angular mesh with a logarithmicly spaced radial mesh.         !
!   Assumes y_range = (0:pi) and z_range = (0:2pi).                                 !
!                                                                                   !
!###################################################################################!
SUBROUTINE Create_Logarithmic_3D_Mesh(  Inner_Radius, Outer_Radius,                 &
                                        nx, ny, nz,                                 &
                                        x_e, x_c, dx_c,                             &
                                        y_e, y_c, dy_c,                             &
                                        z_e, z_c, dz_c                              )

REAL(KIND = idp),   INTENT(IN)                          ::  Inner_Radius,           &
                                                            Outer_Radius

INTEGER,            INTENT(IN)                          ::  nx, ny, nz

REAL(KIND = idp), DIMENSION(0:nx),  INTENT(OUT)         ::  x_e
REAL(KIND = idp), DIMENSION(1:nx),  INTENT(OUT)         ::  x_c
REAL(KIND = idp), DIMENSION(1:nx),  INTENT(OUT)         ::  dx_c

REAL(KIND = idp), DIMENSION(0:ny),  INTENT(OUT)         ::  y_e
REAL(KIND = idp), DIMENSION(1:ny),  INTENT(OUT)         ::  y_c
REAL(KIND = idp), DIMENSION(1:ny),  INTENT(OUT)         ::  dy_c

REAL(KIND = idp), DIMENSION(0:nz),  INTENT(OUT)         ::  z_e
REAL(KIND = idp), DIMENSION(1:nz),  INTENT(OUT)         ::  z_c
REAL(KIND = idp), DIMENSION(1:nz),  INTENT(OUT)         ::  dz_c



! Create Logarithmic Radial Mesh !
Call Create_Logarithmic_1D_Mesh( Inner_Radius, Outer_Radius, nx, x_e, x_c, dx_c )

! Create Uniform Theta Mesh !
CALL Create_Uniform_1D_Mesh( 0.0_idp, pi, ny, y_e, y_c, dy_c )

! Create Uniform Phi Mesh !
CALL Create_Uniform_1D_Mesh( 0.0_idp, 2.0_idp*pi, nz, z_e, z_c, dz_c )



END SUBROUTINE Create_Logarithmic_3D_Mesh










!+202+##############################################################################!
!                                                                                   !
!           Create_Geometric_3D_Mesh                                                !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Creates a uniform angular mesh with a Geometrically spaced radial mesh.         !
!   Assumes y_range = (0:pi) and z_range = (0:2pi).                                 !
!                                                                                   !
!###################################################################################!
SUBROUTINE Create_Geometric_3D_Mesh(  Inner_Radius, Outer_Radius,                   &
                                        nx, ny, nz,                                 &
                                        x_e, x_c, dx_c,                             &
                                        y_e, y_c, dy_c,                             &
                                        z_e, z_c, dz_c,                             &
                                        Zoom                                        )

REAL(KIND = idp),   INTENT(IN)                          ::  Inner_Radius,           &
                                                            Outer_Radius

INTEGER,            INTENT(IN)                          ::  nx, ny, nz

REAL(KIND = idp), DIMENSION(0:nx),  INTENT(OUT)         ::  x_e
REAL(KIND = idp), DIMENSION(1:nx),  INTENT(OUT)         ::  x_c
REAL(KIND = idp), DIMENSION(1:nx),  INTENT(OUT)         ::  dx_c

REAL(KIND = idp), DIMENSION(0:ny),  INTENT(OUT)         ::  y_e
REAL(KIND = idp), DIMENSION(1:ny),  INTENT(OUT)         ::  y_c
REAL(KIND = idp), DIMENSION(1:ny),  INTENT(OUT)         ::  dy_c

REAL(KIND = idp), DIMENSION(0:nz),  INTENT(OUT)         ::  z_e
REAL(KIND = idp), DIMENSION(1:nz),  INTENT(OUT)         ::  z_c
REAL(KIND = idp), DIMENSION(1:nz),  INTENT(OUT)         ::  dz_c

REAL(KIND = idp),   INTENT(IN)                          ::  Zoom

! Create Geometric Radial Mesh !
Call Create_Geometric_1D_Mesh( Inner_Radius, Outer_Radius, Zoom, nx, x_e, x_c, dx_c )

! Create Uniform Theta Mesh !
CALL Create_Uniform_1D_Mesh( 0.0_idp, pi, ny, y_e, y_c, dy_c )

! Create Uniform Phi Mesh !
CALL Create_Uniform_1D_Mesh( 0.0_idp, 2.0_idp*pi, nz, z_e, z_c, dz_c )



END SUBROUTINE Create_Geometric_3D_Mesh












!+203+##############################################################################!
!                                                                                   !
!           Create_Split_3D_Mesh                                                    !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Creates a uniform angular mesh with a split radial mesh.  The nc elements are   !
!   spaced uniformly within Core_Radius while the (nx-nc) elements outside          !
!   Core_Radius are spaced logarithmicly.                                           !
!   Assumes y_range = (0:pi) and z_range = (0:2pi).                                 !
!                                                                                   !
!###################################################################################!
SUBROUTINE Create_Split_3D_Mesh(    Inner_Radius, Core_Radius, Outer_Radius,        &
                                    nx, nc, ny, nz,                                 &
                                    x_e, x_c, dx_c,                                 &
                                    y_e, y_c, dy_c,                                 &
                                    z_e, z_c, dz_c                                  )

REAL(KIND = idp),   INTENT(IN)                          ::  Inner_Radius,           &
                                                            Core_Radius,            &
                                                            Outer_Radius

INTEGER,            INTENT(IN)                          ::  nx, nc, ny, nz

REAL(KIND = idp), DIMENSION(0:nx),  INTENT(OUT)         ::  x_e
REAL(KIND = idp), DIMENSION(1:nx),  INTENT(OUT)         ::  x_c
REAL(KIND = idp), DIMENSION(1:nx),  INTENT(OUT)         ::  dx_c

REAL(KIND = idp), DIMENSION(0:ny),  INTENT(OUT)         ::  y_e
REAL(KIND = idp), DIMENSION(1:ny),  INTENT(OUT)         ::  y_c
REAL(KIND = idp), DIMENSION(1:ny),  INTENT(OUT)         ::  dy_c

REAL(KIND = idp), DIMENSION(0:nz),  INTENT(OUT)         ::  z_e
REAL(KIND = idp), DIMENSION(1:nz),  INTENT(OUT)         ::  z_c
REAL(KIND = idp), DIMENSION(1:nz),  INTENT(OUT)         ::  dz_c



! Create Split Radial Mesh !
CALL Create_Split_1D_Mesh( Inner_Radius, Core_Radius, Outer_Radius, nx, nc, x_e, x_c, dx_c)

! Create Theta Mesh !
CALL Create_Uniform_1D_Mesh( 0.0_idp, pi, ny, y_e, y_c, dy_c )

! Create Phi Mesh !
CALL Create_Uniform_1D_Mesh( 0.0_idp, 2.0_idp*pi, nz, z_e, z_c, dz_c )


END SUBROUTINE Create_Split_3D_Mesh












!+301+##############################################################################!
!                                                                                   !
!           Create_Uniform_1D_Mesh                                                  !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Creates a uniform radial mesh.                                                  !
!                                                                                   !
!###################################################################################!
SUBROUTINE Create_Uniform_1D_Mesh(  Inner_Edge, Outer_Edge,                         &
                                    nx, x_e, x_c, dx_c                              )


REAL(KIND = idp),   INTENT(IN)                          ::  Inner_Edge,           &
                                                            Outer_Edge

INTEGER,            INTENT(IN)                          ::  nx

REAL(KIND = idp), DIMENSION(0:nx),  INTENT(OUT)         ::  x_e
REAL(KIND = idp), DIMENSION(1:nx),  INTENT(OUT)         ::  x_c
REAL(KIND = idp), DIMENSION(1:nx),  INTENT(OUT)         ::  dx_c

INTEGER                                                 ::  i

dx_c(1:nx) = (Outer_Edge - Inner_Edge)/REAL(nx, KIND = idp)

x_e(0) = Inner_Edge
DO i = 1,nx

    x_e(i) = x_e(i-1) + dx_c(i)
    x_c(i) = x_e(i-1) + dx_c(i)/2.0_idp

END DO


END SUBROUTINE Create_Uniform_1D_Mesh




!+302+##############################################################################!
!                                                                                   !
!           Create_Logarithmic_1D_Mesh                                              !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Creates a uniform radial mesh.                                                  !
!                                                                                   !
!###################################################################################!
SUBROUTINE Create_Logarithmic_1D_Mesh(  Inner_Edge, Outer_Edge,                     &
                                        nx, x_e, x_c, dx_c                          )

REAL(KIND = idp),   INTENT(IN)                          ::  Inner_Edge,           &
                                                            Outer_Edge

INTEGER,            INTENT(IN)                          ::  nx

REAL(KIND = idp), DIMENSION(0:nx),  INTENT(OUT)         ::  x_e
REAL(KIND = idp), DIMENSION(1:nx),  INTENT(OUT)         ::  x_c
REAL(KIND = idp), DIMENSION(1:nx),  INTENT(OUT)         ::  dx_c

INTEGER                                                 ::  i
REAL(KIND = idp)                                        ::  tmp_Inner_Edge
REAL(KIND = idp)                                        ::  Scale_Factor
REAL(KIND = idp)                                        ::  Scale_Factor_Power




IF ( Inner_Edge == 0.0_idp ) THEN

    tmp_Inner_Edge = 1.0_idp
!    PRINT*,"Logarithmic mesh requires non-zero inner edge.  Shifted to 1."

ELSE

    tmp_Inner_Edge = Inner_Edge

END IF


Scale_Factor = (Outer_Edge/tmp_Inner_Edge)**(1.0_idp/REAL(nx,KIND = idp))


Scale_Factor_Power = 1.0_idp
x_e(0) = Inner_Edge ! * Scale_Factor_Power

DO i = 1,nx

    Scale_Factor_Power = Scale_Factor_Power*Scale_Factor
    x_e(i) = tmp_Inner_Edge*Scale_Factor_Power

    dx_c(i) = x_e(i) - x_e(i-1)
    x_c(i) = x_e(i-1) + dx_c(i)/2.0_idp

END DO


END SUBROUTINE Create_Logarithmic_1D_Mesh


!+303+##############################################################################!
!                                                                                   !
!           Create_Geometric_1D_Mesh                                              !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!                                                                                   !
!###################################################################################!
SUBROUTINE Create_Geometric_1D_Mesh(  Inner_Edge, Outer_Edge, Zoom,                   &
                                        nx, x_e, x_c, dx_c                          )

REAL(KIND = idp),   INTENT(IN)                          ::  Inner_Edge,           &
                                                            Outer_Edge
REAL(KIND = idp),   INTENT(IN)                          ::  Zoom

INTEGER,            INTENT(IN)                          ::  nx

REAL(KIND = idp), DIMENSION(0:nx),  INTENT(OUT)         ::  x_e
REAL(KIND = idp), DIMENSION(1:nx),  INTENT(OUT)         ::  x_c
REAL(KIND = idp), DIMENSION(1:nx),  INTENT(OUT)         ::  dx_c

INTEGER                                                 ::  i



dx_c(1) = (Outer_Edge - Inner_Edge)*(Zoom - 1.0_idp)/(Zoom**nx - 1.0_idp)
x_e(0)  = Inner_Edge
x_c(1)  = Inner_Edge + 0.5_idp*dx_c(1)
x_e(1)  = x_e(0) + dx_c(1)
DO i = 2,nx
    dx_c(i) = dx_c(i-1)*Zoom
    x_e(i)  = x_e(i-1) + dx_c(i)
    x_c(i)  = 0.5_idp*(x_e(i) + x_e(i-1))
END DO




END SUBROUTINE Create_Geometric_1D_Mesh





!+304+##############################################################################!
!                                                                                   !
!           Create_Split_1D_Mesh                                                    !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Creates a uniform radial mesh.                                                  !
!                                                                                   !
!###################################################################################!
SUBROUTINE Create_Split_1D_Mesh(    Inner_Edge, Core_Edge, Outer_Edge,              &
                                    nx, nc, x_e, x_c, dx_c                          )

REAL(KIND = idp),   INTENT(IN)                          ::  Inner_Edge,           &
                                                            Core_Edge,            &
                                                            Outer_Edge

INTEGER,            INTENT(IN)                          ::  nx, nc

REAL(KIND = idp), DIMENSION(0:nx),  INTENT(OUT)         ::  x_e
REAL(KIND = idp), DIMENSION(1:nx),  INTENT(OUT)         ::  x_c
REAL(KIND = idp), DIMENSION(1:nx),  INTENT(OUT)         ::  dx_c


INTEGER                                                 ::  i
REAL(KIND = idp)                                        ::  tmp_Inner_Edge
REAL(KIND = idp)                                        ::  Scale_Factor
REAL(KIND = idp)                                        ::  Scale_Factor_Power





! Create Uniform Radial Core !
dx_c(1:nc) = (Core_Edge-Inner_Edge)/REAL(nc,KIND = idp)

x_e(0) = Inner_Edge
DO i = 1,nc

    x_e(i) = x_e(i-1) + dx_c(i)
    x_c(i) = x_e(i-1) + dx_c(i)/2.0_idp

END DO


! Create Logarithmic Outer Shell !
Scale_Factor = (Outer_Edge/Core_Edge)**(1.0_idp/(nx-nc))

Scale_Factor_Power = 1.0_idp
DO i = nc+1,nx

    Scale_Factor_Power = Scale_Factor_Power*Scale_Factor
    x_e(i) = Core_Edge*Scale_Factor_Power

    dx_c(i) = x_e(i) - x_e(i-1)
    x_c(i) = x_e(i-1) + dx_c(i)/2.0_idp

END DO



END SUBROUTINE Create_Split_1D_Mesh







!+701+##################################################################################!
!                                                                                       !
!       GENERATE_DEFINED_MESH - Generate the values for the mesh sent in using          !
!                                predefined values for the width of each element.       !
!                                                                                       !
!---------------------------------------------------------------------------------------!
!                                                                                       !
!   Input:  Mesh_Start - Single Real number defining inner boundary location.           !
!                                                                                       !
!           Number_of_Elements - Single Integer value defining the number of elements   !
!                                   in the mesh.                                        !
!                                                                                       !
!           Element_Width_Vector - Real valued vector of length(1:Number_of_Elements)   !
!                                       containing Real numbers defining the width of   !
!                                        each element.                                  !
!                                                                                       !
!---------------------------------------------------------------------------------------!
!                                                                                       !
!   Output: Mesh - Real Vector,length(0:Number_of_Elements) that on output contains     !
!                       values describing the element edge locations.                   !
!                                                                                       !
!#######################################################################################!
SUBROUTINE Generate_Defined_Mesh(Number_of_Elements, Mesh_Start, Element_Width_Vector, Mesh)


INTEGER, INTENT(IN)                                                 ::  Number_of_Elements
REAL(KIND = idp), INTENT(IN)                                        ::  Mesh_Start
REAL(KIND = idp), DIMENSION(1:Number_of_Elements), INTENT(IN)       ::  Element_Width_Vector

REAL(KIND = idp), DIMENSION(0:Number_of_Elements), INTENT(OUT)      ::  Mesh

INTEGER                                                             ::  i


mesh(0) = Mesh_Start
DO i = 1,Number_of_Elements


    mesh(i) = mesh(i-1) + Element_Width_Vector(i)


END DO


END SUBROUTINE Generate_Defined_Mesh






!+802+##################################################################################!
!                                                                                       !
!       GENERATE_DEFINED_COARSE_MESH                                                    !
!                                                                                       !
!#######################################################################################!
SUBROUTINE Generate_Defined_Coarse_Mesh(Input_Number_of_Elements,           &
                                        Output_Number_of_Elements,          &
                                        Coarsen_Factor, Mesh_Start,         &
                                        Dim,                                &
                                        Mesh, dMesh )


INTEGER, INTENT(IN)                                                         ::  Input_Number_of_Elements
INTEGER, INTENT(IN)                                                         ::  Output_Number_of_Elements
INTEGER, INTENT(INOUT)                                                      ::  Coarsen_Factor
REAL(KIND = idp), INTENT(IN)                                                ::  Mesh_Start
INTEGER, INTENT(IN)                                                         ::  Dim

REAL(KIND = idp), DIMENSION(0:Output_Number_of_Elements), INTENT(INOUT)     ::  Mesh
REAL(KIND = idp), DIMENSION(1:Output_Number_of_Elements), INTENT(INOUT)     ::  dMesh

INTEGER                                                                     ::  i
REAL(KIND = idp)                                                            ::  tmp


IF ( Coarsen_Factor .NE. 1 ) THEN
    PRINT*,"!*  Warning in Poseidon *!"
    PRINT*,"  In Generate_Defined_Coarse_Mesh():    "
    PRINT*,"- Grid coarsening unavailable at this time. "
    PRINT*,"- Setting Coarsen_Factor to 1. "
    PRINT*,"!* End of Warning *!"
    Coarsen_Factor = 1
END IF




IF ( Output_Number_of_Elements*Coarsen_Factor == Input_Number_of_Elements ) THEN

    IF ( locs_Set(Dim) .EQV. .FALSE. ) THEN
        mesh(0) = Mesh_Start
        DO i = 1,Output_Number_of_Elements
            mesh(i) = mesh(i-1) + dMesh(i)
        END DO
    END IF

    IF ( dlocs_Set(Dim) .EQV. .FALSE. ) THEN
        
        DO i = 1,Output_Number_of_Elements
            dmesh(i) = Mesh(i)-Mesh(i-1)
        END DO

    END IF

ELSE

    PRINT*,"!* ERROR IN POSEIDON *!"
    PRINT*,"Error in Generate_Defined_Coarse_Mesh(): "
    PRINT*,"- Grid Coarsening Improperly Set Up. "
    PRINT*,"Output_Number_of_Elements*Coarsen_Factor != Input_Number_of_Elements"
    PRINT*,Output_Number_of_Elements,"*",Coarsen_Factor," = ",Output_Number_of_Elements*Coarsen_Factor,"!=",Input_Number_of_Elements
    PRINT*,"!* POSEIDON STOPPING PROGRAM  *!"
    STOP

END IF


END SUBROUTINE Generate_Defined_Coarse_Mesh









!+802+##################################################################################!
!                                                                                       !
!       Define_Refined_Mesh                                                             !
!                                                                                       !
!#######################################################################################!
SUBROUTINE Define_Refined_Mesh( l, n_l,                 &
                                radii,                  &
                                dx_array,               &
                                n_array                 )

INTEGER,                        INTENT(IN)          ::  l
INTEGER,                        INTENT(IN)          ::  n_l
REAL(idp), DIMENSION(1:l+1),    INTENT(INOUT)       ::  radii

REAL(idp), DIMENSION(1:l),      INTENT(OUT)         ::  dx_Array
INTEGER,   DIMENSION(1:l),      INTENT(OUT)         ::  n_Array

INTEGER                                             ::  i, j
REAL(idp)                                           ::  a, b
REAL(idp), DIMENSION(1:l)                           ::  dr_Array


PRINT*,l,n_l

dr_Array(l) = radii(l+1) - radii(l)
n_Array(l) = n_l
dx_Array(l) = dr_Array(l) / n_array(l)


DO  i = l,2,-1
    dr_Array(i-1) = radii(i) - radii(i-1)

    print*,"dr",i,dr_Array(i),dr_array(i-1),dr_array(i-1)/dr_Array(i)
    a = 2.0_idp * n_array(i) * dr_Array(i-1) / dr_Array(i)
    b = 3.0_idp * n_array(i) * dr_Array(i-1) / dr_Array(i)

    print*,a,b, (b-a)/2.0_idp+a
    n_Array(i-1)  = floor( ( b - a ) / 2.0_idp + a )
    PRINT*,"n",n_Array(i-1)
    dx_Array(i-1) = dr_Array(i-1) / n_Array(i-1)
END DO


END SUBROUTINE Define_Refined_Mesh




END MODULE Functions_Mesh
