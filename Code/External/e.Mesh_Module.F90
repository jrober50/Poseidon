   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Mesh_Module                                                                  !##!
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

USE Poseidon_Constants_Module, &
                ONLY :  idp, pi

USE Driver_Parameters, &
                ONLY :  Driver_Zoom


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
SUBROUTINE Create_3D_Mesh(  Mesh_Type,                                              &
                            Inner_Radius, Core_Radius, Outer_Radius,                &
                            nx, nc, ny, nz,                                         &
                            x_e, x_c, dx_c,                                         &
                            y_e, y_c, dy_c,                                         &
                            z_e, z_c, dz_c                                          )

INTEGER,            INTENT(IN)                          ::  Mesh_Type

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

INTEGER                                                 ::  i
REAL(KIND = idp)                                        ::  tmp_Inner_Radius
REAL(KIND = idp)                                        ::  Scale_Factor
REAL(KIND = idp)                                        ::  Scale_Factor_Power




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

    CALL Create_Split_3D_Mesh( Inner_Radius, Core_Radius, Outer_Radius,      &
                               nx, nc, ny, nz,                               &
                               x_e, x_c, dx_c,                               &
                               y_e, y_c, dy_c,                               &
                               z_e, z_c, dz_c                                )

ELSE IF ( Mesh_Type == 4) THEN

    CALL Create_Geometric_3D_Mesh( Inner_Radius, Outer_Radius,         &
                                     nx, ny, nz,                       &
                                     x_e, x_c, dx_c,                   &
                                     y_e, y_c, dy_c,                   &
                                     z_e, z_c, dz_c                    )


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



! Create Geometric Radial Mesh !
Call Create_Geometric_1D_Mesh( Inner_Radius, Outer_Radius, Driver_Zoom, nx, x_e, x_c, dx_c )

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
x_e(0) = tmp_Inner_Edge ! * Scale_Factor_Power

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
    dx_c(i) = dx_c(i-1)*Zoom                ! Width
    x_e(i)  = x_e(i-1) + dx_c(i)            ! Right Edge
!    x_c(i)  = x_e(0) + SUM(dx_c(1:i-1))+dx_c(i)/2.0_idp
    x_c(i)  = 0.5_idp*(x_e(i) + x_e(i-1))   ! Center
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








END MODULE Mesh_Module
