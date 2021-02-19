   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE IO_Print_Results                                                             !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
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
                    ONLY :  pi

USE Units_Module, &
                    ONLY :  Centimeter,     &
                            Shift_Units

USE Variables_Functions, &
                    ONLY :  Calc_3D_Values_At_Location

USE Variables_Mesh, &
                    ONLY :  R_Inner,        &
                            R_Outer


USE Functions_Mesh, &
                    ONLY :  Create_Logarithmic_1D_Mesh







IMPLICIT NONE




!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS




!+302+##########################################################################!
!                                                                               !
!                   Print_Results                                               !
!                                                                               !
!###############################################################################!
SUBROUTINE Print_Results()

INTEGER                                                 ::  i
REAL(KIND = idp)                                        ::  r, theta, phi, deltar
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  x_e
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  x_c, dx_c

REAL(KIND = idp)                                        ::  Return_Psi
REAL(KIND = idp)                                        ::  Return_AlphaPsi
REAL(KIND = idp)                                        ::  Return_Beta1
REAL(KIND = idp)                                        ::  Return_Beta2
REAL(KIND = idp)                                        ::  Return_Beta3


INTEGER                                                 ::  Num_Samples = 20

110 FORMAT (11X,A1,24X,A3,19X,A8,15X,A11,14X,A11,14X,A11)
111 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)

ALLOCATE( x_e(0:Num_Samples) )
ALLOCATE( x_c(1:Num_Samples) )
ALLOCATE( dx_c(1:Num_Samples) )


Call Create_Logarithmic_1D_Mesh( R_Inner,           &
                                 R_Outer,           &
                                 Num_Samples,       &
                                 x_e, x_c, dx_c     )

theta = pi/3.0_idp
phi = pi/2.0_idp


WRITE(*,'(A,F4.2,A,F4.2,A)')"Results taken along ray, theta = ",theta/pi," Pi Radians, Phi = ",phi/pi," Pi Radians"
WRITE(*,110)"r","Psi","AlphaPsi","Beta1 Value","Beta2 Value","Beta3 Value"


DO i = 0,Num_Samples

    r = x_e(i)



    CALL Calc_3D_Values_At_Location( r, theta, phi,                              &
                                    Return_Psi, Return_AlphaPsi,                &
                                    Return_Beta1, Return_Beta2, Return_Beta3    )

    WRITE(*,111) r/Centimeter,              &
                 Return_Psi,                &
                 Return_AlphaPsi,           &
                 Return_Beta1/Shift_Units,  &
                 Return_Beta2,              &
                 Return_Beta3

END DO



END SUBROUTINE Print_Results







END MODULE IO_Print_Results
