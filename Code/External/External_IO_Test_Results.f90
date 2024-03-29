   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE External_IO_Test_Results_Module                                       !##!
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
            ONLY :  idp

USE Poseidon_Numbers_Module, &
            ONLY :  pi

USE Poseidon_Units_Module, &
            ONLY :  C_Square,       &
                    Centimeter,     &
                    Second,         &
                    GravPot_Units

USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,          &
                    iU_X1,          &
                    iVB_X

USE Variables_Mesh, &
            ONLY :  R_Inner,        &
                    R_Outer

USE External_HCT_Solution_Module, &
            ONLY :  HCT_Solution

USE External_MVL_Solution_Module, &
            ONLY :  MVL_Solution,       &
                    Set_MVL_Test_Params

USE External_UST_Solution_Module, &
            ONLY :  UST_Solution

USE External_MLS_Solution_Module, &
            ONLY :  MacLaurin_Potential_Sub

USE Poseidon_Return_Routines_Module, &
            ONLY :  Calc_Var_At_Location

USE Functions_Mesh, &
            ONLY :  Create_Logarithmic_1D_Mesh,     &
                    Create_Uniform_1D_Mesh

USE Variables_Functions, &
            ONLY :  Potential_Solution
            
IMPLICIT NONE


CONTAINS



 !+101+####################################################!
!                                                           !
!          	                                      	    !
!                                                           !
 !#########################################################!
SUBROUTINE Print_HCT_Error()


INTEGER                                                 ::  i
REAL(KIND = idp)                                        ::  r, theta, phi
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  x_e
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  x_c, dx_c


REAL(KIND = idp)                                        ::  Return_Psi

REAL(idp)                                               ::  HCT_Sol
REAL(idp)                                               ::  HCT_Error


INTEGER                                                 ::  Num_Samples = 20

110 FORMAT (11X,A,12X,A,8X,A,14X,A)
111 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)



ALLOCATE( x_e(0:Num_Samples) )
ALLOCATE( x_c(1:Num_Samples) )
ALLOCATE( dx_c(1:Num_Samples) )

IF ( R_Outer - R_Inner > 1000.0_idp ) THEN
    Call Create_Logarithmic_1D_Mesh( R_Inner,           &
                                     R_Outer,           &
                                     Num_Samples,       &
                                     x_e, x_c, dx_c     )
ELSE

    CALL Create_Uniform_1D_Mesh( R_Inner, R_Outer, Num_Samples, x_e, x_c, dx_c )

END IF

theta = 0.5_idp * pi
phi = 0.5_idp * pi



WRITE(*,'(A,F4.2,A,F4.2,A)')"Results taken along ray, theta = ",theta/pi," Pi Radians, Phi = ",phi/pi," Pi Radians"


WRITE(*,110)"r (cm)","Conformal Factor","Analytic Solution","Relative Error"


DO i = 0,Num_Samples

    r = x_e(i)

    Return_Psi = Calc_Var_At_Location(r,theta,phi,iU_CF)
    
    HCT_Sol = HCT_Solution(r)
    HCT_Error = abs(HCT_Sol - Return_Psi)/MAXVAL([abs(HCT_Sol),abs(Return_Psi)])

    WRITE(*,111) r/Centimeter,      &
                 Return_Psi,         &
                 HCT_Sol,           &
                 HCT_Error

END DO






END SUBROUTINE Print_HCT_Error






 !+101+####################################################!
!                                                           !
!                                                            !
!                                                           !
 !#########################################################!
SUBROUTINE Print_MVL_Error()


INTEGER                                                 ::  i
REAL(KIND = idp)                                        ::  r, theta, phi
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  x_e
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  x_c, dx_c


REAL(KIND = idp)                                        ::  Return_X1


REAL(idp)                                               ::  MVL_Sol
REAL(idp)                                               ::  MVL_Error


INTEGER                                                 ::  Num_Samples = 20

110 FORMAT (9X,A,17X,A,13X,A,14X,A)
111 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)

CALL Set_MVL_Test_Params( )

ALLOCATE( x_e(0:Num_Samples) )
ALLOCATE( x_c(1:Num_Samples) )
ALLOCATE( dx_c(1:Num_Samples) )

IF ( R_Outer - R_Inner > 1000.0_idp ) THEN
    Call Create_Logarithmic_1D_Mesh( R_Inner,           &
                                     R_Outer,           &
                                     Num_Samples,       &
                                     x_e, x_c, dx_c     )
ELSE

    CALL Create_Uniform_1D_Mesh( R_Inner, R_Outer, Num_Samples, x_e, x_c, dx_c )

END IF


theta = 0.5_idp * pi
phi = 0.5_idp * pi



WRITE(*,'(A,F4.2,A,F4.2,A)')"Results taken along ray, theta = ",theta/pi," Pi Radians, Phi = ",phi/pi," Pi Radians"


WRITE(*,110)"r (cm)","Radial X","Analytic Solution","Error"


DO i = 0,Num_Samples

    r = x_e(i)

    Return_X1 = Calc_Var_At_Location(r,theta,phi,iU_X1,iVB_X)
    
    MVL_Sol = MVL_Solution(r)
    MVL_Error = abs(MVL_Sol - Return_X1)/ABS(MVL_SoL)

    WRITE(*,111) r/Centimeter,      &
                 Return_X1,         &
                 MVL_Sol,           &
                 MVL_Error

END DO


END SUBROUTINE Print_MVL_Error






 !+101+####################################################!
!                                                           !
!                                                            !
!                                                           !
 !#########################################################!
SUBROUTINE Print_UST_Error()


INTEGER                                                 ::  i
REAL(KIND = idp)                                        ::  r, theta, phi
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  x_e
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  x_c, dx_c


REAL(KIND = idp)                                        ::  Potential

REAL(idp)                                               ::  UST_Sol
REAL(idp)                                               ::  UST_Error


INTEGER                                                 ::  Num_Samples = 20

110 FORMAT (8X,A,18X,A,12X,A,9X,A)
111 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)



ALLOCATE( x_e(0:Num_Samples) )
ALLOCATE( x_c(1:Num_Samples) )
ALLOCATE( dx_c(1:Num_Samples) )

IF ( R_Outer - R_Inner > 1000.0_idp ) THEN
    Call Create_Logarithmic_1D_Mesh( R_Inner,           &
                                     R_Outer,           &
                                     Num_Samples,       &
                                     x_e, x_c, dx_c     )
ELSE

    CALL Create_Uniform_1D_Mesh( R_Inner, R_Outer, Num_Samples, x_e, x_c, dx_c )

END IF


theta = 0.5_idp * pi
phi = 0.5_idp * pi



WRITE(*,'(A,F4.2,A,F4.2,A)')"Results taken along ray, theta = ",theta/pi," Pi Radians, Phi = ",phi/pi," Pi Radians"


WRITE(*,110)"r (cm)","Potential","Analytic Solution","Relative Error"


DO i = 0,Num_Samples

    r = x_e(i)

    Potential = Calc_Var_At_Location(r,theta,phi,iU_CF)
    
    UST_Sol   = UST_Solution(r)
    UST_Error = abs(UST_Sol - Potential)/MAXVAL([abs(UST_Sol),abs(Potential)])

    WRITE(*,111) r/Centimeter,      &
                 Potential,         &
                 UST_Sol,           &
                 UST_Error

END DO


END SUBROUTINE Print_UST_Error





 !+101+####################################################!
!                                                           !
!                                                            !
!                                                           !
 !#########################################################!
SUBROUTINE Print_MacLaurin_Error()


INTEGER                                                 ::  i
REAL(KIND = idp)                                        ::  r, theta, phi
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  x_e
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  x_c, dx_c

REAL(KIND = idp)                                        ::  Psi
REAL(KIND = idp)                                        ::  Potential

REAL(idp)                                               ::  MLS_Sol
REAL(idp)                                               ::  MLS_Error


INTEGER                                                 ::  Num_Samples = 20

110 FORMAT (8X,A,18X,A,12X,A,9X,A)
111 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)



ALLOCATE( x_e(0:Num_Samples) )
ALLOCATE( x_c(1:Num_Samples) )
ALLOCATE( dx_c(1:Num_Samples) )

IF ( R_Outer - R_Inner > 1000.0_idp ) THEN
    Call Create_Logarithmic_1D_Mesh( R_Inner,           &
                                     R_Outer,           &
                                     Num_Samples,       &
                                     x_e, x_c, dx_c     )
ELSE

    CALL Create_Uniform_1D_Mesh( R_Inner, R_Outer, Num_Samples, x_e, x_c, dx_c )

END IF


theta = 0.5_idp * pi
phi = 0.5_idp * pi



WRITE(*,'(A,F4.2,A,F4.2,A)')"Results taken along ray, theta = ",theta/pi," Pi Radians, Phi = ",phi/pi," Pi Radians"


WRITE(*,110)"r (cm)","Potential","Analytic Solution","Relative Error"


DO i = 0,Num_Samples

    r = x_e(i)

    Psi = Calc_Var_At_Location(r,theta,phi,iU_CF)
!    Potential = 2.0_idp*C_Square*(1.0_idp - Psi)
    Potential = Psi
    
    CALL MacLaurin_Potential_Sub(r, theta, phi, MLS_Sol)
    MLS_Sol = 1.0_idp - 0.5_idp*MLS_Sol/C_Square
    
    
    MLS_Error = abs(MLS_Sol - Potential)/MAXVAL([abs(MLS_Sol),abs(Potential)])

    WRITE(*,111) r/Centimeter,      &
                 Potential,         &
                 MLS_Sol,           &
                 MLS_Error

END DO


END SUBROUTINE Print_MacLaurin_Error




 !+101+####################################################!
!                                                           !
!                                                            !
!                                                           !
 !#########################################################!
SUBROUTINE Print_Yahil_Error()


INTEGER                                                 ::  i
REAL(KIND = idp)                                        ::  r, theta, phi
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  x_e
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  x_c, dx_c


REAL(KIND = idp)                                        ::  Potential

REAL(idp)                                               ::  Yahil_Sol
REAL(idp)                                               ::  Yahil_Error


INTEGER                                                 ::  Num_Samples = 20

110 FORMAT (8X,A,18X,A,12X,A,9X,A)
111 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)



ALLOCATE( x_e(0:Num_Samples) )
ALLOCATE( x_c(1:Num_Samples) )
ALLOCATE( dx_c(1:Num_Samples) )

IF ( R_Outer - R_Inner > 1000.0_idp ) THEN
    Call Create_Logarithmic_1D_Mesh( R_Inner,           &
                                     R_Outer,           &
                                     Num_Samples,       &
                                     x_e, x_c, dx_c     )
ELSE

    CALL Create_Uniform_1D_Mesh( R_Inner, R_Outer, Num_Samples, x_e, x_c, dx_c )

END IF


theta = 0.5_idp * pi
phi = 0.5_idp * pi



WRITE(*,'(A,F4.2,A,F4.2,A)')"Results taken along ray, theta = ",theta/pi," Pi Radians, Phi = ",phi/pi," Pi Radians"


WRITE(*,110)"r (cm)","Potential","Analytic Solution","Relative Error"


DO i = 0,Num_Samples

    r = x_e(i)

    Potential = Calc_Var_At_Location(r,theta,phi,iU_CF)
    
    Yahil_Sol = Potential_Solution(r,theta,phi)

    Yahil_Error = abs(Yahil_Sol - Potential)/MAXVAL([abs(Yahil_Sol),abs(Potential)])

    WRITE(*,111) r/Centimeter,      &
                 Potential,         &
                 Yahil_Sol,           &
                 Yahil_Error

END DO


END SUBROUTINE Print_Yahil_Error



END MODULE External_IO_Test_Results_Module
