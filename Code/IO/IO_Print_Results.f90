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

USE Poseidon_Units_Module, &
            ONLY :  C_Square,           &
                    Centimeter,         &
                    Shift_Units

USE Poseidon_Parameters, &
            ONLY :  DEGREE,             &
                    L_Limit

USE Poseidon_IO_Parameters, &
            ONLY :  CFA_Var_Names

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,     &
                    rlocs

USE Variables_Derived, &
            ONLY :  Num_R_Nodes

USE Variables_MPI, &
            ONLY :  myID_Poseidon,      &
                    MasterID_Poseidon

USE Variables_Functions, &
            ONLY :  Calc_3D_Values_At_Location

USE Poseidon_Return_Routines_Module, &
            ONLY :  Calc_Var_At_Location

USE Variables_Mesh, &
            ONLY :  R_Inner,        &
                    R_Outer

USE Variables_Vectors, &
            ONLY :  cVA_Coeff_Vector

USE Functions_Mesh, &
            ONLY :  Create_Logarithmic_1D_Mesh,     &
                    Create_Uniform_1D_Mesh

USE Functions_Math, &
            ONLY :  Lagrange_Poly,      &
                    Spherical_Harmonic

USE Maps_X_Space, &
            ONLY :  Map_To_X_Space

USE Maps_Domain, &
            ONLY :  Map_To_LM

USE Functions_Quadrature, &
            ONLY : Initialize_LGL_Quadrature

USE Timer_Routines_Module, &
            ONLY :  TimerStart,     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Core_PrintResults

USE Flags_Core_Module, &
            ONLY :  lPF_Core_Flags,         &
                    iPF_Core_Newtonian_Mode

IMPLICIT NONE


INTERFACE Print_Single_Var_Results
    MODULE PROCEDURE Print_Single_Var_Results_A
    MODULE PROCEDURE Print_Single_Var_Results_B
END INTERFACE Print_Single_Var_Results

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
REAL(KIND = idp)                                        ::  r, theta, phi
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  x_e
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  x_c, dx_c

REAL(KIND = idp)                                        ::  Return_Psi
REAL(KIND = idp)                                        ::  Return_AlphaPsi
REAL(KIND = idp)                                        ::  Return_Beta1
REAL(KIND = idp)                                        ::  Return_Beta2
REAL(KIND = idp)                                        ::  Return_Beta3

REAL(idp)                                               ::  Potential

INTEGER                                                 ::  Num_Samples = 20

110 FORMAT (11X,A,17X,A,10X,A,12X,A,13X,A,14X,A)
111 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)

210 FORMAT (11X,A1,16X,A12,16X,A12)
211 FORMAT (4X,A16,7X,A16,7X,A16)
212 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15)



CALL TimerStart( Timer_Core_PrintResults )


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





IF ( lPF_Core_Flags(iPF_Core_Newtonian_Mode) ) THEN


    CALL Print_Single_Var_Results_A( 1 )

ELSE

    WRITE(*,'(A,F4.2,A,F4.2,A)')"Results taken along ray, theta = ",theta/pi," Pi Radians, Phi = ",phi/pi," Pi Radians"


    WRITE(*,110)"r","Conformal Factor","Lapse Function","Radial Shift","Theta Shift","Phi Shift"


    DO i = 0,Num_Samples

        r = x_e(i)



        CALL Calc_3D_Values_At_Location( r, theta, phi,                              &
                                        Return_Psi, Return_AlphaPsi,                &
                                        Return_Beta1, Return_Beta2, Return_Beta3    )


        Potential = 2.0_idp*C_Square*(1.0_idp - Return_Psi)

!        WRITE(*,212)x_e(i)/Centimeter,potential,Return_Psi


        WRITE(*,111) r/Centimeter,                  &
                     Return_Psi,                    &
                     Return_AlphaPsi/Return_Psi,    &
                     Return_Beta1/Shift_Units,      &
                     Return_Beta2,                  &
                     Return_Beta3

    END DO

END IF ! lPF_Core_Flags(iPF_Core_Newtonian_Mode)


CALL TimerStop( Timer_Core_PrintResults )


END SUBROUTINE Print_Results

















!+302+##########################################################################!
!                                                                               !
!                   Print_Results                                               !
!                                                                               !
!###############################################################################!
SUBROUTINE Print_Single_Var_Results_A( iU )

INTEGER, INTENT(IN)                                     ::  iU

INTEGER                                                 ::  i
REAL(KIND = idp)                                        ::  r, theta, phi
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  x_e
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  x_c, dx_c

REAL(idp)                                               ::  Value

INTEGER                                                 ::  Num_Samples = 20


110 FORMAT (11X,A1,17X,A)
111 FORMAT (ES22.15,3X,ES22.15,3X)



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
WRITE(*,110)"r",CFA_Var_Names(iU)


DO i = 0,Num_Samples

    r = x_e(i)

    Value = Calc_Var_At_Location(r, theta, phi, iU)

    WRITE(*,111) r/Centimeter, Value

END DO




END SUBROUTINE Print_Single_Var_Results_A





!+302+##########################################################################!
!                                                                               !
!                   Print_Results                                               !
!                                                                               !
!###############################################################################!
SUBROUTINE Print_Single_Var_Results_B( iU, iVB)

INTEGER, INTENT(IN)                                     ::  iU
INTEGER, INTENT(IN)                                     ::  iVB

INTEGER                                                 ::  i
REAL(KIND = idp)                                        ::  r, theta, phi
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  x_e
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  x_c, dx_c


REAL(idp)                                               ::  Value

INTEGER                                                 ::  Num_Samples = 20



110 FORMAT (11X,A1,17X,A)
111 FORMAT (ES22.15,3X,ES22.15,3X)



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
WRITE(*,110)"r",CFA_Var_Names(iU)


DO i = 0,Num_Samples

    r = x_e(i)

    Value = Calc_Var_At_Location(r, theta, phi, iU, iVB)

    WRITE(*,111) r/Centimeter, Value

END DO


END SUBROUTINE Print_Single_Var_Results_B



END MODULE IO_Print_Results
