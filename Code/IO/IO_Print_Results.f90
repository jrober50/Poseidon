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
                    L_Limit,            &
                    Poisson_Mode

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,     &
                    rlocs

USE Variables_Derived, &
            ONLY :  Num_R_Nodes

USE Variables_Poisson, &
            ONLY :  Coefficient_Vector



USE Variables_Functions, &
            ONLY :  Calc_3D_Values_At_Location

USE Variables_Mesh, &
            ONLY :  R_Inner,        &
                    R_Outer


USE Functions_Mesh, &
            ONLY :  Create_Logarithmic_1D_Mesh,     &
                    Create_Uniform_1D_Mesh

USE Functions_Math, &
            ONLY :  Lagrange_Poly,      &
                    Spherical_Harmonic
USE Functions_Mapping, &
            ONLY :  Map_To_X_Space

USE Functions_Quadrature, &
            ONLY : Initialize_LGL_Quadrature

USE Timer_Routines_Module, &
            ONLY :  TimerStart,     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Core_PrintResults

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

110 FORMAT (11X,A1,24X,A3,19X,A8,15X,A11,14X,A11,14X,A11)
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


WRITE(*,'(A,F4.2,A,F4.2,A)')"Results taken along ray, theta = ",theta/pi," Pi Radians, Phi = ",phi/pi," Pi Radians"



IF ( Poisson_Mode ) THEN


    WRITE(*,210)"r","Code Results","Lapse Function"
    WRITE(*,211)"----------------","----------------","----------------"

    DO i = 0,NUM_SAMPLES

        Potential =  Calculate_Potential_At_Location( x_e(i) )

        Return_Psi = 1.0_idp    &
                    - 0.5_idp*Potential/C_Square

        WRITE(*,212)x_e(i)/Centimeter,potential,Return_Psi


    END DO

ELSE

    WRITE(*,110)"r","Psi","AlphaPsi","Beta1 Value","Beta2 Value","Beta3 Value"


    DO i = 0,Num_Samples

        r = x_e(i)



        CALL Calc_3D_Values_At_Location( r, theta, phi,                              &
                                        Return_Psi, Return_AlphaPsi,                &
                                        Return_Beta1, Return_Beta2, Return_Beta3    )


        Potential = 2.0_idp*C_Square*(1.0_idp - Return_Psi)

!        WRITE(*,212)x_e(i)/Centimeter,potential,Return_Psi


        WRITE(*,111) r/Centimeter,              &
                     Return_Psi,                &
                     Return_AlphaPsi,           &
                     Return_Beta1/Shift_Units,  &
                     Return_Beta2,              &
                     Return_Beta3

    END DO

END IF ! Poisson_Mode


CALL TimerStop( Timer_Core_PrintResults )


END SUBROUTINE Print_Results









!+300+##################################################################!
!                                                                       !
!       CALC_POTENTIAL - Use coefficents to calculate the potential.    !
!                                                                       !
!#######################################################################!
FUNCTION Calculate_Potential_At_Location( r_input, theta_input, phi_input)

REAL(KIND = idp)                                               ::  Calculate_Potential_At_Location


REAL(KIND = idp),               INTENT(IN)                      ::  r_input
REAL(KIND = idp),               INTENT(IN), OPTIONAL            ::  theta_input, phi_input




INTEGER                                 ::  l, m, re, rd
REAL(KIND = idp)                        ::  r, theta, phi
REAL(KIND = idp)                        ::  r_tmp
REAL(KIND = idp), DIMENSION(0:DEGREE)   ::  LagP
REAL(KIND = idp), DIMENSION(0:DEGREE)   ::  xlocP, weightP


COMPLEX(KIND = idp)                     ::  potential



r = r_input

if ( PRESENT(theta_input) ) THEN

    theta = theta_input

ELSE

    theta = 0.0_idp

END IF


IF (PRESENT(phi_input) ) THEN

    phi = phi_input

ELSE

    phi = 0.0_idp

END IF






potential = 0.0_idp

IF (r == rlocs(0)) THEN

    DO l = 0,L_LIMIT

        DO m = -l,l

            potential = potential + Coefficient_Vector(0,m,l)*Spherical_Harmonic(l,m,theta,phi)

        END DO

    END DO

ELSE

    CALL Initialize_LGL_Quadrature(DEGREE,xlocP,weightP)


    DO re = 0, NUM_R_ELEMENTS-1

        IF (r > rlocs(re) .AND. r <= rlocs(re+1)) THEN

            r_tmp = Map_To_X_Space(rlocs(re),rlocs(re+1),r)

            LagP = Lagrange_Poly(r_tmp,DEGREE,xlocP)


            DO l = 0,L_LIMIT

                DO m = -l,l

                    DO rd = 0,DEGREE


                        potential = potential + Coefficient_Vector(re*DEGREE + rd, m, l)*Spherical_Harmonic(l,m,theta,phi)*LagP(rd)


                    END DO

                END DO

            END DO


        END IF

    END DO



    IF (    r > rlocs(NUM_R_ELEMENTS)   ) THEN



        r_tmp = Map_To_X_Space(rlocs(NUM_R_ELEMENTS-1),rlocs(NUM_R_ELEMENTS),r)

        LagP = Lagrange_Poly(r_tmp,DEGREE,xlocP)


        DO l = 0,L_LIMIT

            DO m = -l,l

                DO rd = 0,DEGREE


                    potential = potential + Coefficient_Vector(NUM_R_NODES-1, m, l)    &
                                            * Spherical_Harmonic(l,m,theta,phi)                      &
                                            * LagP(rd)


                END DO

            END DO

        END DO



    END IF

END IF



Calculate_Potential_At_Location = REAL( potential )




END FUNCTION Calculate_Potential_At_Location


END MODULE IO_Print_Results
