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

USE Parameters_Variable_Indices, &
            ONLY :  iVB_X,                      &
                    iVB_S,                      &
                    iU_CF,                      &
                    iU_LF,                      &
                    iU_S1,                      &
                    iU_S2,                      &
                    iU_S3,                      &
                    iU_X1,                      &
                    iU_X2,                      &
                    iU_X3

USE Poseidon_Units_Module, &
            ONLY :  C_Square,           &
                    Centimeter,         &
                    Shift_Units

USE Poseidon_IO_Parameters, &
            ONLY :  CFA_Var_Names
            
USE Variables_MPI, &
            ONLY :  myID_Poseidon,      &
                    MasterID_Poseidon,  &
                    nPROCS_Poseidon

USE Variables_Functions, &
            ONLY :  Calc_3D_Values_At_Location

USE Poseidon_Return_Routines_Module, &
            ONLY :  Calc_Var_At_Location

USE Variables_Mesh, &
            ONLY :  R_Inner,        &
                    R_Outer

USE Functions_Mesh, &
            ONLY :  Create_Logarithmic_1D_Mesh,     &
                    Create_Uniform_1D_Mesh

USE Flags_Core_Module, &
            ONLY :  iPF_Core_Flags,         &
                    iPF_Core_Method_Mode,   &
                    iPF_Core_Method_Newtonian

USE MPI
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

INTEGER                                                 ::  ierr, id

INTEGER                                                 ::  Num_Samples = 20

110 FORMAT (11X,A,17X,A,10X,A,12X,A,13X,A,14X,A)
111 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)


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



DO id = 0,nPROCS_Poseidon-1

IF ( id == myID_Poseidon ) THEN

IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_Newtonian ) THEN


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

END IF ! iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_Newtonian

END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
END DO



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






!+302+##########################################################################!
!                                                                               !
!                   Print_Results                                               !
!                                                                               !
!###############################################################################!
SUBROUTINE Print_Vector_Var_Results( iVB)

INTEGER, INTENT(IN)                                     ::  iVB

INTEGER                                                 ::  i
REAL(KIND = idp)                                        ::  r, theta, phi
REAL(KIND = idp),   DIMENSION(:),   ALLOCATABLE         ::  x_e
REAL(KIND = idp),   DIMENSION(:),   ALLOCATABLE         ::  x_c, dx_c

INTEGER,            DIMENSION(3)                        ::  iU
REAL(idp),          DIMENSION(3)                        ::  Values

INTEGER                                                 ::  Num_Samples = 20

IF ( iVB == iVB_X ) THEN
    iU = [iU_X1, iU_X2, iU_X3]
ElSE if (iVB == iVB_S) THEN
    iU = [iU_S1, iU_S2, iU_S3]
END IF

142 FORMAT (11X,A1,16X,A,8X,A,8X,A)
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
WRITE(*,142)"r",                            &
            TRIM(CFA_Var_Names(iU(1))),     &
            TRIM(CFA_Var_Names(iU(2))),     &
            TRIM(CFA_Var_Names(iU(3)))


DO i = 0,Num_Samples

    r = x_e(i)

    Values(1) = Calc_Var_At_Location(r, theta, phi, iU(1), iVB)
    Values(2) = Calc_Var_At_Location(r, theta, phi, iU(2), iVB)
    Values(3) = Calc_Var_At_Location(r, theta, phi, iU(3), iVB)

    WRITE(*,111) r/Centimeter, Values

END DO


END SUBROUTINE Print_Vector_Var_Results



END MODULE IO_Print_Results
