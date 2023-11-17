   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE External_MLS_Profile_Module                                                  !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Initialize_MacLaurin_Sources                                        !##!
!##!                                                                                !##!
!##!    +201+   Test_Solution_MacLaurin                                             !##!
!##!                                                                                !##!
!##!    +301+   MacLaurin_Radius                                                    !##!
!##!    +302+   MacLaurin_Root_Finder                                               !##!
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

USE Poseidon_Units_Module, &
            ONLY :  Grav_Constant_G,        &
                    Gram,                  &
                    Centimeter,             &
                    C_Square

USE Poseidon_Numbers_Module, &
            ONLY :  pi

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Variables_MPI, &
            ONLY :  myID_Poseidon,          &
                    MasterID_Poseidon
                    
USE Variables_External, &
            ONLY :  MLS_SemiMinor,          &
                    MLS_SemiMajor,          &
                    MLS_Ecc,                &
                    MLS_SphereType,         &
                    MLS_SphereName,         &
                    MLS_Rho,                &
                    iMLS_Oblate,            &
                    iMLS_Prolate

USE Maps_X_Space, &
            ONLY :  Map_To_X_Space,         &
                    Map_From_X_Space

USE Maps_Quadrature,   &
            ONLY :  Quad_Map

USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,           &
                    iPF_IO_Print_Setup


IMPLICIT NONE



CONTAINS



 !+101+################################################################!
!                                                                       !
!          Initialize_MacLaurin_Sources                                 !
!                                                                       !
 !#####################################################################!
SUBROUTINE Initialize_MacLaurin_Sources( Rho_O, Spheroid_Type,              &
                                         SemiMajor, SemiMinor,              &
                                         Num_Quad,                          &
                                         R_Quad, T_Quad, P_Quad,            &
                                         Left_Limit, Right_Limit,           &
                                         Num_Elem,                          &
                                         R_Locs, T_Locs, P_Locs,            &
                                         Output_E, Output_S, Output_Si      )


REAL(idp),  INTENT(IN)                                                      :: Rho_O
CHARACTER(LEN=7), INTENT(IN)                                                :: Spheroid_Type

REAL(idp),  INTENT(IN)                                                      :: SemiMajor
REAL(idp),  INTENT(IN)                                                      :: SemiMinor

INTEGER,    INTENT(IN), DIMENSION( 1:3 )                                    :: Num_Quad
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Quad(1) )                          :: R_Quad
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Quad(2) )                          :: T_Quad
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Quad(3) )                          :: P_Quad
REAL(idp),  INTENT(IN)                                                      :: Left_Limit
REAL(idp),  INTENT(IN)                                                      :: Right_Limit

INTEGER,    INTENT(IN), DIMENSION( 1:3 )                                    :: Num_Elem
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Elem(1)+1 )                        :: R_Locs
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Elem(2)+1 )                        :: T_Locs
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Elem(3)+1 )                        :: P_Locs

REAL(idp),  INTENT(OUT),DIMENSION( 1:Num_Quad(1)*Num_Quad(2)*Num_Quad(3),   &
                                   0:Num_Elem(1)-1,                         &
                                   0:Num_Elem(2)-1,                         &
                                   0:Num_Elem(3)-1  )                       ::  Output_E,   &
                                                                                Output_S

REAL(idp), INTENT(OUT), DIMENSION( 1:Num_Quad(1)*Num_Quad(2)*Num_Quad(3),   &
                                   0:Num_Elem(1)-1,                         &
                                   0:Num_Elem(2)-1,                         &
                                   0:Num_Elem(3)-1,                         &
                                   1:3                  )                   ::  Output_Si



INTEGER                                                     ::  re, te, pe
INTEGER                                                     ::  rd, td, pd
INTEGER                                                     ::  Tot_Quad

REAL(idp),  DIMENSION( 1:Num_Quad(1) )                      ::  RQ_Loc
REAL(idp),  DIMENSION( 1:Num_Quad(2) )                      ::  TQ_Loc
REAL(idp),  DIMENSION( 1:Num_Quad(3) )                      ::  PQ_Loc

REAL(idp),  DIMENSION( 1:Num_Quad(1) )                      ::  rsqr

REAL(idp),  DIMENSION( 1:Num_Quad(2) )                      ::  sinsqr_t
REAL(idp),  DIMENSION( 1:Num_Quad(2) )                      ::  cossqr_t

REAL(idp),  DIMENSION( 1:Num_Quad(3) )                      ::  sinsqr_p
REAL(idp),  DIMENSION( 1:Num_Quad(3) )                      ::  cossqr_p


REAL(idp),  DIMENSION( 1:Num_Quad(1) )                      ::  Cur_R_Locs
REAL(idp),  DIMENSION( 1:Num_Quad(2) )                      ::  Cur_T_Locs
REAL(idp),  DIMENSION( 1:Num_Quad(3) )                      ::  Cur_P_Locs

INTEGER                                                     ::  Here

REAL(idp)                                                   ::  Energy
REAL(idp)                                                   ::  Pressure
REAL(idp)                                                   ::  Spec_Ent

REAL(idp)                                                   ::  Value

REAL(idp)                                                   ::  A,  B,  C
REAL(idp)                                                   ::  AA, BB, CC

Tot_Quad = Num_Quad(1) * Num_Quad(2) * Num_Quad(3)


RQ_Loc = Map_To_X_Space(Left_Limit, Right_Limit, R_Quad )
TQ_Loc = Map_To_X_Space(Left_Limit, Right_Limit, T_Quad )
PQ_Loc = Map_To_X_Space(Left_Limit, Right_Limit, P_Quad )

MLS_SemiMinor = SemiMinor
MLS_SemiMajor = SemiMajor

CALL Set_MLS_Parameters( Spheroid_Type, SemiMajor, SemiMinor, Rho_O )

CALL Calc_MLS_ABCs( A, B, C, AA, BB, CC )

!IF ( myID_Poseidon == MasterID_Poseidon ) THEN
!IF ( lPF_IO_Flags(iPF_IO_Print_Setup)   ) THEN
!    WRITE(*,'(A)')'------------- Test Parameters ----------------'
!    WRITE(*,'(A)')' Source Configuration : MacLaurin Spheroid'
!    WRITE(*,'(A,A)')      ' - Spheroid Type  : ', MLS_SphereName(MLS_SphereType)
!    WRITE(*,'(A,ES12.5)') ' - Semimajor Axis : ', MLS_SemiMajor
!    WRITE(*,'(A,ES12.5)') ' - Semiminor Axis : ', MLS_SemiMinor
!    WRITE(*,'(A,ES12.5)') ' - Density        : ', MLS_Rho
!    WRITE(*,'(/)')
!END IF
!END IF

MLS_Rho  = Rho_O
Pressure = 0.0_idp
Energy   = 0.0_idp

!Pressure = kappa * MLS_Rho**Gamma
!Energy   = Pressure / (Gamma - 1.0_dip )
Spec_Ent = C_Square + (Energy + Pressure)/MLS_Rho





DO pe = 1,Num_Elem(3)
DO te = 1,Num_Elem(2)
DO re = 1,Num_Elem(1)
    
    Cur_R_Locs(:) = Map_From_X_Space(R_locs(re),R_locs(re+1),RQ_Loc)
    Cur_T_Locs(:) = Map_From_X_Space(T_locs(te),T_locs(te+1),TQ_Loc)
    Cur_P_Locs(:) = Map_From_X_Space(P_locs(pe),P_locs(pe+1),PQ_Loc)

    rsqr(:)     = Cur_R_Locs(:) * Cur_R_Locs(:)
    sinsqr_t(:) = SIN( Cur_T_Locs(:) ) * SIN( CUR_T_LOCS(:) )
    cossqr_t(:) = COS( Cur_T_Locs(:) ) * COS( CUR_T_LOCS(:) )
    sinsqr_p(:) = SIN( Cur_P_Locs(:) ) * SIN( CUR_P_LOCS(:) )
    cossqr_p(:) = COS( Cur_P_Locs(:) ) * COS( CUR_P_LOCS(:) )


    DO rd = 1,Num_Quad(1)
    DO td = 1,Num_Quad(2)
    DO pd = 1,Num_Quad(3)

        Here = Quad_Map(rd, td, pd, Num_Quad )


        Value = ( rsqr(rd) * cossqr_p(pd) * sinsqr_t(td) ) / AA      &
              + ( rsqr(rd) * sinsqr_p(pd) * sinsqr_t(td) ) / BB      &
              + ( rsqr(rd) * cossqr_t(td) ) / CC


        
        IF ( Value .LE. 1.0_idp ) THEN
             
            Output_E(Here,re-1,te-1,pe-1) = MLS_Rho * Spec_Ent - Pressure
        ELSE
            Output_E(Here,re-1,te-1,pe-1) = 0.0_idp
        END IF

    END DO ! pd
    END DO ! td
    END DO ! rd
END DO ! pe
END DO ! te
END DO ! re



END SUBROUTINE Initialize_MacLaurin_Sources







!+101+###########################################################################!
!                                                                                !
!                  Initialize_MacLaurin_Sources                                  !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_MacLaurin_Density( Rho_O, Spheroid_Type,              &
                                         SemiMajor, SemiMinor,              &
                                         Num_Quad,                          &
                                         R_Quad, T_Quad, P_Quad,            &
                                         Left_Limit, Right_Limit,           &
                                         Num_Elem,                          &
                                         R_Locs, T_Locs, P_Locs,            &
                                         Output_Rho      )


REAL(idp),  INTENT(IN)                                                      :: Rho_O
CHARACTER(LEN=7), INTENT(IN)                                                :: Spheroid_Type

REAL(idp),  INTENT(IN)                                                      :: SemiMajor
REAL(idp),  INTENT(IN)                                                      :: SemiMinor

INTEGER,    INTENT(IN), DIMENSION( 1:3 )                                    :: Num_Quad
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Quad(1) )                          :: R_Quad
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Quad(2) )                          :: T_Quad
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Quad(3) )                          :: P_Quad
REAL(idp),  INTENT(IN)                                                      :: Left_Limit
REAL(idp),  INTENT(IN)                                                      :: Right_Limit

INTEGER,    INTENT(IN), DIMENSION( 1:3 )                                    :: Num_Elem
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Elem(1)+1 )                        :: R_Locs
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Elem(2)+1 )                        :: T_Locs
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Elem(3)+1 )                        :: P_Locs

REAL(idp),  INTENT(OUT),DIMENSION( 1:Num_Quad(1)*Num_Quad(2)*Num_Quad(3),   &
                                   0:Num_Elem(1)-1,                         &
                                   0:Num_Elem(2)-1,                         &
                                   0:Num_Elem(3)-1  )                       ::  Output_Rho



INTEGER                                                     ::  re, te, pe
INTEGER                                                     ::  rd, td, pd
INTEGER                                                     ::  Tot_Quad

REAL(idp),  DIMENSION( 1:Num_Quad(1) )                      ::  RQ_Loc
REAL(idp),  DIMENSION( 1:Num_Quad(2) )                      ::  TQ_Loc
REAL(idp),  DIMENSION( 1:Num_Quad(3) )                      ::  PQ_Loc

REAL(idp),  DIMENSION( 1:Num_Quad(1) )                      ::  rsqr

REAL(idp),  DIMENSION( 1:Num_Quad(2) )                      ::  sinsqr_t
REAL(idp),  DIMENSION( 1:Num_Quad(2) )                      ::  cossqr_t

REAL(idp),  DIMENSION( 1:Num_Quad(3) )                      ::  sinsqr_p
REAL(idp),  DIMENSION( 1:Num_Quad(3) )                      ::  cossqr_p


REAL(idp),  DIMENSION( 1:Num_Quad(1) )                      ::  Cur_R_Locs
REAL(idp),  DIMENSION( 1:Num_Quad(2) )                      ::  Cur_T_Locs
REAL(idp),  DIMENSION( 1:Num_Quad(3) )                      ::  Cur_P_Locs

INTEGER                                                     ::  Here

REAL(idp)                                                   ::  Value

REAL(idp)                                                   ::  A,  B,  C
REAL(idp)                                                   ::  AA, BB, CC

Tot_Quad = Num_Quad(1) * Num_Quad(2) * Num_Quad(3)

RQ_Loc = Map_To_X_Space(Left_Limit, Right_Limit, R_Quad )
TQ_Loc = Map_To_X_Space(Left_Limit, Right_Limit, T_Quad )
PQ_Loc = Map_To_X_Space(Left_Limit, Right_Limit, P_Quad )


CALL Set_MLS_Parameters( Spheroid_Type, SemiMajor, SemiMinor, Rho_O)
CALL Calc_MLS_ABCs( A, B, C, AA, BB, CC )

CALL Calc_MLS_Ecc( AA, CC )



!IF ( myID_Poseidon == MasterID_Poseidon ) THEN
!IF ( lPF_IO_Flags(iPF_IO_Print_Setup)   ) THEN
!    WRITE(*,'(A)')'------------- Test Parameters ----------------'
!    WRITE(*,'(A)')' Source Configuration : MacLaurin Spheroid'
!    WRITE(*,'(A,A)')      ' - Spheroid Type  : ', MLS_SphereName(MLS_SphereType)
!    WRITE(*,'(A,ES12.5)') ' - Semimajor Axis : ', MLS_SemiMajor
!    WRITE(*,'(A,ES12.5)') ' - Semiminor Axis : ', MLS_SemiMinor
!    WRITE(*,'(A,ES12.5)') ' - Density        : ', MLS_Rho
!    WRITE(*,'(/)')
!END IF
!END IF




DO pe = 1,Num_Elem(3)
DO te = 1,Num_Elem(2)
DO re = 1,Num_Elem(1)
!    PRINT*,R_Locs(re)
    
    Cur_R_Locs(:) = Map_From_X_Space(R_locs(re),R_locs(re+1),RQ_Loc)
    Cur_T_Locs(:) = Map_From_X_Space(T_locs(te),T_locs(te+1),TQ_Loc)
    Cur_P_Locs(:) = Map_From_X_Space(P_locs(pe),P_locs(pe+1),PQ_Loc)

    rsqr(:)     = Cur_R_Locs(:) * Cur_R_Locs(:)
    sinsqr_t(:) = SIN( Cur_T_Locs(:) ) * SIN( CUR_T_LOCS(:) )
    cossqr_t(:) = COS( Cur_T_Locs(:) ) * COS( CUR_T_LOCS(:) )
    sinsqr_p(:) = SIN( Cur_P_Locs(:) ) * SIN( CUR_P_LOCS(:) )
    cossqr_p(:) = COS( Cur_P_Locs(:) ) * COS( CUR_P_LOCS(:) )


    DO rd = 1,Num_Quad(1)
    DO td = 1,Num_Quad(2)
    DO pd = 1,Num_Quad(3)

        Here = Quad_Map(rd, td, pd, Num_Quad )


        Value = ( rsqr(rd) * cossqr_p(pd) * sinsqr_t(td) ) / AA      &
              + ( rsqr(rd) * sinsqr_p(pd) * sinsqr_t(td) ) / BB      &
              + ( rsqr(rd) * cossqr_t(td) ) / CC


        
        IF ( Value .LE. 1.0_idp ) THEN
             
            Output_Rho(Here,re-1,te-1,pe-1) = MLS_Rho
        ELSE
            Output_Rho(Here,re-1,te-1,pe-1) = 0.0_idp
        END IF

!        PRINT*,Here,re-1,te-1,pe-1,Output_E(Here,re-1,te-1,pe-1)

    END DO ! pd
    END DO ! td
    END DO ! rd
END DO ! pe
END DO ! te
END DO ! re


END SUBROUTINE Initialize_MacLaurin_Density





                                                       
!+301+###########################################################################!
!                                                                                !
!                               MacLaurin_Radius                                 !
!                                                                                !
!################################################################################!
PURE FUNCTION MacLaurin_Radius( theta, phi, A, B, C, AA, BB, CC )


REAL(idp),  INTENT(IN)          ::  theta, phi
REAL(idp),  INTENT(IN)          ::  A,  B,  C
REAL(idp),  INTENT(IN)          ::  AA, BB, CC

REAL(idp)                       ::  RADIUS_HERE

REAL(idp)                       ::  MacLaurin_Radius



IF (MLS_SphereType == 1 ) THEN

    RADIUS_HERE = A * C                                     &
                / SQRT( CC * SIN(theta) * SIN(theta)        &
                      + AA * COS(theta) * COS(theta)     )

ELSE IF (MLS_SphereType == 2) THEN

    RADIUS_HERE = A * B                                                     &
                * SQRT( BB                                                  &
                        * COS(theta) * COS(theta) * SIN(phi) * SIN(phi)     &
                      + AA                                                  &
                        * ( SIN(phi) * SIN(phi) * SIN(theta) * SIN(theta)   &
                            + COS(theta) * COS(theta)   )   )

END IF


MacLaurin_Radius = RADIUS_HERE

END FUNCTION MacLaurin_Radius




!+302+#################################################################
!
!   MacLaurin_Root_Finder
!
!#######################################################################
PURE FUNCTION MacLaurin_Root_Finder( r, theta, phi, AA, CC )


REAL(KIND = idp),   INTENT(IN)                      ::  r, theta, phi
REAL(KIND = idp),   INTENT(IN)                      ::  AA, CC

REAL(KIND = idp)                                    ::  MacLaurin_Root_Finder

REAL(KIND = idp)                                    ::  rsqr, xsqr, ysqr, zsqr

REAL(KIND = idp)                                    ::  LenA, LenB

!! r^2  !!
rsqr = r * r


!! x^2 = r^2 * cos(theta)^2 * sin(phi)^2 !!
xsqr = rsqr * cos(phi) * cos(phi)* sin(theta) * sin(theta)

!! y^2 = r^2 * sin(theta)^2 * sin(phi)^2 !!
ysqr = rsqr * sin(phi) * sin(phi)* sin(theta) * sin(theta)

!! z^2 = r^2 * cos(phi)^2  !!
zsqr = rsqr * cos(theta) * cos(theta)



LenA = AA + CC - xsqr - ysqr - zsqr


IF ( MLS_SphereType == iMLS_Oblate ) THEN

    LenB = AA*(CC - zsqr) - CC*(xsqr + ysqr)

ELSE IF ( MLS_SphereType == iMLS_Prolate ) THEN

    LenB = CC*(AA - xsqr) - AA*(ysqr + zsqr)

END IF


MacLaurin_Root_Finder = 0.5 * ( - LenA + sqrt(LenA*LenA - 4*LenB))




END FUNCTION MacLaurin_Root_Finder










 !+301+################################################################!
!                                                                       !
!          Set_MLS_Parameters                                           !
!                                                                       !
 !#####################################################################!
SUBROUTINE Set_MLS_Parameters( SphereType, SemiMajor, SemiMinor, Rho )


CHARACTER(LEN=7),   INTENT(IN)      :: SphereType
REAL(idp),          INTENT(IN)      :: SemiMajor
REAL(idp),          INTENT(IN)      :: SemiMinor
REAL(idp),          INTENT(IN)      :: Rho

REAL(idp)                           ::  A, B, C
REAL(idp)                           ::  AA, BB, CC

IF ( SphereType == 'Prolate') THEN
    MLS_SphereType  = iMLS_Prolate
ELSE
    MLS_SphereType  = iMLS_Oblate
END IF

MLS_SemiMajor = SemiMajor
MLS_SemiMinor = SemiMinor
MLS_Rho       = Rho


CALL Calc_MLS_ABCs( A, B, C, AA, BB, CC )

print*,"Beofer ECC"
CALL Calc_MLS_Ecc( AA, CC )

IF ( myID_Poseidon == MasterID_Poseidon ) THEN
IF ( lPF_IO_Flags(iPF_IO_Print_Setup)   ) THEN
    WRITE(*,'(A)')'------------- Test Parameters ----------------'
    WRITE(*,'(A)')' Source Configuration : MacLaurin Spheroid'
    WRITE(*,'(A,A)')      ' - Spheroid Type  : ', MLS_SphereName(MLS_SphereType)
    WRITE(*,'(A,ES12.5)') ' - Semimajor Axis : ', MLS_SemiMajor
    WRITE(*,'(A,ES12.5)') ' - Semiminor Axis : ', MLS_SemiMinor
    WRITE(*,'(A,ES12.5)') ' - Density        : ', MLS_Rho
    WRITE(*,'(/)')
END IF
END IF

END SUBROUTINE Set_MLS_Parameters



 !+302+################################################################!
!                                                                       !
!          Calc_MLS_ABCs                                                !
!                                                                       !
 !#####################################################################!
SUBROUTINE Calc_MLS_ABCs( A, B, C, AA, BB, CC )

REAL(idp), INTENT(OUT)          ::  A,  B,  C,  &
                                    AA, BB, CC


IF ( MLS_SphereType == iMLS_Prolate ) THEN ! Prolate

    A = MLS_SemiMajor
    B = MLS_SemiMinor
    C = B
    
ELSE                            ! Oblate

    A = MLS_SemiMajor
    B = A
    C = MLS_SemiMinor

END IF


AA = A*A
BB = B*B
CC = C*C

END SUBROUTINE Calc_MLS_ABCs



 !+302+################################################################!
!                                                                       !
!          Calc_MLS_ABCs                                                !
!                                                                       !
 !#####################################################################!
SUBROUTINE Calc_MLS_Ecc( AA, CC )

REAL(idp), INTENT(IN)          ::  AA, CC

PRINT*,"AA,CC",AA,CC
print*,MLS_SphereType,iMLS_Oblate,imLS_prolate
IF ( MLS_SphereType == iMLS_Oblate ) THEN
    print*,"IN"
    MLS_Ecc = sqrt(1.0_idp - CC/AA)
ELSE IF (MLS_SphereType == iMLS_Prolate) THEN
    print*,"Out"
    MLS_Ecc = sqrt(1.0_idp - AA/CC)
END IF

PRINT*,MLS_Ecc
END SUBROUTINE Calc_MLS_Ecc

END MODULE External_MLS_Profile_Module
