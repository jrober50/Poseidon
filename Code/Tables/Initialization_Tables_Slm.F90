   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Initialization_Tables_Slm                                                    !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +201+   Initialize_Ylm_Tables                                               !##!
!##!    +202+   Initialize_Lagrange_Poly_Tables                                     !##!
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

USE Poseidon_Message_Routines_Module, &
            ONLY :  Init_Message,           &
                    Warning_Message

USE Poseidon_Parameters, &
            ONLY :  Domain_Dim,             &
                    Degree,                 &
                    L_Limit,                &
                    Verbose_Flag

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Max_Grid_Size,          &
                    AMReX_Num_Levels

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Num_T_Quad_Points,      &
                    Num_P_Quad_Points,      &
                    Int_R_Locations,        &
                    Int_T_Locations,        &
                    Int_P_Locations

USE Variables_Mesh, &
            ONLY :  R_Inner,                &
                    R_Outer,                &
                    Num_R_Elements,         &
                    Num_T_Elements,         &
                    Num_P_Elements,         &
                    rlocs,                  &
                    tlocs,                  &
                    plocs

USE Variables_Derived, &
            ONLY :  LM_Length,              &
                    LM_Short_Length


USE Variables_Tables, &
            ONLY :  Slm_Elem_Values,        &
                    Slm_Elem_dt_Values,     &
                    Slm_Elem_dp_Values,     &
                    Plm_Values,             &
                    Plm_dt_Values,          &
                    Nlm_Values,             &
                    Am_Values,              &
                    AM_dp_Values,           &
                    Level_dx,               &
                    Level_Ratios


USE Functions_Math, &
            ONLY :  Lagrange_Poly,          &
                    Lagrange_Poly_Deriv,    &
                    Legendre_Poly,          &
                    Legendre_Poly_Array,    &
                    Norm_Factor,            &
                    Sqrt_Factor

USE Maps_Quadrature, &
            ONLY :  Map_To_tpd

USE Maps_X_Space, &
            ONLY :  Map_From_X_Space

USE Maps_Domain, &
            ONLY :  Map_To_lm,              &
                    Map_To_Short_lm



IMPLICIT NONE

CONTAINS

!+201+##########################################################################!
!                                                                               !
!                  Initialize_Slm_Tables                                        !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_Slm_Tables()

IF ( Verbose_Flag ) CALL Init_Message('Intializing Slm Tables.')

! Initialize Nlm Table
CALL Initialize_Nlm_Table(  L_Limit,                &
                            LM_Short_Length,        &
                            Nlm_Values              )

! Initialize Am Table
CALL Initialize_Am_Tables(  Num_P_Quad_Points,      &
                            Int_P_Locations,        &
                            L_Limit,                &
                            Num_P_Elements,         &
                            [0, Num_P_Elements-1],  &
                            plocs,                  &
                            Am_Values,              &
                            Am_dp_Values            )

! Initialize Plm Table
CALL Initialize_Plm_Tables( Num_T_Quad_Points,      &
                            Int_T_Locations,        &
                            L_Limit,                &
                            LM_Short_Length,        &
                            Num_T_Elements,         &
                            [0, Num_T_Elements-1],  &
                            tlocs,                  &
                            Plm_Values,             &
                            Plm_dt_Values           )


END SUBROUTINE Initialize_Slm_Tables





!+203+##########################################################################!
!                                                                               !
!          Initialize_Slm_Tables_on_Elem                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_Slm_Tables_on_Elem( te, pe,               &
                                          nTQ, nPQ,             &
                                          E_Num,                &
                                          E_Lo,                 &
                                          Plm_Table,            &
                                          Plm_dt_Table,         &
                                          Am_Table,             &
                                          Am_dp_Table,          &
                                          Slm_Elem_Table,       &
                                          Slm_Elem_dt_Table,    &
                                          Slm_Elem_dp_Table     )

INTEGER, INTENT(IN)                                                 ::  te, pe
INTEGER, INTENT(IN)                                                 ::  nTQ, nPQ

INTEGER,    DIMENSION(3),                           INTENT(IN)      ::  E_Num
INTEGER,    DIMENSION(3),                           INTENT(IN)      ::  E_Lo

REAL(idp),  DIMENSION( 1:nTQ,                   &
                       1:LM_Short_Length,       &
                       0:E_Num(2)-1             ),  INTENT(IN)      ::  Plm_Table

REAL(idp),  DIMENSION( 1:nTQ,                   &
                       1:LM_Short_Length,       &
                       0:E_Num(2)-1             ),  INTENT(IN)      ::  Plm_dt_Table
                      
REAL(idp),  DIMENSION( 1:nPQ,                   &
                       -L_Limit:L_Limit,        &
                       0:E_Num(3)-1             ),  INTENT(IN)      ::  Am_Table
                      
REAL(idp),  DIMENSION( 1:nPQ,                   &
                       -L_Limit:L_Limit,        &
                       0:E_Num(3)-1             ),  INTENT(IN)      ::  Am_dp_Table
                      
                      
REAL(idp),  DIMENSION( 1:LM_Length,             &
                       1:nTQ*nPQ                ),  INTENT(OUT)     ::  Slm_Elem_Table
                      
REAL(idp),  DIMENSION( 1:LM_Length,             &
                       1:nTQ*nPQ                ),  INTENT(OUT)     ::  Slm_Elem_dt_Table
                      
REAL(idp),  DIMENSION( 1:LM_Length,             &
                       1:nTQ*nPQ                ),  INTENT(OUT)     ::  Slm_Elem_dp_Table


INTEGER                                 ::  l, m, td, pd
INTEGER                                 ::  tpd, lm, short_lm
INTEGER                                 ::  te_Offset
INTEGER                                 ::  pe_Offset

te_Offset = te-E_Lo(2)
pe_Offset = pe-E_Lo(3)


DO l = 0,L_Limit
DO m = -l,l
DO td = 1,nTQ
DO pd = 1,nPQ

    tpd = Map_To_tpd(td,pd,nPQ)
    lm = Map_To_lm(l,m)
    short_lm = Map_To_Short_LM(l,abs(m))
    Slm_Elem_Table(lm,tpd) = Nlm_Values(short_lm)                      &
                            * Plm_Table(td,short_lm,te_Offset)         &
                            * Am_Table(pd,m,pe_Offset)
                            
    Slm_Elem_dt_Table(lm,tpd) = Nlm_Values(short_lm)                   &
                               * Plm_dt_Table(td,short_lm,te_Offset)   &
                               * Am_Table(pd,m,pe_Offset)
                               
    Slm_Elem_dp_Table(lm,tpd) = Nlm_Values(short_lm)                   &
                               * Plm_Table(td,short_lm,te_Offset)      &
                               * Am_dp_Table(pd,m,pe_Offset)


!    PRINT*,te,pe,tpd,lm,Nlm_Values(short_lm) ,Plm_Table(td,short_lm,te_Offset),Am_Table(pd,m,pe_Offset)

!    PRINT*,te,pe,tpd,lm,Plm_Table(td,short_lm,te_Offset),       &
!                        Plm_dt_Table(td,short_lm,te_Offset),    &
!                        Am_Table(pd,m,pe_Offset),               &
!                        Am_dp_Table(pd,m,pe_Offset)
                        
!    PRINT*,te,pe,tpd,lm,Slm_Elem_dp_Table(lm,tpd),Nlm_Values(short_lm),                   &
!                               Plm_Table(td,short_lm,te_Offset),      &
!                               Am_dp_Table(pd,m,pe_Offset)

END DO ! m
END DO ! l
END DO ! pd
END DO ! td




!STOP "At end of Initialize_Slm_Tables_on_Elem"
END SUBROUTINE Initialize_Slm_Tables_on_Elem




!+203+##########################################################################!
!                                                                               !
!          Initialize_Ylm_Tables_on_Elem                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_Slm_Table_on_Elem(  te, pe,               &
                                          nTQ, nPQ,             &
                                          E_Num,                &
                                          E_Lo,                 &
                                          Plm_Table,            &
                                          Am_Table,             &
                                          Slm_Elem_Table        )

INTEGER,                                            INTENT(IN)      ::  te, pe
INTEGER,                                            INTENT(IN)      ::  nTQ, nPQ
INTEGER,    DIMENSION(3),                           INTENT(IN)      ::  E_Num
INTEGER,    DIMENSION(3),                           INTENT(IN)      ::  E_Lo

REAL(idp),  DIMENSION( 1:nTQ,                   &
                       1:LM_Short_Length,       &
                       0:E_Num(2)-1             ),  INTENT(IN)      ::  Plm_Table
                      
                      
REAL(idp),  DIMENSION( 1:nPQ,                   &
                       -L_Limit:L_Limit,        &
                       0:E_Num(3)-1             ),  INTENT(IN)      ::  Am_Table

REAL(idp),  DIMENSION( 1:LM_Length,             &
                       1:nTQ*nPQ                ),  INTENT(INOUT)   ::  Slm_Elem_Table

INTEGER                                 ::  l, m, td, pd
INTEGER                                 ::  tpd, lm, short_lm

INTEGER                                 ::  te_Offset
INTEGER                                 ::  pe_Offset

te_Offset = te-E_Lo(2)
pe_Offset = pe-E_Lo(3)

DO l = 0,L_Limit
DO m = -l,l
DO td = 1,nTQ
DO pd = 1,nPQ

    tpd = Map_To_tpd(td,pd)
    lm = Map_To_lm(l,m)
    short_lm = Map_To_Short_LM(l,abs(m))
                            
                            
    Slm_Elem_Table(lm,tpd) = Nlm_Values(short_lm)               &
                           * Plm_Table(td,short_lm,te_Offset)   &
                           * Am_Table(pd,m,pe_Offset)
                    
    
END DO ! m
END DO ! l
END DO ! pd
END DO ! td


END SUBROUTINE Initialize_Slm_Table_on_Elem








!+201+##########################################################################!
!                                                                               !
!                  Initialize_Alm_Table                                         !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_Am_Tables(    QP_Num,     &
                                    QP_xLocs,   &
                                    L_Max,      &
                                    E_Num,      &
                                    E_LoHi,     &
                                    E_Locs,     &
                                    Am_Table,   &
                                    Am_dp_Table )

INTEGER,                                                INTENT(IN)  ::  QP_Num
REAL(idp),  DIMENSION(1:QP_Num),                        INTENT(IN)  ::  QP_xLocs
INTEGER,                                                INTENT(IN)  ::  L_Max
INTEGER,                                                INTENT(IN)  ::  E_Num
INTEGER,    DIMENSION(1:2),                             INTENT(IN)  ::  E_LoHi
REAL(idp),  DIMENSION(0:E_Num),                         INTENT(IN)  ::  E_Locs

REAL(idp),  DIMENSION(1:QP_Num,-L_Max:L_Max,0:E_Num-1), INTENT(OUT) ::  Am_Table
REAL(idp),  DIMENSION(1:QP_Num,-L_Max:L_Max,0:E_Num-1), INTENT(OUT) ::  AM_dp_Table



REAL(idp),  DIMENSION(1:QP_Num,0:L_Max,0:E_Num-1)                   ::  cosm
REAL(idp),  DIMENSION(1:QP_Num,0:L_Max,0:E_Num-1)                   ::  sinm


INTEGER                                                         ::  Elem, Quad
INTEGER                                                         ::  m

REAL(idp)                                                       ::  sqrt_two



! Initialize temporary trig tables
CALL Initialize_Trig_Tables(QP_Num,     &
                            QP_xLocs,   &
                            L_Max,      &
                            E_Num,      &
                            E_LoHi,     &
                            E_Locs,     &
                            Cosm,       &
                            Sinm        )


sqrt_two = SQRT(2.0_idp)

Am_Table(:,0,:) = 1.0_idp
Am_dp_Table(:,0,:) = 0.0_idp

IF ( L_Max > 0 ) THEN

    DO Elem = 0,E_Num-1
    
        DO m  = -L_Max,-1
        DO Quad = 1,QP_Num
            Am_Table(Quad,m,Elem) = sqrt_two*sinm(Quad,-m,Elem)
        END DO ! pd Loop
        END DO ! m  Loop



        DO m  = 1,L_Max
        DO Quad = 1,QP_Num
            Am_Table(Quad,m,Elem) = sqrt_two*cosm(Quad,m,Elem)
        END DO ! Quad Loop
        END DO ! m  Loop

    END DO ! Elem Loop


    DO Elem = 0,E_Num-1
    
        DO m  = -L_Max,-1
        DO Quad = 1,QP_Num
            Am_dp_Table(Quad,m,Elem) = sqrt_two*m*cosm(Quad,-m,Elem)
        END DO ! Quad Loop
        END DO ! m  Loop

        DO m  = 1,L_Max
        DO Quad = 1,QP_Num
            Am_dp_Table(Quad,m,Elem) = -sqrt_two*abs(m)*sinm(Quad,m,Elem)
        END DO ! Quad Loop
        END DO ! m  Loop

    END DO ! Elem Loop

END IF

!STOP "At end of Initialize_Am_Tables"

END SUBROUTINE Initialize_Am_Tables



!+201+##########################################################################!
!                                                                               !
!                  Initialize_Alm_Table                                         !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_Am_Table(     QP_Num,     &
                                    QP_xLocs,   &
                                    L_Max,      &
                                    E_Num,      &
                                    E_LoHi,     &
                                    E_Locs,     &
                                    Am_Table    )

INTEGER,                                                INTENT(IN)  ::  QP_Num
REAL(idp),  DIMENSION(1:QP_Num),                        INTENT(IN)  ::  QP_xLocs
INTEGER,                                                INTENT(IN)  ::  L_Max
INTEGER,                                                INTENT(IN)  ::  E_Num
INTEGER,    DIMENSION(1:2),                             INTENT(IN)  ::  E_LoHi
REAL(idp),  DIMENSION(0:E_Num),                         INTENT(IN)  ::  E_Locs

REAL(idp),  DIMENSION(1:QP_Num,-L_Max:L_Max,0:E_Num-1), INTENT(OUT) ::  Am_Table



REAL(idp),  DIMENSION(1:QP_Num,0:L_Max,0:E_Num)                   ::  cosm
REAL(idp),  DIMENSION(1:QP_Num,0:L_Max,0:E_Num)                   ::  sinm


INTEGER                                                         ::  Elem, Quad
INTEGER                                                         ::  m

REAL(idp)                                                       ::  sqrt_two



! Initialize temporary trig tables
CALL Initialize_Trig_Tables(QP_Num,     &
                            QP_xLocs,   &
                            L_Max,      &
                            E_Num,      &
                            E_LoHi,     &
                            E_Locs,     &
                            Cosm,       &
                            Sinm        )




sqrt_two = SQRT(2.0_idp)

Am_Table(:,0,:) = 1.0_idp

IF ( L_Max > 0 ) THEN

    DO Elem = 0,E_Num
    
        DO m  = -L_Max,-1,-1
        DO Quad = 1,QP_Num
            Am_Table(Quad,m,Elem) = sqrt_two*sinm(Quad,-m,Elem)
        END DO ! pd Loop
        END DO ! m  Loop



        DO m  = 1,L_Max
        DO Quad = 1,QP_Num
            Am_Table(Quad,m,Elem) = sqrt_two*cosm(Quad,m,Elem)
        END DO ! Quad Loop
        END DO ! m  Loop

    END DO ! Elem Loop

END IF

END SUBROUTINE Initialize_Am_Table







!+201+##########################################################################!
!                                                                               !
!                  Initialize_Plm_Table                                         !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_Plm_Tables(   QP_Num,         &
                                    QP_xLocs,       &
                                    L_Max,          &
                                    LM_Length,      &
                                    E_Num,          &
                                    E_LoHi,         &
                                    E_Locs,         &
                                    Plm_Table,      &
                                    Plm_dt_Table    )

INTEGER,                                                INTENT(IN)  ::  QP_Num
REAL(idp),  DIMENSION(1:QP_Num),                        INTENT(IN)  ::  QP_xLocs
INTEGER,                                                INTENT(IN)  ::  L_Max
INTEGER,                                                INTENT(IN)  ::  LM_Length
INTEGER,                                                INTENT(IN)  ::  E_Num
INTEGER,    DIMENSION(1:2),                             INTENT(IN)  ::  E_LoHi
REAL(idp),  DIMENSION(0:E_Num),                         INTENT(IN)  ::  E_Locs

REAL(idp),  DIMENSION(1:QP_Num,1:LM_Length,0:E_Num-1),  INTENT(OUT) ::  Plm_Table
REAL(idp),  DIMENSION(1:QP_Num,1:LM_Length,0:E_Num-1),  INTENT(OUT) ::  Plm_dt_Table

REAL(idp),  DIMENSION(1:QP_Num,0:L_Max,0:E_Num-1)                   ::  cosm
REAL(idp),  DIMENSION(1:QP_Num,0:L_Max,0:E_Num-1)                   ::  sinm
                    
INTEGER                 ::  Elem
INTEGER                 ::  l, m
INTEGER                 ::  short_lma
INTEGER                 ::  short_lmb
INTEGER                 ::  short_lmc



CALL Initialize_Trig_Tables(QP_Num,     &
                            QP_xLocs,   &
                            L_Max,      &
                            E_Num,      &
                            E_LoHi,     &
                            E_Locs,     &
                            Cosm,       &
                            Sinm        )

Plm_Table(:,:,:) = 0.0_idp
Plm_Table(:,1,:) = 1.0_idp
Plm_dt_Table(:,:,:) = 0.0_idp


IF ( L_Max > 0 ) THEN

    DO Elem = 0,E_Num-1
    DO l = 1,L_Max
    
        Short_lma = Map_To_Short_lm(l,l)
        Short_lmb = Map_To_Short_lm(l-1,l-1)
        
        Plm_Table(:,Short_lma,Elem) = (1.0_idp - 2.0_idp*l)*sinm(:,1,Elem)*Plm_Table(:,short_lmb,Elem)
    END DO  ! l Loop
    END DO  ! Elem Loop


    DO Elem = 0,E_Num-1
    DO l = 0,L_Max-1
        Short_lma = Map_To_Short_lm(l+1,l)
        Short_lmb = Map_To_Short_lm(l,l)
        
        
        Plm_Table(:,Short_lma,Elem) = (2.0_idp*l + 1.0_idp)*cosm(:,1,Elem)*Plm_Table(:,short_lmb,Elem)
    END DO  ! l Loop
    END DO  ! Elem Loop


    DO Elem = 0,E_Num-1
    DO l = 2,L_Max
    DO m = 0,l-1

        Short_lma = Map_To_Short_lm(l,m)
        Short_lmb = Map_To_Short_lm(l-1,m)
        Short_lmc = Map_To_Short_lm(l-2,m)
        
        IF ( m <= l-2 ) THEN
            Plm_Table(:,Short_lma,Elem) = REAL(2.0*l - 1.0_idp, KIND = idp)     &
                                        / REAL(l-m, KIND = idp)                 &
                                        * Cosm(:,1,Elem)                        &
                                        * Plm_Table(:,Short_lmb,Elem)           &
                                        - REAL( l + m - 1.0_idp, Kind = idp)    &
                                        / REAL(l-m, KIND = idp)                 &
                                        * Plm_Table(:,Short_lmc,Elem)
        END IF

    END DO ! m Loop
    END DO ! l Loop
    END DO ! pe Loop


    ! Derivatives of the Associated Legendre Polynomials
    DO Elem = 0,E_Num-1
    DO l = 0,L_Max

        Short_lma = Map_To_Short_lm(l,l)
        Plm_dt_Table(:,Short_lma,Elem) =REAL(l,kind=idp)                &
                                        * Cosm(:,1,Elem)                &
                                        / Sinm(:,1,Elem)                &
                                        * Plm_Table(:,Short_lma,Elem)

    END DO ! l Loop
    END DO ! pe Loop


    DO Elem = 0,E_Num-1
    DO l = 0,L_Max
    DO m = 0,l-1

        Short_lma = Map_To_Short_lm(l,m)
        Short_lmb = Map_To_Short_lm(l-1,m)

        Plm_dt_Table(:,Short_lma,Elem) = REAL(l,kind=idp)                   &
                                          * Cosm(:,1,Elem)                  &
                                          / Sinm(:,1,Elem)                  &
                                          * Plm_Table(:,Short_lma,Elem)     &
                                        - REAL(l+m,kind=idp)                &
                                          / Sinm(:,1,Elem)                  &
                                          * Plm_Table(:,Short_lmb,Elem)
 
    END DO ! m Loop
    END DO ! l Loop
    END DO ! pe Loop

END IF

!DO Elem = 0,E_Num-1
!DO l = 0,L_Max
!DO m = 0,l
!
!    Short_lma = Map_To_Short_lm(l,m)
!
!    PRINT*,E_LoHi(1)+Elem,l,m,                      &
!                cosm(:,1,Elem),                     &
!                sinm(:,1,Elem),                     &
!                Plm_Table(:,Short_lma,Elem),        &
!                Plm_dt_Table(:,Short_lma,Elem)
!
!END DO
!END DO
!END DO

!STOP "At end of Initialize_Plm_Tables"
END SUBROUTINE Initialize_Plm_Tables










!+201+##########################################################################!
!                                                                               !
!                  Initialize_Plm_Table                                         !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_Plm_Table(    QP_Num,         &
                                    QP_xLocs,       &
                                    L_Max,          &
                                    LM_Length,      &
                                    E_Num,          &
                                    E_LoHi,         &
                                    E_Locs,         &
                                    Plm_Table       )

INTEGER,                                                INTENT(IN)  ::  QP_Num
REAL(idp),  DIMENSION(1:QP_Num),                        INTENT(IN)  ::  QP_xLocs
INTEGER,                                                INTENT(IN)  ::  L_Max
INTEGER,                                                INTENT(IN)  ::  LM_Length
INTEGER,                                                INTENT(IN)  ::  E_Num
INTEGER,    DIMENSION(1:2),                             INTENT(IN)  ::  E_LoHi
REAL(idp),  DIMENSION(0:E_Num),                         INTENT(IN)  ::  E_Locs

REAL(idp),  DIMENSION(1:QP_Num,1:LM_Length,0:E_Num-1),  INTENT(OUT) ::  Plm_Table

REAL(idp),  DIMENSION(1:QP_Num,0:L_Max,0:E_Num-1)                   ::  cosm
REAL(idp),  DIMENSION(1:QP_Num,0:L_Max,0:E_Num-1)                   ::  sinm
                    
INTEGER                 ::  Elem
INTEGER                 ::  l, m
INTEGER                 ::  short_lma
INTEGER                 ::  short_lmb
INTEGER                 ::  short_lmc


CALL Initialize_Trig_Tables(QP_Num,     &
                            QP_xLocs,   &
                            1,          &
                            E_Num,      &
                            E_LoHi,     &
                            E_Locs,     &
                            Cosm,       &
                            Sinm        )



Plm_Table(:,1,:) = 1.0_idp


IF ( L_Max > 0 ) THEN

    DO Elem = 0,E_Num-1
    DO l = 1,L_Max
        Short_lma = Map_To_Short_lm(l,l)
        Short_lmb = Map_To_Short_lm(l-1,l-1)
        
        Plm_Table(:,Short_lma,Elem) = (2.0_idp*l - 1.0_idp)*sinm(:,1,Elem)*Plm_Table(:,short_lmb,Elem)
    END DO  ! l Loop
    END DO  ! Elem Loop



    DO Elem = 0,E_Num-1
    DO l = 1,L_Max-1
        Short_lma = Map_To_Short_lm(l+1,l)
        Short_lmb = Map_To_Short_lm(l,l)
        
        Plm_Table(:,Short_lma,Elem) = (2.0_idp*l - 1.0_idp)*cosm(:,1,Elem)*Plm_Table(:,short_lmb,Elem)
    END DO  ! l Loop
    END DO  ! Elem Loop


    DO Elem = 0,E_Num-1
    DO l = 2,L_Max
    DO m = 0,L_Max-2

        Short_lma = Map_To_Short_lm(l,m)
        Short_lmb = Map_To_Short_lm(l-1,m)
        Short_lmc = Map_To_Short_lm(l-2,m)
    
        
        Plm_Table(:,Short_lma,Elem) = REAL(2.0*l - 1.0_idp, KIND = idp)     &
                                    / REAL(l-m, KIND = idp)                 &
                                    * Cosm(:,1,Elem)                        &
                                    * Plm_Table(:,Short_lmb,Elem)           &
                                    - REAL( l + m - 1.0_idp, Kind = idp)    &
                                    / REAL(l-m, KIND = idp)                 &
                                    * Plm_Table(:,Short_lmc,Elem)
        


    END DO ! m Loop
    END DO ! l Loop
    END DO ! pe Loop

END IF

END SUBROUTINE Initialize_Plm_Table






!+201+##########################################################################!
!                                                                               !
!                  Initialize_Nlm_Table                                         !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_Nlm_Table( L_Max,     &
                                 LM_Length, &
                                 Nlm_Table  )

INTEGER,                                    INTENT(IN)  ::  L_Max
INTEGER,                                    INTENT(IN)  ::  LM_Length
REAL(idp),  DIMENSION(1:LM_Length),         INTENT(OUT) ::  Nlm_Table

INTEGER                 :: l, m
INTEGER                 :: Short_lm

DO l = 0,L_Max
DO m = 0,l

    Short_lm = Map_To_Short_lm(l,m)
    Nlm_Table(Short_lm) = Norm_Factor(l,m)

END DO ! m Loop
END DO ! l Loop

END SUBROUTINE Initialize_Nlm_Table








!+201+##########################################################################!
!                                                                               !
!                  Initialize_Trig_Table                                        !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_Trig_Tables(  QP_Num,     &
                                    QP_xLocs,   &
                                    L_Max,      &
                                    E_Num,      &
                                    E_LoHi,     &
                                    E_Locs,     &
                                    Cosm,       &
                                    Sinm        )

INTEGER,                                            INTENT(IN)  ::  QP_Num
REAL(idp),  DIMENSION(1:QP_Num),                    INTENT(IN)  ::  QP_xLocs
INTEGER,                                            INTENT(IN)  ::  L_Max
INTEGER,                                            INTENT(IN)  ::  E_Num
INTEGER,    DIMENSION(1:2),                         INTENT(IN)  ::  E_LoHi
REAL(idp),  DIMENSION(0:E_Num),                     INTENT(IN)  ::  E_Locs

REAL(idp),  DIMENSION(1:QP_Num,0:L_Max,0:E_Num-1),  INTENT(OUT) ::  Cosm
REAL(idp),  DIMENSION(1:QP_Num,0:L_Max,0:E_Num-1),  INTENT(OUT) ::  Sinm



INTEGER                                                         ::  Elem, Quad,Here
INTEGER                                                         ::  m

REAL(idp), DIMENSION(1:QP_Num)                                  ::  Cur_Locs



Cosm(:,0,:) = 1.0_idp
Sinm(:,0,:) = 0.0_idp


IF (L_Max > 0 ) THEN
    DO Elem = E_LoHi(1),E_LoHi(2)
    
        Here = Elem-E_LoHi(1)
        Cur_Locs = Map_From_X_Space(E_Locs(Here), E_Locs(Here + 1), QP_xLocs)
        
        DO Quad = 1,QP_Num

            Cosm(Quad,1,Here) = DCOS(Cur_Locs(Quad))
            Sinm(Quad,1,Here) = DSIN(Cur_Locs(Quad))

        END DO ! pd Loop
    END DO ! pe Loop


    DO Elem = 0,E_Num-1
    DO m    = 2,L_Max
    DO Quad = 1,QP_Num

        Cosm(Quad,m,Elem) = cosm(Quad,1,Elem)*cosm(Quad,m-1,Elem) - sinm(Quad,1,Elem)*sinm(Quad,m-1,Elem)
        Sinm(Quad,m,Elem) = sinm(Quad,1,Elem)*cosm(Quad,m-1,Elem) + cosm(Quad,1,elem)*sinm(Quad,m-1,Elem)

    END DO ! pd Loop
    END DO ! m  Loop
    END DO ! pe Loop
ELSE IF ( L_Max < 0 ) THEN

    CALL Warning_Message('Attempting to intialize Slm Trig tables with L_Max < 0')

END IF

END SUBROUTINE Initialize_Trig_Tables







END MODULE Initialization_Tables_Slm 
