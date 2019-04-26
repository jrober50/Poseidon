!+102+###########################################################################!
!                                                                                !
!                  CREATE_3D_RHS_VECTOR                                          !
!                                                                                !
!################################################################################!
SUBROUTINE CREATE_3D_RHS_VECTOR(   re, te, pe,                                      &
                                    DELTAR_OVERTWO, DELTAT_OVERTWO, DELTAP_OVERTWO, &
                                    SIN_VAL, R_SQUARE,                              &
                                    RHS_TERMS                                       )



INTEGER, INTENT(IN)                                                     ::  re, te, pe


REAL(KIND = idp), INTENT(IN)                                            ::  DELTAR_OVERTWO, &
                                                                            DELTAT_OVERTWO, &
                                                                            DELTAP_OVERTWO


REAL(KIND = idp), INTENT(IN), DIMENSION(1:Source_Degrees(1))            ::  R_SQUARE

REAL(KIND = idp), INTENT(IN), DIMENSION(1:Source_Degrees(2))            ::  SIN_VAL


REAL(KIND = idp), INTENT(IN), DIMENSION( 1:5,                        &
                                         1:Source_Degrees(1),        &
                                         1:Source_Degrees(2),        &
                                         1_Source_Degrees(3)         )  ::  RHS_TERMS







REAL(KIND = idp), INTENT(IN), DIMENSION(-L_LIMIT:L_LIMIT)               :: M_POWER_TABLE




COMPLEX(KIND = idp)                                                     :: Ylm_Conjugate


INTEGER                                                                 ::  pd, td, rd,     &
                                                                            l, m, d

INTEGER                                                                 ::  L_SHIFT
INTEGER                                                                 ::  Current_i_Location




REAL(KIND = idp), DIMENSION(1:Source_Degrees(1),        &
                            1:Source_Degrees(2),        &
                            1_Source_Degrees(3)         )               ::  Int_Factor


COMPLEX(KIND = idp), DIMENSION(1:5)                                     ::  RHS_TMP





L_SHIFT = 5*(L_LIMIT+1)*(L_LIMIT+1)




!$OMP PARALLEL DEFAULT(none)                                            &
!$OMP PRIVATE( pd, td, rd, d, l, m                                      &
!$OMP           Current_i_Location, RHS_TMP                         )   &
!$OMP SHARED( re, te, pe,                                               &
!$OMP           Int_Factor,                                             &
!$OMP           DEGREE, L_LIMIT,                                        &
!$OMP           M_POWER_TABLE, Ylm_Table, Lagrange_Poly_Table,          &
!$OMP           RHS_Vector                                          )   &







!$OMP DO SCHEDULE(dynamic), COLLAPSE(3)
DO d = 0,DEGREE
    DO l = 0,L_LIMIT
        DO m = -l,l

            
            Current_i_Location = (re*DEGREE + d)*L_SHIFT  &
                                + 5 *(l*(l+1) + m)

            DO pd = 1,Source_Degrees(3)
                DO td = 1,Source_Degrees(2)
                    DO rd = 1,Source_Degrees(1)


                        RHS_TMP(1:5) = RHS_TERMS(1:5, rd, td, pd)                           &
                                            * M_POWER_TABLE(m)                              &
                                            * Ylm_Table(-m, l, td, pd, te, pe)              &
                                            * Lagrange_Poly_Table(d, rd, 0)                 &
                                            * Int_Factor(rd, td, pd)


                    END DO  ! rd Loop
                END DO  ! td Loop
            END DO  ! pd Loop


            RHS_Vector(Current_i_Location:Current_i_Location+4) = RHS_TMP(1:5)


        END DO  ! m Loop
    END DO  ! l Loop
END DO  ! d Loop
!$OMP END DO




END SUBROUTINE CREATE_3D_RHS_VECTOR
