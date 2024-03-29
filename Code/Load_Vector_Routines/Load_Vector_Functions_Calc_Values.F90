   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Load_Vector_Functions_Calc_Values_Module                              !##!
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

USE Poseidon_Parameters, &
            ONLY :  DEGREE

USE Variables_Derived, &
            ONLY :  LM_Length

USE Variables_Mesh, &
            ONLY :  Num_P_Elements
            

USE Variables_Quadrature, &
            ONLY :  NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    NUM_TP_QUAD_POINTS,         &
                    INT_R_LOCATIONS,            &
                    INT_T_LOCATIONS,            &
                    INT_P_LOCATIONS,            &
                    INT_R_WEIGHTS,              &
                    INT_T_WEIGHTS,              &
                    INT_P_WEIGHTS,              &
                    INT_TP_WEIGHTS

USE Variables_Tables, &
            ONLY :  Slm_Elem_Values,            &
                    Slm_Elem_dt_Values,         &
                    Slm_Elem_dp_Values,         &
                    Lagrange_Poly_Table
                    
USE Variables_Vectors, &
            ONLY :  dVA_Coeff_Vector,          &
                    dVB_Coeff_Vector,           &
                    dVA_Coeff_Vector,          &
                    dVB_Coeff_Vector

USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,                      &
                    iU_LF,                      &
                    iU_S1,                      &
                    iU_S2,                      &
                    iU_S3,                      &
                    iU_X1,                      &
                    iU_X2,                      &
                    iU_X2

USE Maps_Fixed_Point, &
            ONLY :  FP_Array_Map_TypeB

USE Maps_Domain, &
            ONLY :  Map_To_FEM_Node,            &
                    FEM_Elem_Map

USE Maps_Quadrature, &
            ONLY :  Map_To_tpd

IMPLICIT NONE


CONTAINS


 !+101+####################################################!
!                                                           !
!                  Calc_Int_Weights                         !
!                                                           !
 !#########################################################!
SUBROUTINE Calc_Int_Weights( DROT, DTOT,            &
                             R_Square, Sin_Val,     &
                             rWeights, tpWeights    )


REAL(idp), INTENT(IN)                                       ::  DROT, DTOT
REAL(idp), INTENT(IN),    DIMENSION(1:Num_R_Quad_Points)    ::  R_Square
REAL(idp), INTENT(IN),    DIMENSION(1:Num_TP_Quad_Points)   ::  Sin_Val
REAL(idp), INTENT(INOUT), DIMENSION(1:Num_R_Quad_Points)    ::  rWeights
REAL(idp), INTENT(INOUT), DIMENSION(1:Num_TP_Quad_Points)   ::  tpWeights

INTEGER                                                     ::  tpd, td, pd


                         !                                                 !
                        !!                                                 !!
                       !!!          Initialize Local Quadratures           !!!
                        !!                                                 !!
                         !                                                 !
rWeights(:) = DROT * R_SQUARE(:) * INT_R_WEIGHTS(:)


DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS

    
   tpd = Map_To_tpd(td,pd)
   tpWeights( tpd ) = SIN_VAL(tpd)                   &
                    * DTOT * INT_T_WEIGHTS(td)      &
                    * INT_P_WEIGHTS(pd)



END DO


END DO

END SUBROUTINE Calc_Int_Weights





 !+02+####################################################!
!                                                           !
!                  Calc_Int_Weights                         !
!                                                           !
 !#########################################################!
SUBROUTINE Calc_Int_Weights_AMReX( DROT, DTOT,              &
                                    R_Square, Sin_Val,      &
                                    rWeights, tpWeights,    &
                                    Level                   )


REAL(idp), INTENT(IN)                                       ::  DROT, DTOT
REAL(idp), INTENT(IN),    DIMENSION(1:Num_R_Quad_Points)    ::  R_Square
REAL(idp), INTENT(IN),    DIMENSION(1:Num_TP_Quad_Points)   ::  Sin_Val
REAL(idp), INTENT(INOUT), DIMENSION(1:Num_R_Quad_Points)    ::  rWeights
REAL(idp), INTENT(INOUT), DIMENSION(1:Num_TP_Quad_Points)   ::  tpWeights
INTEGER,   INTENT(IN)                                       ::  Level


INTEGER                                                     ::  tpd, td, pd
!REAL(idp)                                                   ::  Int_P_Weight

!Int_P_Weight = (2*pi/(Num_P_Elements*2**Level))/Num_P_Quad_Points


                         !                                                 !
                        !!                                                 !!
                       !!!          Initialize Local Quadratures           !!!
                        !!                                                 !!
                         !                                                 !
rWeights(:) = DROT * R_SQUARE(:) * INT_R_WEIGHTS(:)


DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS

    
    
   tpd = Map_To_tpd(td,pd)
   tpWeights( tpd ) = SIN_VAL(tpd)                   &
                    * DTOT * INT_T_WEIGHTS(td)      &
                    * Int_P_Weights(pd)

END DO
END DO

END SUBROUTINE Calc_Int_Weights_AMReX










 !+201+####################################################!
!                                                           !
!          Calc_Val_On_Elem_TypeA                           !
!                                                           !
 !#########################################################!
SUBROUTINE Calc_Val_On_Elem_TypeA( iE, Val, iU, Level )

INTEGER,    INTENT(IN), DIMENSION(3)        ::  iE

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points),   INTENT(OUT)  :: Val
INTEGER,    INTENT(IN)                      :: iU
INTEGER,    INTENT(IN)                      :: Level

REAL(idp)                                   :: TMP_Val
INTEGER                                     :: rd, tpd, d
INTEGER,        DIMENSION(0:DEGREE)         :: Here
INTEGER                                     :: iRE


#ifdef POSEIDON_AMREX_FLAG
    iRE = FEM_Elem_Map(iE(1),Level)
#else
    iRE = iE(1)
#endif





DO d = 0,DEGREE
    Here(d) = Map_To_FEM_Node(iRE,d)
END DO ! d



DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS
    
    TMP_Val = 0.0_idp
    DO d = 0,DEGREE


#ifdef POSEIDON_AMREX_FLAG

        TMP_Val = TMP_Val                               &
            + SUM( dVA_Coeff_Vector( Here(d), :, iU )  &
                   * Slm_Elem_Values( :, tpd )   )      &
            * Lagrange_Poly_Table( d, rd, 0)
            
#else


!        TMP_Val = TMP_Val                                   &
!            + SUM( dVA_Coeff_Vector( Here(d), :, iU )      &
!                   * Slm_Elem_Values( :, tpd ) )   &
!            * Lagrange_Poly_Table( d, rd, 0 )
            
            
        TMP_Val = TMP_Val                                   &
                + SUM( dVA_Coeff_Vector(Here(d),:,iU)       &
                        * Slm_Elem_Values(:,tpd)      )     &
                * Lagrange_Poly_Table(d,rd,0)

#endif


    END DO  ! d

    Val(tpd,rd)         = TMP_Val

END DO ! td
END DO ! rd




END SUBROUTINE Calc_Val_On_Elem_TypeA







 !+202+####################################################!
!                                                           !
!          Calc_Drv_On_Elem_TypeA                           !
!                                                           !
 !#########################################################!
SUBROUTINE Calc_Drv_On_Elem_TypeA( iE, DROT, Drv, iU, Level )

INTEGER,    INTENT(IN), DIMENSION(3)        ::  iE
REAL(idp), INTENT(IN)                       ::  DROT

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3),  INTENT(OUT)  :: Drv
INTEGER,   INTENT(IN)                       ::  iU
INTEGER,    INTENT(IN)                      ::  Level

REAL(idp), DIMENSION(3)                     ::  TMP_Drv
INTEGER                                     ::  rd, tpd, d
INTEGER                                     ::  Here
INTEGER                                     ::  iRE



#ifdef POSEIDON_AMREX_FLAG
    iRE = FEM_Elem_Map(iE(1),Level)
#else
    iRE = iE(1)
#endif




DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS

   
    TMP_Drv = 0.0_idp
    DO d = 0,DEGREE
        Here = Map_To_FEM_Node(iRE,d)


        TMP_Drv(1) = TMP_Drv(1)                                 &
                   + SUM( dVA_Coeff_Vector(Here,:,iU)           &
                        * Slm_Elem_Values(:,tpd)      )         &
                   * Lagrange_Poly_Table(d,rd,1)                &
                   / DROT

        TMP_Drv(2) = TMP_Drv(2)                                 &
                   + SUM( dVA_Coeff_Vector(Here,:,iU)           &
                        * Slm_Elem_dt_Values(:,tpd)   )         &
                   * Lagrange_Poly_Table(d,rd,0)

        TMP_Drv(3) = TMP_Drv(3)                                 &
                   + SUM( dVA_Coeff_Vector(Here,:,iU)           &
                        * Slm_Elem_dp_Values(:,tpd)   )         &
                   * Lagrange_Poly_Table(d,rd,0)



   END DO  ! d

   Drv(tpd,rd,1)       = TMP_Drv(1)
   Drv(tpd,rd,2)       = TMP_Drv(2)
   Drv(tpd,rd,3)       = TMP_Drv(2)

END DO ! td
END DO ! rd




END SUBROUTINE Calc_Drv_On_Elem_TypeA







 !+203+####################################################!
!                                                           !
!          Calc_Val_And_Drv_On_Elem_TypeA                   !
!                                                           !
 !#########################################################!
SUBROUTINE Calc_Val_And_Drv_On_Elem_TypeA( iE, DROT, Val, Drv, iU, Level)

INTEGER,   INTENT(IN), DIMENSION(3)         ::  iE
REAL(idp), INTENT(IN)                       ::  DROT

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points),    INTENT(OUT)  :: Val
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3),  INTENT(OUT)  :: Drv
INTEGER,   INTENT(IN)                       ::  iU
INTEGER,   INTENT(IN)                       ::  Level

REAL(idp)                                   ::  TMP_Val
REAL(idp), DIMENSION(3)                     ::  TMP_Drv

INTEGER                                     ::  rd, tpd, d
INTEGER                                     ::  Here
INTEGER                                     ::  iRE


#ifdef POSEIDON_AMREX_FLAG
    iRE = FEM_Elem_Map(iE(1),Level)
#else
    iRE = iE(1)
#endif


DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS

   TMP_Val = 0.0_idp
   TMP_Drv = 0.0_idp
   DO d = 0,DEGREE
       Here = Map_To_FEM_Node(iRE,d)

           
        TMP_Val = TMP_Val                                       &
                + SUM( dVA_Coeff_Vector(Here,:,iU)              &
                        * Slm_Elem_Values(:,tpd)      )         &
                * Lagrange_Poly_Table(d,rd,0)

        TMP_Drv(1) = TMP_Drv(1)                                 &
                   + SUM( dVA_Coeff_Vector(Here,:,iU)           &
                        * Slm_Elem_Values(:,tpd)      )         &
                   * Lagrange_Poly_Table(d,rd,1)                &
                   / DROT

        TMP_Drv(2) = TMP_Drv(2)                                 &
                   + SUM( dVA_Coeff_Vector(Here,:,iU)           &
                        * Slm_Elem_dt_Values(:,tpd)   )         &
                   * Lagrange_Poly_Table(d,rd,0)

        TMP_Drv(3) = TMP_Drv(3)                                 &
                   + SUM( dVA_Coeff_Vector(Here,:,iU)           &
                        * Slm_Elem_dp_Values(:,tpd)   )         &
                   * Lagrange_Poly_Table(d,rd,0)
                  



   END DO  ! d

   Val(tpd,rd)         = TMP_Val
   Drv(tpd,rd,1)       = TMP_Drv(1)
   Drv(tpd,rd,2)       = TMP_Drv(2)
   Drv(tpd,rd,3)       = TMP_Drv(2)

END DO ! td
END DO ! rd






END SUBROUTINE Calc_Val_And_Drv_On_Elem_TypeA















 !+301+####################################################!
!                                                           !
!          Calc_Val_On_Elem_TypeB                           !
!                                                           !
 !#########################################################!
SUBROUTINE Calc_Val_On_Elem_TypeB( iE, Val, iU, iVB, Level )

INTEGER,    INTENT(IN), DIMENSION(3)        ::  iE

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points),   INTENT(OUT)  :: Val
INTEGER,   INTENT(IN)                       ::  iU, iVB
INTEGER,    INTENT(IN)                      ::  Level

REAL(idp)                                   ::  TMP_Val
INTEGER                                     ::  rd, tpd, d
INTEGER                                     ::  Here, There
INTEGER                                     ::  iRE



#ifdef POSEIDON_AMREX_FLAG
    iRE = FEM_Elem_Map(iE(1),Level)
#else
    iRE = iE(1)
#endif




DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS

   
    TMP_Val = 0.0_idp
    DO d = 0,DEGREE

        Here  = FP_Array_Map_TypeB(iU,iVB,iRE,d,1)
        There = FP_Array_Map_TypeB(iU,iVB,iRE,d,LM_Length)


        TMP_Val = TMP_Val                                   &
                + SUM( dVB_Coeff_Vector(Here:There, iVB)    &
                        * Slm_Elem_Values(:,tpd)      )     &
                * Lagrange_Poly_Table(d,rd,0)


    END DO  ! d

    Val(tpd,rd)         = TMP_Val


END DO ! td
END DO ! rd




END SUBROUTINE Calc_Val_On_Elem_TypeB




 !+302+####################################################!
!                                                           !
!          Calc_Drv_On_Elem_TypeB                           !
!                                                           !
 !#########################################################!
SUBROUTINE Calc_Drv_On_Elem_TypeB( iE, DROT, Drv, iU, iVB, Level )

INTEGER,    INTENT(IN), DIMENSION(3)        ::  iE
REAL(idp), INTENT(IN)                       :: DROT

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3),  INTENT(OUT)  :: Drv
INTEGER,   INTENT(IN)                       :: iU, iVB
INTEGER,    INTENT(IN)                      :: Level


REAL(idp), DIMENSION(3)                     ::  TMP_Drv
INTEGER                                     ::  rd, tpd, d
INTEGER                                     ::  Here, There
INTEGER                                     ::  iRE



#ifdef POSEIDON_AMREX_FLAG
    iRE = FEM_Elem_Map(iE(1),Level)
#else
    iRE = iE(1)
#endif




DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS

   
    TMP_Drv = 0.0_idp
    DO d = 0,DEGREE
        Here  = FP_Array_Map_TypeB(iU,iVB,iRE,d,1)
        There = FP_Array_Map_TypeB(iU,iVB,iRE,d,LM_Length)


        TMP_Drv(1) = TMP_Drv(1)                                 &
                   + SUM( dVB_Coeff_Vector(Here:There, iVB)     &
                        * Slm_Elem_Values(:,tpd)      )         &
                   * Lagrange_Poly_Table(d,rd,1)                &
                   / DROT
                    
        TMP_Drv(2) = TMP_Drv(2)                                &
                   + SUM( dVB_Coeff_Vector(Here:There, iVB)    &
                        * Slm_Elem_dt_Values(:,tpd)   )        &
                   * Lagrange_Poly_Table(d,rd,0)

        TMP_Drv(3) = TMP_Drv(3)                                &
                   + SUM( dVB_Coeff_Vector(Here:There, iVB)    &
                        * Slm_Elem_dp_Values(:,tpd)   )        &
                   * Lagrange_Poly_Table(d,rd,0)


    END DO  ! d

   Drv(tpd,rd,1)       = TMP_Drv(1)
   Drv(tpd,rd,2)       = TMP_Drv(2)
   Drv(tpd,rd,3)       = TMP_Drv(2)
   
END DO ! td
END DO ! rd




END SUBROUTINE Calc_Drv_On_Elem_TypeB



 !+301+####################################################!
!                                                           !
!          Calc_Val_And_Drv_On_Elem_TypeB                   !
!                                                           !
 !#########################################################!
SUBROUTINE Calc_Val_And_Drv_On_Elem_TypeB( iE, DROT, Val, Drv, iU, iVB, Level )

INTEGER,    INTENT(IN), DIMENSION(3)        ::  iE
REAL(idp), INTENT(IN)                       ::  DROT

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points),    INTENT(OUT)  :: Val
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3),  INTENT(OUT)  :: Drv
INTEGER,   INTENT(IN)                       ::  iU, iVB
INTEGER,    INTENT(IN)                      ::  Level



REAL(idp)                                   ::  TMP_Val
REAL(idp), DIMENSION(3)                     ::  TMP_Drv

INTEGER                                     ::  rd, tpd, d
INTEGER                                     ::  Here, There
INTEGER                                     ::  iRE



#ifdef POSEIDON_AMREX_FLAG
    iRE = FEM_Elem_Map(iE(1),Level)
#else
    iRE = iE(1)
#endif



DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS


    TMP_Val = 0.0_idp
    TMP_Drv = 0.0_idp
    DO d = 0,DEGREE
        Here  = FP_Array_Map_TypeB(iU,iVB,iRE,d,1)
        There = FP_Array_Map_TypeB(iU,iVB,iRE,d,LM_Length)

        TMP_Val = TMP_Val                                       &
               + SUM( dVB_Coeff_Vector( Here:There, iVB )      &
                       * Slm_Elem_Values( :, tpd )   )       &
               * Lagrange_Poly_Table( d, rd, 0 )



        TMP_Drv(1) = TMP_Drv(1)                                 &
                  + SUM( dVB_Coeff_Vector( Here:There, iVB )   &
                        * Slm_Elem_Values( :, tpd )     )    &
                  * Lagrange_Poly_Table( d, rd, 1 )             &
                  / DROT

        TMP_Drv(2) = TMP_Drv(2)                                 &
                  + SUM( dVB_Coeff_Vector( Here:There, iVB )   &
                        * Slm_Elem_dt_Values( :, tpd )   )    &
                  * Lagrange_Poly_Table( d, rd, 0)



        TMP_Drv(3) = TMP_Drv(3)                                 &
                  + SUM( dVB_Coeff_Vector( Here:There, iVB )   &
                        * Slm_Elem_dp_Values( :, tpd )   )    &
                  * Lagrange_Poly_Table( d, rd, 0)


   END DO  ! d

   Val(tpd,rd)         = TMP_Val
   Drv(tpd,rd,1)       = TMP_Drv(1)
   Drv(tpd,rd,2)       = TMP_Drv(2)
   Drv(tpd,rd,3)       = TMP_Drv(2)

END DO ! td
END DO ! rd







END SUBROUTINE Calc_Val_And_Drv_On_Elem_TypeB
























END MODULE Load_Vector_Functions_Calc_Values_Module
