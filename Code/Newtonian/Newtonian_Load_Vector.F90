   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poisson_Load_Vector                                                 !##!
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

USE Poseidon_Message_Routines_Module, &
            ONLY :  Run_Message

USE Poseidon_Numbers_Module, &
            ONLY :  pi,                         &
                    TwoPi

USE Poseidon_Units_Module, &
            ONLY :  Grav_Constant_G

USE Poseidon_Parameters, &
            ONLY :  Degree,             &
                    L_Limit,            &
                    Verbose_Flag

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,             &
                    Num_T_Elements,             &
                    Num_P_Elements,             &
                    rlocs,                      &
                    tlocs,                      &
                    plocs

USE Variables_Derived, &
            ONLY :  LM_LENGTH,                  &
                    Num_R_Nodes

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

USE Variables_Vectors, &
            ONLY :  cVA_Load_Vector

USE Variables_Source, &
            ONLY :  Source_Rho

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature,       &
                    Initialize_LGL_Quadrature,      &
                    Initialize_Trapezoid_Quadrature

USE Functions_Math, &
            ONLY :  Legendre_Poly,                  &
                    Lagrange_Poly,                  &
                    Norm_Factor

USE Maps_X_Space, &
            ONLY :  Map_From_X_Space

USE Maps_Quadrature, &
            ONLY :  Quad_Map

USE Maps_Domain, &
            ONLY :  Map_To_FEM_Node,                &
                    Map_To_lm

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Poisson_SourceVector_Subparts,          &
                    Timer_Poisson_SourceVector_Main



IMPLICIT NONE


CONTAINS



 !+201b+####################################################################!
!                                                                           !
!                      Calculate_Poisson_Load_Vector                      !
!                                                                           !
!---------------------------------------------------------------------------!
!                                                                           !
!   This function performs the triple intergral over the source function    !
!       used to defined values in the source vector of the linear system    !
!       established by the mixed spectral/fintie element method.            !
!                                                                           !
 !#########################################################################!
SUBROUTINE Calculate_Poisson_Load_Vector()

COMPLEX(idp)                                        ::  Tmp_Val


INTEGER                                             ::  lm,l,m, re, rd, te, td, pe, pd, d


REAL(idp)                                           ::  drot, dtot
REAL(idp)                                           ::  SphereHarm_NormFactor

REAL(idp), ALLOCATABLE, DIMENSION(:)                ::  R_locs, P_locs, T_locs
                                                        


REAL(KIND = idp), DIMENSION(0:DEGREE)               ::  Poly_xlocs
REAL(KIND = idp), DIMENSION(0:DEGREE)               ::  Poly_weights
REAL(KIND = idp), DIMENSION(0:DEGREE)               ::  LagP

INTEGER                                             ::  There, Here


REAL(idp),           DIMENSION(:,:,:), ALLOCATABLE  ::  R_Pre
REAL(idp),           DIMENSION(:,:,:), ALLOCATABLE  ::  T_Pre
COMPLEX(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE  ::  P_Pre

ALLOCATE(R_locs(1:Num_R_Quad_Points ) )
ALLOCATE(T_locs(1:Num_T_Quad_Points ) )
ALLOCATE(P_locs(1:Num_P_Quad_Points ) )



ALLOCATE( R_Pre(1:Num_R_Quad_Points,0:Degree,0:Num_R_Elements-1) )
ALLOCATE( T_Pre(1:Num_T_Quad_Points,0:Num_T_Elements-1,1:LM_Length) )
ALLOCATE( P_Pre(1:Num_P_Quad_Points,0:Num_P_Elements-1,1:LM_Length) )


IF ( Verbose_Flag ) CALL Run_Message("Calculating Load Vector.")



CALL Initialize_LGL_Quadrature(DEGREE, Poly_xlocs, Poly_weights)


P_locs = Map_From_X_Space(0.0_idp, TwoPi, Int_P_Locations)
T_locs = Map_From_X_Space(0.0_idp, Pi, Int_T_Locations)


!#if defined(POSEIDON_OPENMP_FLAG)
!!$OMP PARALLEL
!PRINT*,"Hello From Process ",OMP_GET_THREAD_NUM()
!!$OMP END PARALLEL
!#endif

#if defined(POSEIDON_OPENMP_OL_FLAG)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:     ) &
    !$OMP MAP( alloc: iErr )
#elif defined(POSEIDON_OPENACC_FLAG)
    !$ACC ENTER DATA &
    !$ACC COPYIN(      ) &
    !$ACC CREATE(     iErr )
#endif

CALL TimerStart( Timer_Poisson_SourceVector_Subparts )

DO re = 0, Num_R_Elements -1
DO  d = 0, Degree

    There = Map_To_FEM_Node(re,d)
    R_locs = Map_From_X_Space(rlocs(re), rlocs(re+1), Int_R_Locations)
    drot = 0.5_idp *(rlocs(re+1) - rlocs(re))

    DO rd = 1, Num_R_Quad_Points

        LagP = Lagrange_Poly(Int_R_Locations(rd), DEGREE, Poly_xlocs)

        R_Pre(rd,d,re) = 4.0_idp * pi          &
                        * Grav_Constant_G       &
                        * R_locs(rd)            &
                        * R_locs(rd)            &
                        * drot                  &
                        * Int_R_weights(rd)     &
                        * LagP(d)


    END DO
END DO
END DO




DO l = 0, L_LIMIT
DO m = -l,l
DO te = 0, NUM_T_ELEMENTS - 1

    dtot = 0.5_idp *(tlocs(te+1) - tlocs(te))
    T_locs = Map_From_X_Space(tlocs(te), tlocs(te+1), Int_T_Locations)
    SphereHarm_NormFactor = (-1.0_idp**m)*Norm_Factor(l,-m)
    lm = Map_To_lm( l, m )

    T_Pre(:,te,lm) = sin(T_locs(:))                                     &
                    * SphereHarm_NormFactor                             &
                    * Legendre_Poly(l, -m, Num_T_Quad_Points, T_locs)   &
                    * dtot                                              &
                    * Int_T_weights(:)

END DO ! te
END DO ! lm
END DO



DO l = 0, L_LIMIT
DO m = -l,l
DO pe = 0, Num_P_Elements - 1
DO pd = 1,Num_P_Quad_Points
    lm = Map_To_lm( l, m )
    P_locs = Map_From_X_Space(plocs(pe), plocs(pe+1), Int_P_Locations)
    P_PRE(pd,pe,lm) = CDEXP(CMPLX(0.0_idp,-m * P_locs(pd),idp)) * Int_P_weights(pd)
END DO
END DO
END DO
END DO



CALL TimerStop( Timer_Poisson_SourceVector_Subparts)
CALL TimerStart( Timer_Poisson_SourceVector_Main )






cVA_Load_Vector = 0.0_idp
Tmp_Val = 0.0_idp

#if defined(POSEIDON_OPENMP_OL_FLAG)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE(  ) &
    !$OMP REDUCTION( MIN: TimeStep )
#elif defined(POSEIDON_OPENACC_FLAG)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE(  ) &
    !$ACC PRESENT( ) &
    !$ACC REDUCTION( MIN: TimeStep )
#elif defined(POSEIDON_OPENMP_FLAG)
    !$OMP PARALLEL DO SIMD COLLAPSE(3) REDUCTION(+:Source_Vector) &
    !$OMP PRIVATE(  re,te,pe,lm,p,              &
    !$OMP           Here, There, TMP_Val, Int_Term   )
#endif


DO lm = 1, LM_Length
DO re = 0,NUM_R_ELEMENTS - 1
DO d = 0,DEGREE

    
    There = Map_To_FEM_Node( re, d )
    
    DO pe = 0,NUM_P_ELEMENTS - 1
    DO te = 0,NUM_T_ELEMENTS - 1

    DO pd = 1, Num_P_Quad_Points
    DO td = 1, Num_T_Quad_Points
    DO rd = 1, Num_R_Quad_Points


        Here = Quad_Map(rd,td,pd)


        Tmp_Val = Tmp_Val                       &
                - Source_Rho(Here,re,te,pe)   &
                * R_PRE(rd,d,re)               &
                * T_PRE(td,te,lm)               &
                * P_PRE(pd,pe,lm)

    END DO ! rd Loop
    END DO ! td Loop
    END DO ! pd Loop

    END DO  ! te Loop
    END DO  ! pe Loop
    
    cVA_Load_Vector(There,lm,1) = cVA_Load_Vector(There,lm,1) + Tmp_Val

    TMP_Val = 0.0_idp

END DO  ! d Loop
END DO  ! re Loop
END DO  ! lm Loop




#if defined(POSEIDON_OPENMP_OL_FLAG)
    !$OMP END PARALLEL DO SIMD
#elif defined(POSEIDON_OPENACC_FLAG)
    !$ACC END PARALLEL LOOP
#elif defined(POSEIDON_OPENMP_FLAG)
    !$OMP END PARALLEL DO SIMD
#endif

CALL TimerStop( Timer_Poisson_SourceVector_Main )


END SUBROUTINE Calculate_Poisson_Load_Vector









END MODULE Poisson_Load_Vector
