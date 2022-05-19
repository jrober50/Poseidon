   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Initialization_Subroutines                                            !##!
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

USE Poseidon_Parameters, &
            ONLY :  Degree,                         &
                    L_Limit,                        &
                    Verbose_Flag,                   &
                    Convergence_Criteria,           &
                    Convergence_Criteria_Default,   &
                    Max_Iterations,                 &
                    Max_Iterations_Default

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    Num_T_Elements,         &
                    Num_P_Elements

USE Variables_IO, &
            ONLY :  Report_Flags,           &
                    iRF_Setup,              &
                    iRF_Time,               &
                    iWF_Source,             &
                    iWF_Results,            &
                    Write_Flags,            &
                    File_Suffix

USE Variables_FP, &
            ONLY :  FP_Anderson_M,          &
                    FP_Anderson_M_Default

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Max_Grid_Size,    &
                    AMReX_Max_Level,        &
                    AMReX_Num_Levels,       &
                    iFRL,                   &
                    iIRL

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Num_T_Quad_Points,      &
                    Num_P_Quad_Points,      &
                    Num_TP_Quad_Points,     &
                    Local_Quad_DOF,         &
                    Int_R_Locations,        &
                    Int_T_Locations,        &
                    Int_P_Locations,        &
                    Int_R_Weights,          &
                    Int_T_Weights,          &
                    Int_P_Weights,          &
                    Int_TP_Weights


USE Variables_Interface, &
            ONLY :  Caller_Set,                     &
                    Caller_nLevels,                 &
                    Caller_NQ,                      &
                    Caller_Quad_DOF,                &
                    Caller_xL,                      &
                    Caller_RQ_xlocs,                &
                    Caller_TQ_xlocs,                &
                    Caller_PQ_xlocs,                &
                    Translation_Matrix

USE Functions_Math, &
            ONLY :  Lagrange_Poly

USE Maps_X_Space, &
            ONLY :  Map_To_X_Space

IMPLICIT NONE


CONTAINS



 !+101+############################################################################!
!                                                                                   !
!       Init_IO_Params                                                              !
!                                                                                   !
 !#################################################################################!
SUBROUTINE Init_IO_Params(  WriteAll_Option,            &
                            Print_Setup_Option,         &
                            Write_Setup_Option,         &
                            Print_Results_Option,       &
                            Write_Results_Option,       &
                            Print_Timetable_Option,     &
                            Write_Timetable_Option,     &
                            Write_Sources_Option,       &
                            Suffix_Flag_Option,         &
                            Suffix_Tail_Option,         &
                            Frame_Option                )


LOGICAL,            INTENT(IN), OPTIONAL               ::  WriteAll_Option
LOGICAL,            INTENT(IN), OPTIONAL               ::  Print_Setup_Option
LOGICAL,            INTENT(IN), OPTIONAL               ::  Write_Setup_Option
LOGICAL,            INTENT(IN), OPTIONAL               ::  Print_Results_Option
LOGICAL,            INTENT(IN), OPTIONAL               ::  Write_Results_Option
LOGICAL,            INTENT(IN), OPTIONAL               ::  Print_Timetable_Option
LOGICAL,            INTENT(IN), OPTIONAL               ::  Write_Timetable_Option
LOGICAL,            INTENT(IN), OPTIONAL               ::  Write_Sources_Option

CHARACTER(LEN=10),  INTENT(IN), OPTIONAL               ::  Suffix_Flag_Option
CHARACTER(LEN=1),   INTENT(IN), OPTIONAL               ::  Suffix_Tail_Option
INTEGER,            INTENT(IN), OPTIONAL               ::  Frame_Option



IF ( Verbose_Flag ) THEN
    Report_Flags = 1
    Write_Flags  = 1
ELSE
    Report_Flags = 0
    Write_Flags  = 0
END IF

IF ( PRESENT( WriteAll_Option) ) THEN
    Report_Flags = Report_Flags + 2
    Write_Flags = Write_Flags + 2
END IF




IF ( PRESENT(Print_Setup_Option) ) THEN
    IF ( Print_Setup_Option ) THEN
        IF (Report_Flags(iRF_Setup) > 1 ) THEN
            Report_Flags(iRF_Setup) = 3
        ELSE
            Report_Flags(iRF_Setup) = 1
        END IF
    ELSE
        IF (Report_Flags(iRF_Setup) > 1 ) THEN
            Report_Flags(iRF_Setup) = 2
        ELSE
            Report_Flags(iRF_Setup) = 0
        END IF
    END IF
END IF


IF ( PRESENT(Write_Setup_Option) ) THEN
    IF ( Write_Setup_Option ) THEN
        IF (Report_Flags(iRF_Setup) < 2 ) THEN
            Report_Flags(iRF_Setup) = Report_Flags(4) + 2
        END IF
    ELSE
        IF ( Report_Flags(iRF_Setup) > 1 ) THEN
            Report_Flags(iRF_Setup) = Report_Flags(4) - 2
        END IF
    END IF
END IF





IF ( PRESENT(Print_Results_Option) ) THEN
    IF ( Print_Results_Option ) THEN
        IF (Write_Flags(iWF_Results) > 1 ) THEN
            Write_Flags(iWF_Results) = 3
        ELSE
            Write_Flags(iWF_Results) = 1
        END IF
    ELSE
        IF (Write_Flags(iWF_Results) > 1 ) THEN
            Write_Flags(iWF_Results) = 2
        ELSE
            Write_Flags(iWF_Results) = 0
        END IF
    END IF
END IF


IF ( PRESENT(Write_Results_Option) ) THEN
    IF ( Write_Results_Option ) THEN
        IF (Write_Flags(iWF_Results) < 2 ) THEN
            Write_Flags(iWF_Results) = Write_Flags(5) + 2
        END IF
    ELSE
        IF ( Write_Flags(iWF_Results) > 1 ) THEN
            Write_Flags(iWF_Results) = Write_Flags(5) - 2
        END IF
    END IF
END IF







IF ( PRESENT(Print_Timetable_Option) ) THEN
    IF ( Print_Timetable_Option ) THEN
        IF (Report_Flags(iRF_Time) > 1 ) THEN
            Report_Flags(iRF_Time) = 3
        ELSE
            Report_Flags(iRF_Time) = 1
        END IF
    ELSE
        IF (Report_Flags(iRF_Time) > 1 ) THEN
            Report_Flags(iRF_Time) = 2
        ELSE
            Report_Flags(iRF_Time) = 0
        END IF
    END IF
END IF



IF ( PRESENT(Write_Timetable_Option) ) THEN
    IF ( Write_Timetable_Option ) THEN
        IF (Report_Flags(iRF_Time) < 2 ) THEN
            Report_Flags(iRF_Time) = Report_Flags(iRF_Time) + 2
        END IF
    ELSE
        IF ( Report_Flags(iRF_Time) > 1 ) THEN
            Report_Flags(iRF_Time) = Report_Flags(iRF_Time) - 2
        END IF
    END IF
END IF



IF ( PRESENT(Write_Sources_Option) ) THEN
    IF ( Write_Sources_Option ) THEN
        Write_Flags(iWF_Source) = 2
    ELSE
        Write_Flags(iWF_Source) = 0
    END IF
END IF


IF ( PRESENT(Suffix_Flag_Option) ) THEN


    IF ( Suffix_Flag_Option == "Params") THEN

        WRITE(File_Suffix,'(A,I5.5,A,I3.3,A,I2.2,A,I2.2)')             &
            "RE",Num_R_Elements,"_TE",Num_T_Elements,"_D",Degree,"_L",L_Limit

    ELSEIF ( SUffix_Flag_Option == "Frame") THEN
        
        IF ( PRESENT(Frame_Option) ) THEN
            WRITE(File_Suffix,'(I5.5)') Frame_Option
        ELSE
            WRITE(File_Suffix,'(I5.5)') 1
        END IF
    END IF
ELSE
    WRITE(File_Suffix,'(I5.5)') 1
END IF


IF ( PRESENT(Suffix_Tail_Option) ) THEN
    WRITE(File_Suffix,'(A,A,A)') TRIM(File_Suffix),"_",Suffix_Tail_Option

END IF




END SUBROUTINE Init_IO_Params








 !+102+############################################################################!
!                                                                                   !
!       Init_Fixed_Point_Params                                                     !
!                                                                                   !
 !#################################################################################!
SUBROUTINE Init_Fixed_Point_Params( Max_Iterations_Option,          &
                                    Convergence_Criteria_Option,    &
                                    Anderson_M_Option               )

INTEGER,                 INTENT(IN), OPTIONAL               ::  Max_Iterations_Option
REAL(idp),               INTENT(IN), OPTIONAL               ::  Convergence_Criteria_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  Anderson_M_Option





IF ( PRESENT( Max_Iterations_Option ) ) THEN
    Max_Iterations = Max_Iterations_Option
ELSE
    Max_Iterations = Max_Iterations_Default
END IF


IF ( PRESENT( Convergence_Criteria_Option) ) THEN
    Convergence_Criteria = Convergence_Criteria_Option
ELSE
    Convergence_Criteria = Convergence_Criteria_Default
END IF


IF ( PRESENT(Anderson_M_Option) ) THEN
    FP_Anderson_M = Anderson_M_Option
ELSE
    FP_Anderson_M = FP_Anderson_M_Default
END IF


END SUBROUTINE Init_Fixed_Point_Params





 !+102+####################################################################!
!                                                                           !
!       Init_AMReX_Params                                                   !
!                                                                           !
 !#########################################################################!
SUBROUTINE Init_AMReX_Params(   AMReX_Max_Level_Option,             &
                                AMReX_Max_Grid_Size_Option,         &
                                AMReX_FEM_Refinement_Option,        &
                                AMReX_Integral_Refinement_Option    )


INTEGER,                 INTENT(IN), OPTIONAL               ::  AMReX_Max_Level_Option
INTEGER,   DIMENSION(3), INTENT(IN), OPTIONAL               ::  AMReX_Max_Grid_Size_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  AMReX_FEM_Refinement_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  AMReX_Integral_Refinement_Option


IF ( PRESENT(AMReX_Max_Level_Option) ) THEN
    AMReX_Max_Level  = AMReX_Max_Level_Option
ELSE
    AMReX_Max_Level  = 0
END IF
AMReX_Num_Levels = AMReX_Max_Level+1




IF ( PRESENT(AMReX_Max_Level_Option) ) THEN
    AMReX_Max_Grid_Size = AMReX_Max_Grid_Size_Option
ELSE
    AMReX_Max_Grid_Size = 4
END IF



IF ( PRESENT(AMReX_FEM_Refinement_Option) ) THEN
    iFRL = AMReX_FEM_Refinement_Option
ELSE
    iFRL = AMReX_Max_Level
END IF

IF ( PRESENT(AMReX_Integral_Refinement_Option) ) THEN
    iIRL = AMReX_Integral_Refinement_Option
ELSE
    iIRL = AMReX_Max_Level
END IF


END SUBROUTINE Init_AMReX_Params










 !+102+####################################################################!
!                                                                           !
!       Init_AMReX_Params                                                   !
!                                                                           !
 !#########################################################################!
SUBROUTINE Set_Caller_Quadrature( Source_NQ,        &
                                  Source_xL,        &
                                  Source_RQ_xlocs,  &
                                  Source_TQ_xlocs,  &
                                  Source_PQ_xlocs   )


INTEGER,    DIMENSION(3),               INTENT(IN)      ::  Source_NQ
REAL(idp),  DIMENSION(2),               INTENT(IN)      ::  Source_xL
REAL(idp),  DIMENSION(1:Source_NQ(1)),  INTENT(IN)      ::  Source_RQ_xlocs
REAL(idp),  DIMENSION(1:Source_NQ(2)),  INTENT(IN)      ::  Source_TQ_xlocs
REAL(idp),  DIMENSION(1:Source_NQ(3)),  INTENT(IN)      ::  Source_PQ_xlocs


INTEGER                                             ::  Local_R
INTEGER                                             ::  Local_T
INTEGER                                             ::  Local_P

INTEGER                                             ::  Caller_T
INTEGER                                             ::  Caller_P

INTEGER                                             ::  Here
INTEGER                                             ::  There
INTEGER                                             ::  Local_Here

REAL(idp),  DIMENSION( 1:Source_NQ(1) )              ::  Scaled_R_Quad
REAL(idp),  DIMENSION( 1:Source_NQ(2) )              ::  Scaled_T_Quad
REAL(idp),  DIMENSION( 1:Source_NQ(3) )              ::  Scaled_P_Quad

REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  R_Lag_Poly_Values
REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  T_Lag_Poly_Values
REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  P_Lag_Poly_Values


Caller_nLevels = AMReX_Num_Levels
Caller_NQ = Source_NQ
Caller_xL = Source_xL

ALLOCATE( Caller_RQ_xlocs(1:Caller_NQ(1)) )
ALLOCATE( Caller_TQ_xlocs(1:Caller_NQ(2)) )
ALLOCATE( Caller_PQ_xlocs(1:Caller_NQ(3)) )

Caller_RQ_xlocs = Source_RQ_xlocs
Caller_TQ_xlocs = Source_TQ_xlocs
Caller_PQ_xlocs = Source_PQ_xlocs

Caller_Quad_DOF = Caller_NQ(1)*Caller_NQ(2)*Caller_NQ(3)

Caller_Set = .TRUE.


PRINT*,Caller_Quad_DOF, Local_Quad_DOF

ALLOCATE( Translation_Matrix(1:Caller_Quad_DOF, 1:Local_Quad_DOF) )
ALLOCATE( R_Lag_Poly_Values(1:Caller_NQ(1),1:NUM_R_QUAD_POINTS) )
ALLOCATE( T_Lag_Poly_Values(1:Caller_NQ(2),1:NUM_T_QUAD_POINTS) )
ALLOCATE( P_Lag_Poly_Values(1:Caller_NQ(3),1:NUM_P_QUAD_POINTS) )


Scaled_R_Quad = Map_To_X_Space(Caller_xL(1),Caller_xL(2),Caller_RQ_xlocs)
Scaled_T_Quad = Map_To_X_Space(Caller_xL(1),Caller_xL(2),Caller_TQ_xlocs)
Scaled_P_Quad = Map_To_X_Space(Caller_xL(1),Caller_xL(2),Caller_PQ_xlocs)




DO Local_R = 1,NUM_R_QUAD_POINTS
    R_Lag_Poly_Values(:,Local_R) = Lagrange_Poly( Int_R_Locations(Local_R), &
                                                  Caller_NQ(1)-1,            &
                                                  Scaled_R_Quad              )

END DO

DO Local_T = 1,NUM_T_QUAD_POINTS
    T_Lag_Poly_Values(:,Local_T) = Lagrange_Poly( Int_T_Locations(Local_T), &
                                                  Caller_NQ(2)-1,            &
                                                  Scaled_T_Quad              )

END DO

DO Local_P = 1,NUM_P_QUAD_POINTS
    P_Lag_Poly_Values(:,Local_P) = Lagrange_Poly( Int_P_Locations(Local_P), &
                                                  Caller_NQ(3)-1,            &
                                                  Scaled_P_Quad              )

END DO



DO Local_P = 1,NUM_P_QUAD_POINTS
DO Local_T = 1,NUM_T_QUAD_POINTS
DO Local_R = 1,NUM_R_QUAD_POINTS

    Local_Here = (Local_P-1) * NUM_T_QUAD_POINTS * NUM_R_QUAD_POINTS        &
               + (Local_T-1) * NUM_R_QUAD_POINTS                            &
               + Local_R

    DO Caller_P = 1,Caller_NQ(3)
    DO Caller_T = 1,Caller_NQ(2)

            Here = (Caller_P-1) * Caller_NQ(2) * Caller_NQ(1)   &
                 + (Caller_T-1) * Caller_NQ(1)

            There = Here + Caller_NQ(1)

            Translation_Matrix(Here+1:There, Local_Here)  =                 &
                              R_Lag_Poly_Values(1:Caller_NQ(1),Local_R)    &
                            * T_Lag_Poly_Values(Caller_T,Local_T)            &
                            * P_Lag_Poly_Values(Caller_P,Local_P)

    END DO  !   Caller_T Loop
    END DO  !   Caller_P Loop
END DO  !   Local_R Loop
END DO  !   Local_T Loop
END DO  !   Local_P Loop



DEALLOCATE( R_Lag_Poly_Values )
DEALLOCATE( T_Lag_Poly_Values )
DEALLOCATE( P_Lag_Poly_Values )

END SUBROUTINE Set_Caller_Quadrature




END MODULE Initialization_Subroutines
