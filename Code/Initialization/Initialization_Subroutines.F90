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










END MODULE Initialization_Subroutines
