   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE IO_Setup_Report_Module                                                        !##!
!##!                                                                                !##!
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

USE Poseidon_Numbers_Module, &
                ONLY :  pi

USE Units_Module, &
                ONLY :  Set_Units

USE Poseidon_Parameters, &
                ONLY :  Domain_Dim,             &
                        Degree,                 &
                        L_Limit,                &
                        Method_Flag,            &
                        Verbose_Flag,           &
                        Convergence_Criteria,   &
                        Num_CFA_Vars,           &
                        Max_Iterations

USE Poseidon_IO_Parameters, &
                ONLY :  Method_Names,           &
                        Poseidon_Reports_Dir

USE Variables_MPI, &
                ONLY :  myID_Poseidon,          &
                        MasterID_Poseidon

USE Variables_Mesh, &
                ONLY :  Num_R_Elements,         &
                        Num_T_Elements,         &
                        Num_P_Elements

USE Variables_Quadrature, &
                ONLY :  Num_R_Quad_Points,      &
                        Num_T_Quad_Points,      &
                        Num_P_Quad_Points

USE Variables_FP,   &
                ONLY :  FP_Anderson_M

USE Variables_IO, &
                ONLY :  Report_Flags,           &
                        Report_IDs,             &
                        File_Suffix

USE IO_File_Routines_Module, &
                ONLY :  Open_New_File

IMPLICIT NONE

PUBLIC :: Output_Setup_Report


CONTAINS


!+101+######################################################################################!
!                                                                                           !
!       Output_Setup_Report                                                                 !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Output_Setup_Report()

CHARACTER(LEN = 100)                                ::  Report_Name
INTEGER                                             ::  Suggested_Number = 400

PRINT*, "In Output_Setup_Report",myID_Poseidon
IF (myID_Poseidon == MasterID_Poseidon ) THEN


    IF ( Report_Flags(4) > 1 ) THEN
        WRITE(Report_Name,'(A,A,A,A)') Poseidon_Reports_Dir,"Setup_Report_",trim(File_Suffix),".out"
        CALL Open_New_File( Report_Name, Report_IDs(4), Suggested_Number)
    END IF

    CALL Output_Params_Report( Report_IDs(4) )
    CALL Output_Method_Report( Report_IDs(4) )


END IF


END SUBROUTINE Output_Setup_Report





!+201+##################################################################!
!                                                                       !
!       OUTPUT_SETUP_TABLE                                              !
!                                                                       !
!#######################################################################!
SUBROUTINE Output_Params_Report( Report_ID )

INTEGER, INTENT(IN)                                 ::  Report_ID


1400 FORMAT(/)
1401 FORMAT('------------- POSEIDON PARAMETERS --------------')

1402 FORMAT(/'------------- Expansion Parameters ------------')
1403 FORMAT('      Finite Element Degree = ',I3.1)
1404 FORMAT(' Spherical Harmonic L Limit = ',I3.1)

1405 FORMAT(/'------------- Mesh Parameters -----------------')
1406 FORMAT(' # Radial Elements = ',I5.1)
1407 FORMAT(' # Theta Elements  = ',I5.1)
1408 FORMAT(' # Phi Elements    = ',I5.1)


1410 FORMAT(/'------------- Quadrature Parameters -----------')
1411 FORMAT(' # Radial Quad Points = ',I3.1)
1412 FORMAT(' # Theta Quad Points  = ',I3.1)
1413 FORMAT(' # Phi Quad Points    = ',I3.1)


IF ( (Report_Flags(4) == 1) .OR. (Report_Flags(4) == 3)) THEN
    WRITE(*,1400)
    WRITE(*,1401)
    WRITE(*,1402)
    WRITE(*,1403)Degree
    WRITE(*,1404)L_LIMIT
    WRITE(*,1405)
    WRITE(*,1406)Num_R_Elements
    WRITE(*,1407)Num_T_Elements
    WRITE(*,1408)Num_P_Elements
    WRITE(*,1410)
    WRITE(*,1411)Num_R_Quad_Points
    WRITE(*,1412)Num_T_Quad_Points
    WRITE(*,1413)Num_P_Quad_Points
    WRITE(*,1400)
END IF

IF ( Report_ID .NE. -1 ) THEN

    WRITE(Report_ID,1400)
    WRITE(Report_ID,1401)
    WRITE(Report_ID,1402)
    WRITE(Report_ID,1403)Degree
    WRITE(Report_ID,1404)L_LIMIT
    WRITE(Report_ID,1405)
    WRITE(Report_ID,1406)Num_R_Elements
    WRITE(Report_ID,1407)Num_T_Elements
    WRITE(Report_ID,1408)Num_P_Elements
    WRITE(Report_ID,1410)
    WRITE(Report_ID,1411)Num_R_Quad_Points
    WRITE(Report_ID,1412)Num_T_Quad_Points
    WRITE(Report_ID,1413)Num_P_Quad_Points
    WRITE(Report_ID,1400)


END IF



END SUBROUTINE Output_Params_Report









!+201+##################################################################!
!                                                                       !
!       OUTPUT_SETUP_TABLE                                              !
!                                                                       !
!#######################################################################!
SUBROUTINE Output_Method_Report( Report_ID )

INTEGER, INTENT(IN)                                 ::  Report_ID

1500 FORMAT(/)
1501 FORMAT('------------- Method Parameters --------------')
1502 FORMAT(' Solver Method : ',A /)
1503 FORMAT(' Maximum Iterations   = ',I4.0)
1504 FORMAT(' Convergence Criteria = ',ES12.4E3)

1505 FORMAT(' Anderson Acceleration Memory Parameter = ',I3.1)



IF ( (Report_Flags(4) == 1) .OR. (Report_Flags(4) == 3)) THEN
    WRITE(*,1501)
    WRITE(*,1502)Method_Names(Method_Flag)
    WRITE(*,1503)Max_Iterations
    WRITE(*,1504)Convergence_Criteria
    IF ( Method_Flag >= 2 ) THEN
        WRITE(*,1505) FP_Anderson_M
    END IF
    WRITE(*,1500)
END IF



IF ( Report_ID .NE. -1 ) THEN
    WRITE(Report_ID,1501)
    WRITE(Report_ID,1502)Method_Names(Method_Flag)
    WRITE(Report_ID,1503)Max_Iterations
    WRITE(Report_ID,1504)Convergence_Criteria
    IF ( Method_Flag >= 2 ) THEN
        WRITE(Report_ID,1505) FP_Anderson_M
    END IF
END IF


END SUBROUTINE Output_Method_Report








END MODULE IO_Setup_Report_Module
