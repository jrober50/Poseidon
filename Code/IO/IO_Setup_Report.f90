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

USE Poseidon_Units_Module, &
            ONLY :  Set_Units,              &
                    Centimeter

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
                    Num_P_Elements,         &
                    R_Inner,                &
                    R_Outer

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Num_T_Quad_Points,      &
                    Num_P_Quad_Points

USE Variables_FP,   &
            ONLY :  FP_Anderson_M

USE Variables_IO, &
            ONLY :  Report_IDs,             &
                    File_Suffix,            &
                    iRF_Setup

USE Poseidon_File_Routines_Module, &
            ONLY :  Open_New_File


USE Variables_AMReX_Core, &
            ONLY :  AMReX_Max_Level,        &
                    iNumLeafElements

USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,               &
                    iPF_IO_Print_Setup,         &
                    iPF_IO_Write_Setup

USE Flags_Core_Module, &
            ONLY :  lPF_Core_Flags,         &
                    iPF_Core_Newtonian_Mode

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


IF (myID_Poseidon == MasterID_Poseidon ) THEN


    IF ( lPF_IO_Flags(iPF_IO_Write_Setup) ) THEN
        WRITE(Report_Name,'(A,A,A,A)') Poseidon_Reports_Dir,"Setup_Report_",trim(File_Suffix),".out"
        CALL Open_New_File( Report_Name, Report_IDs(iRF_Setup), Suggested_Number)
    END IF

    CALL Output_Params_Report( Report_IDs(iRF_Setup) )

    IF ( .NOT. lPF_Core_Flags(iPF_Core_Newtonian_Mode) ) THEN
        CALL Output_NL_Solver_Report( Report_IDs(iRF_Setup) )
    END IF

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
1301 FORMAT('         Running in Poisson Solver Mode         ')

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


1415 FORMAT(/'----------------- Domain Limits ---------------')
1416 FORMAT(' Inner Radius = ',ES12.4E3,' cm')
1417 FORMAT(' Outer Radius = ',ES12.4E3,' cm')


IF ( lPF_IO_Flags(iPF_IO_Print_Setup) ) THEN
    WRITE(*,1400)
    WRITE(*,1401)
    IF ( lPF_Core_Flags(iPF_Core_Newtonian_Mode) ) WRITE(*,1301)
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
    WRITE(*,1415)
    WRITE(*,1416)R_Inner/Centimeter
    WRITE(*,1417)R_Outer/Centimeter
    WRITE(*,1400)
END IF

IF ( lPF_IO_Flags(iPF_IO_Write_Setup) ) THEN

    WRITE(Report_ID,1400)
    WRITE(Report_ID,1401)
    IF ( lPF_Core_Flags(iPF_Core_Newtonian_Mode) ) WRITE(Report_ID,1301)
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
!        Output_NL_Solver_Report                                        !
!                                                                       !
!#######################################################################!
SUBROUTINE Output_NL_Solver_Report( Report_ID )

INTEGER, INTENT(IN)                                 ::  Report_ID

1500 FORMAT(/)
1501 FORMAT('--------- Non-Linear Solver Parameters ---------')
1502 FORMAT(' Solver Method : ',A /)
1503 FORMAT(' Maximum Iterations   = ',I4.0)
1504 FORMAT(' Convergence Criteria = ',ES12.4E3)

1505 FORMAT(' Anderson Acceleration Memory Parameter = ',I3.1)



IF ( lPF_IO_Flags(iPF_IO_Print_Setup)  ) THEN
    WRITE(*,1501)
    WRITE(*,1502)Method_Names(Method_Flag)
    WRITE(*,1503)Max_Iterations
    WRITE(*,1504)Convergence_Criteria
    IF ( Method_Flag >= 2 ) THEN
        WRITE(*,1505) FP_Anderson_M
    END IF
    WRITE(*,1500)
END IF



IF ( lPF_IO_Flags(iPF_IO_Write_Setup) ) THEN
    WRITE(Report_ID,1501)
    WRITE(Report_ID,1502)Method_Names(Method_Flag)
    WRITE(Report_ID,1503)Max_Iterations
    WRITE(Report_ID,1504)Convergence_Criteria
    IF ( Method_Flag >= 2 ) THEN
        WRITE(Report_ID,1505) FP_Anderson_M
    END IF
END IF


END SUBROUTINE Output_NL_Solver_Report







!+201+##################################################################!
!                                                                       !
!       PRINT_AMReX_Setup                                              !
!                                                                       !
!#######################################################################!
SUBROUTINE PRINT_AMReX_Setup( )

1400 FORMAT(/)
1401 FORMAT('--------------- AMReX PARAMETERS ---------------')
1402 FORMAT(' Maximum Level of Refinement = ',I5.1)
!1403 FORMAT(' # of Coarse Radial Elements = ',I5.1)
1404 FORMAT(' # of Leaf Radial Elements   = ',I5.1)

WRITE(*,1400)
WRITE(*,1401)
WRITE(*,1402)AMReX_Max_Level
!WRITE(*,1403)
WRITE(*,1404)iNumLeafElements
WRITE(*,1400)


END SUBROUTINE PRINT_AMReX_Setup




END MODULE IO_Setup_Report_Module
