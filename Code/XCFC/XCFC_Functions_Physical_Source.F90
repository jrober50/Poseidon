   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE XCFC_Functions_Physical_Source_Module                                 !##!
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


USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,              &
                    iU_LF,              &
                    iU_S1,              &
                    iU_S2,              &
                    iU_S3,              &
                    iU_X1,              &
                    iU_X2,              &
                    iU_X3,              &
                    iVB_S,              &
                    iVB_X,              &
                    iS_E,               &
                    iS_S,               &
                    iS_S1,              &
                    iS_S2,              &
                    iS_S3


USE Variables_Source, &
            ONLY :  Block_Source_E,     &
                    Block_Source_S,     &
                    Block_Source_Si

USE Variables_Quadrature, &
            ONLY :  NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    NUM_TP_QUAD_POINTS,         &
                    Num_Quad_DOF

USE Functions_Domain_Maps, &
            ONLY :  Map_To_tpd

USE Quadrature_Mapping_Functions, &
            ONLY :  Quad_Map

USE Variables_AMReX_Source, &
            ONLY :  Source_PTR

USE Variables_MPI, &
                ONLY :  myID_Poseidon,      &
                        nPROCS_Poseidon,    &
                        Poseidon_Comm_World

USE Poseidon_MPI_Utilities_Module, &
                ONLY :  STOP_MPI,               &
                        MPI_Master_Print,       &
                        MPI_All_Print

#ifdef POSEIDON_AMREX_FLAG
use amrex_fort_module, &
                ONLY :  amrex_spacedim
#endif

IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!                                                     				!
!                                                                               !
!###############################################################################!
SUBROUTINE Get_Physical_Source( Source, iU, iE)

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points), INTENT(OUT)    :: Source
INTEGER,                                                     INTENT(IN)     :: iU
INTEGER,   DIMENSION(3),                                     INTENT(IN)     :: iE

#ifdef POSEIDON_AMREX_FLAG

    CALL Get_Physical_Source_AMReX( Source, iU, iE)

#else

    CALL Get_Physical_Source_Native( Source, iU, iE)

#endif

END SUBROUTINE Get_Physical_Source




!+101+##########################################################################!
!                                                                               !
!                                                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Get_Physical_Source_Native( Source, iU, iE)

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points), INTENT(OUT)    :: Source
INTEGER,                                                     INTENT(IN)     :: iU
INTEGER,   DIMENSION(3),                                     INTENT(IN)     :: iE


INTEGER                                                         :: rd, td, pd, tpd
INTEGER                                                         :: Here

IF ( iU == 1 ) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd  = Map_To_tpd(td,pd)
        Here = Quad_Map(rd,td,pd)
        Source(tpd,rd) = Block_Source_E(Here,iE(1),iE(2),iE(3))

    END DO ! pd
    END DO ! td
    END DO ! rd


ELSEIF ( iU == 2 ) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd  = Map_To_tpd(td,pd)
        Here = Quad_Map(rd,td,pd)
        Source(tpd,rd) = Block_Source_E(Here,iE(1),iE(2),iE(3))                &
                       + 2.0_idp*Block_Source_S(Here,iE(1),iE(2),iE(3))

!        PRINT*,iE, rd,td,pd
!        PRINT*,Block_Source_E(Here,iE(1),iE(2),iE(3)),                &
!                2.0_idp*Block_Source_S(Here,iE(1),iE(2),iE(3))

    END DO ! pd
    END DO ! td
    END DO ! rd

ELSEIF ( iU == 3) THEN

   DO rd = 1,NUM_R_QUAD_POINTS
   DO td = 1,NUM_T_QUAD_POINTS
   DO pd = 1,NUM_P_QUAD_POINTS

       tpd  = Map_To_tpd(td,pd)
       Here = Quad_Map(rd,td,pd)

!       PRINT*,iE,rd,td,pd,Block_Source_Si(Here,iE(1),iE(2),iE(3),1)
       Source(tpd,rd) = Block_Source_Si(Here,iE(1),iE(2),iE(3),1)

   END DO ! pd
   END DO ! td
!
!    PRINT*,iE,rd
!    PRINT*,Source(:,rd)

   END DO ! rd

ELSE IF ( iU == 4) THEN

   DO rd = 1,NUM_R_QUAD_POINTS
   DO td = 1,NUM_T_QUAD_POINTS
   DO pd = 1,NUM_P_QUAD_POINTS

       tpd  = Map_To_tpd(td,pd)
       Here = Quad_Map(rd,td,pd)
       Source(tpd,rd) = Block_Source_Si(Here,iE(1),iE(2),iE(3),2)


   END DO ! pd
   END DO ! td
   END DO ! rd

ELSE IF ( iU == 5 ) THEN

   DO rd = 1,NUM_R_QUAD_POINTS
   DO td = 1,NUM_T_QUAD_POINTS
   DO pd = 1,NUM_P_QUAD_POINTS

       tpd  = Map_To_tpd(td,pd)
       Here = Quad_Map(rd,td,pd)
       Source(tpd,rd) = Block_Source_Si(Here,iE(1),iE(2),iE(3),3)


   END DO ! pd
   END DO ! td
   END DO ! rd

END IF

END SUBROUTINE Get_Physical_Source_Native




!+101+##########################################################################!
!                                                                               !
!                                                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Get_Physical_Source_AMReX( Source, iU, iE)

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points), INTENT(OUT)    :: Source
INTEGER,                                                     INTENT(IN)     :: iU
INTEGER,   DIMENSION(3),                                     INTENT(IN)     :: iE
INTEGER                                                         :: rd, td, pd, tpd
INTEGER                                                         :: Here, There

INTEGER                                                         :: ierr
INTEGER                                                         ::  iEOff(3)


#ifdef POSEIDON_AMREX_FLAG
IF ( amrex_spacedim == 1 ) THEN
    iEoff(2:3) = 1
ELSEIF ( amrex_spacedim == 2) THEN
    iEoff(2)   = iE(2)
    iEoff(3)   = 1
ELSEIF ( amrex_spacedim == 3 ) THEN
    iEoff(2:3) = iE(2:3)
END IF
#else
    iEOff(2:3) = iE(2:3)
    PRINT*,"Warining Get_Physical Source_AMReX is being called with the POSEIDON_AMREX_FLAG = False."
#endif

!PRINT*,"In Physical_Source_AMReX, MyID ",myID_Poseidon," - ",iE(1),iE(2),iE(3)
!CALL MPI_Barrier(Poseidon_Comm_World, ierr )

IF ( iU == 1 ) THEN


    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = Map_To_tpd(td,pd)
        Here = (iS_E-1)*NUM_Quad_DOF + Quad_Map(rd,td,pd)
        Source(tpd,rd) = Source_PTR(iE(1),iEOff(2),iEOff(3),Here)
!        PRINT*,Source_PTR(iE(1),iE(2),iE(3),Here)

    END DO ! pd
    END DO ! td
    END DO ! rd

ELSEIF ( iU == 2 ) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = Map_To_tpd(td,pd)
        Here = (iS_E-1)*NUM_Quad_DOF + Quad_Map(rd,td,pd)
        There = (iS_S-1)*NUM_Quad_DOF + Quad_Map(rd,td,pd)
        Source(tpd,rd) = Source_PTR(iE(1),iEOff(2),iEOff(3),Here)  &
                    + 2.0_idp*Source_PTR(iE(1),iEOff(2),iEOff(3),There)

!        PRINT*,iE, rd,td,pd
!        PRINT*,Source_PTR(iE(1),iE(2),iE(3),Here),              &
!                2.0_idp*Source_PTR(iE(1),iE(2),iE(3),There)
    END DO ! pd
    END DO ! td
    END DO ! rd

ELSE IF ( iU == 3) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = Map_To_tpd(td,pd)
        Here = (iS_S1-1)*NUM_Quad_DOF + Quad_Map(rd,td,pd)

        Source(tpd,rd) = Source_PTR(iE(1),iEOff(2),iEOff(3),Here)
    END DO ! pd
    END DO ! td
    END DO ! rd

ELSE IF ( iU == 4) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = Map_To_tpd(td,pd)
        Here = (iS_S2-1)*NUM_Quad_DOF + Quad_Map(rd,td,pd)

        Source(tpd,rd) = Source_PTR(iE(1),iEOff(2),iEOff(3),Here)


    END DO ! pd
    END DO ! td
    END DO ! rd

ELSE IF ( iU == 5 ) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = Map_To_tpd(td,pd)
        Here = (iS_S3-1)*NUM_Quad_DOF + Quad_Map(rd,td,pd)

        Source(tpd,rd) = Source_PTR(iE(1),iEOff(2),iEOff(3),Here)
        

    END DO ! pd
    END DO ! td
    END DO ! rd

END IF

END SUBROUTINE Get_Physical_Source_AMReX












END MODULE XCFC_Functions_Physical_Source_Module
