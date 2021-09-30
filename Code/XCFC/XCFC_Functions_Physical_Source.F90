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

USE FP_Functions_Mapping, &
            ONLY :  FP_FEM_Node_Map,    &
                    FP_tpd_Map

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

IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!                                                     				!
!                                                                               !
!###############################################################################!
SUBROUTINE Get_Physical_Source( Source, iU, RE, TE, PE)

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points), INTENT(OUT)    :: Source
INTEGER,                                                     INTENT(IN)     :: iU
INTEGER,                                                     INTENT(IN)     :: RE, TE, PE

#ifdef POSEIDON_AMREX_FLAG

    CALL Get_Physical_Source_AMReX( Source, iU, RE, TE, PE)

#else

    CALL Get_Physical_Source_Native( Source, iU, RE, TE, PE)

#endif

END SUBROUTINE Get_Physical_Source




!+101+##########################################################################!
!                                                                               !
!                                                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Get_Physical_Source_Native( Source, iU, RE, TE, PE)

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points), INTENT(OUT)    :: Source
INTEGER,                                                     INTENT(IN)     :: iU
INTEGER,                                                     INTENT(IN)     :: RE, TE, PE


INTEGER                                                         :: rd, td, pd, tpd


IF ( iU == 1 ) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = FP_tpd_Map(td,pd)
        Source(tpd,rd) = Block_Source_E(rd,td,pd,re,te,pe)

    END DO ! pd
    END DO ! td
    END DO ! rd


ELSEIF ( iU == 2 ) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = FP_tpd_Map(td,pd)
        Source(tpd,rd) = Block_Source_E(rd,td,pd,re,te,pe)                &
                       + 2.0_idp*Block_Source_S(rd,td,pd,re,te,pe)



    END DO ! pd
    END DO ! td
    END DO ! rd

ELSEIF ( iU == 3) THEN

   DO rd = 1,NUM_R_QUAD_POINTS
   DO td = 1,NUM_T_QUAD_POINTS
   DO pd = 1,NUM_P_QUAD_POINTS

       tpd = FP_tpd_Map(td,pd)
       Source(tpd,rd) = Block_Source_Si(rd,td,pd,re,te,pe,1)

   END DO ! pd
   END DO ! td
   END DO ! rd

ELSE IF ( iU == 4) THEN

   DO rd = 1,NUM_R_QUAD_POINTS
   DO td = 1,NUM_T_QUAD_POINTS
   DO pd = 1,NUM_P_QUAD_POINTS

       tpd = FP_tpd_Map(td,pd)
       Source(tpd,rd) = Block_Source_Si(rd,td,pd,re,te,pe,2)


   END DO ! pd
   END DO ! td
   END DO ! rd

ELSE IF ( iU == 5 ) THEN

   DO rd = 1,NUM_R_QUAD_POINTS
   DO td = 1,NUM_T_QUAD_POINTS
   DO pd = 1,NUM_P_QUAD_POINTS

       tpd = FP_tpd_Map(td,pd)
       Source(tpd,rd) = Block_Source_Si(rd,td,pd,re,te,pe,3)


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
SUBROUTINE Get_Physical_Source_AMReX( Source, iU, RE, TE, PE)

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points), INTENT(OUT)    :: Source
INTEGER,                                                     INTENT(IN)     :: iU
INTEGER,                                                     INTENT(IN)     :: RE, TE, PE
INTEGER                                                         :: rd, td, pd, tpd
INTEGER                                                         :: Here, There

INTEGER                                                         :: ierr

!PRINT*,"In Physical_Source_AMReX, MyID ",myID_Poseidon," - ",re,te,pe
!CALL MPI_Barrier(Poseidon_Comm_World, ierr )

IF ( iU == 1 ) THEN


    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = FP_tpd_Map(td,pd)
        Here = (iS_E-1)*NUM_Quad_DOF                       &
             + (pd-1)*Num_R_Quad_Points*Num_T_Quad_Points   &
             + (td-1)*Num_R_Quad_Points                     &
             + rd

        Source(tpd,rd) = Source_PTR(re+1,te+1,pe+1,Here)


    END DO ! pd
    END DO ! td
    END DO ! rd

ELSEIF ( iU == 2 ) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = FP_tpd_Map(td,pd)
        Here = (iS_E-1)*NUM_Quad_DOF                       &
             + (pd-1)*Num_R_Quad_Points*Num_T_Quad_Points   &
             + (td-1)*Num_R_Quad_Points                     &
             + rd
        There = (iS_S-1)*NUM_Quad_DOF                       &
                + (pd-1)*Num_R_Quad_Points*Num_T_Quad_Points   &
                + (td-1)*Num_R_Quad_Points                     &
                + rd
        Source(tpd,rd) = Source_PTR(re+1,te+1,pe+1,Here)  &
                    + 2.0_idp*Source_PTR(re+1,te+1,pe+1,There)


    END DO ! pd
    END DO ! td
    END DO ! rd

ELSE IF ( iU == 3) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = FP_tpd_Map(td,pd)
        Here = (iS_S1-1)*NUM_Quad_DOF                       &
             + (pd-1)*Num_R_Quad_Points*Num_T_Quad_Points   &
             + (td-1)*Num_R_Quad_Points                     &
             + rd

!        PRINT*,"myID ",myID_Poseidon," - ",Here,re+1,te+1,pe+1
        Source(tpd,rd) = Source_PTR(re+1,te+1,pe+1,Here)

!        PRINT*,"myID ",myID_Poseidon," - ",Source(tpd,rd)

    END DO ! pd
    END DO ! td
    END DO ! rd

ELSE IF ( iU == 4) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = FP_tpd_Map(td,pd)
        Here = (iS_S2-1)*NUM_Quad_DOF                       &
             + (pd-1)*Num_R_Quad_Points*Num_T_Quad_Points   &
             + (td-1)*Num_R_Quad_Points                     &
             + rd

        Source(tpd,rd) = Source_PTR(re+1,te+1,pe+1,Here)


    END DO ! pd
    END DO ! td
    END DO ! rd

ELSE IF ( iU == 5 ) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = FP_tpd_Map(td,pd)
        Here = (iS_S3-1)*NUM_Quad_DOF                       &
             + (pd-1)*Num_R_Quad_Points*Num_T_Quad_Points   &
             + (td-1)*Num_R_Quad_Points                     &
             + rd

        Source(tpd,rd) = Source_PTR(re+1,te+1,pe+1,Here)
        

    END DO ! pd
    END DO ! td
    END DO ! rd

END IF

END SUBROUTINE Get_Physical_Source_AMReX












END MODULE XCFC_Functions_Physical_Source_Module
