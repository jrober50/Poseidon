   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Remesh_Module                                          	     !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!    +101+   Fill_Coeff_Vector_From_Copy                                  !##!
!##!                                                                         !##!
!##!    +201+   Make_Remesh_Copies                                           !##!
!##!    +202+   Destroy_Remesh_Copies                                        !##!
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
            ONLY :  iU_CF,      &
                    iU_LF,      &
                    iU_S1,      &
                    iU_S2,      &
                    iU_S3,      &
                    iU_X1,      &
                    iU_X2,      &
                    iU_X3,      &
                    iVB_S,      &
                    iVB_X

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag,           &
                    Degree

USE Poseidon_Message_Routines_Module, &
            ONLY :  Run_Message

USE Variables_Vectors, &
            ONLY :  dVA_Coeff_Vector,        &
                    dVB_Coeff_Vector

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    Num_T_Elements,         &
                    Num_P_Elements,         &
                    rlocs

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Int_R_Locations

USE Variables_Derived, &
            ONLY :  Num_R_Nodes,            &
                    LM_Length,              &
                    iVB_Prob_Dim

USE Variables_FEM_Module, &
            ONLY :  FEM_Node_xlocs

USE Maps_X_Space, &
            ONLY :  Map_To_X_Space

USE Maps_Domain, &
            ONLY :  Map_To_FEM_Node

USE Maps_Fixed_Point, &
            ONLY :  FP_Array_Map_TypeB

USE Functions_Math, &
            ONLY :  Lagrange_Poly

USE Timer_Routines_Module, &
            ONLY :  TimerStart,     &
                    TimerSTop

USE Timer_Variables_Module, &
            ONLY :  Timer_Remesh,                   &
                    Timer_Remesh_MakeCopies,        &
                    Timer_Remesh_FillTotal,         &
                    Timer_Remesh_MakeLambdaArray,   &
                    Timer_Remesh_FillTypeA,         &
                    Timer_Remesh_FillX,             &
                    Timer_Remesh_FillS,             &
                    Timer_Remesh_DestroyCopies

USE Flags_Core_Module, &
            ONLY :  iPF_Core_Flags,             &
                    iPF_Core_Method_Mode,       &
                    iPF_Core_Method_Newtonian

IMPLICIT NONE


INTEGER                                         ::  iVB_Prob_Dim_Old
INTEGER                                         ::  Num_R_Nodes_Old
INTEGER,        DIMENSION(3)                    ::  NE_Old


REAL(idp),      DIMENSION(:),       ALLOCATABLE ::  rlocs_Old

REAL(idp),   DIMENSION(:,:,:),   ALLOCATABLE ::  dVA_Coeff_Old
REAL(idp),   DIMENSION(:,:),     ALLOCATABLE ::  dVB_Coeff_Old


REAL(idp),      DIMENSION(:,:),     ALLOCATABLE ::  lm_at_rn

CONTAINS



 !+101+####################################################!
!                                                           !
!          Fill_Coeff_Vector_From_Copy                      !
!                                                           !
 !#########################################################!
SUBROUTINE Fill_Coeff_Vector_From_Copy()

INTEGER                                     ::  re
INTEGER                                     ::  d, d_old, lm
INTEGER                                     ::  iU, iVB
INTEGER                                     ::  New_Node
INTEGER                                     ::  FEM_Node
INTEGER                                     ::  re_Old

INTEGER                                     ::  Here_New
INTEGER                                     ::  Here_Old, There_Old

REAL(idp)                                   ::  x
REAL(idp)                                   ::  DROT
REAL(idp)                                   ::  TMP

REAL(idp),  DIMENSION(0:Degree)             ::  Cur_R_Locs


INTEGER,    DIMENSION(:),   ALLOCATABLE     ::  rn_in_re
REAL(idp),  DIMENSION(:,:), ALLOCATABLE     ::  lm_at_rn


IF ( Verbose_Flag ) CALL Run_Message('Beginning Remeshing of Coefficient Vector')
CALL TimerStart(Timer_Remesh_FillTotal)
CALL TimerStart(Timer_Remesh_MakeLambdaArray)

ALLOCATE( lm_at_rn(0:Degree,1:Num_R_Nodes) )
ALLOCATE( rn_in_re(1:Num_R_Nodes) )

DO re = 0,Num_R_Elements-1  ! Cycle through new elements

    DROT = 0.5_idp * (rlocs(re+1) - rlocs(re))
    Cur_R_Locs(:) = DROT * (FEM_Node_xlocs(:)+1.0_idp) + rlocs(re)

    DO d = 0,Degree                 ! Cycle through each node in new element.
    DO re_old = 0,NE_Old(1)-1       ! Cycle through old elements

        IF ( Cur_R_Locs(d) .GE. rlocs_Old(re_old) ) THEN       ! Check if current location is
        IF ( Cur_R_Locs(d) .LE. rlocs_Old(re_old+1) ) THEN     ! within old element, re_old.

            ! Map value to the old elements reference element
            x = Map_To_X_Space(rlocs_Old(re_old),rlocs_Old(re_old+1),Cur_R_Locs(d))

            FEM_Node = Map_To_FEM_Node(re,d)   ! Determine new node number

            ! Evaluate and Store the Lagrange Polynomials defined by the
            ! FEM nodes on refence element, at the location of the new node, FEM_Node,
            ! within the element; Lambda_m at location r_n.
            
            lm_at_rn(:,FEM_Node) = Lagrange_Poly(x,Degree,FEM_Node_xlocs)


            ! Also needed to reconstruct the solution is the number
            ! of the old element containing the new node.
            rn_in_re(FEM_Node) = re_old

        END IF  ! If below right edge
        END IF  ! If above left edge

    END DO ! re_old
    END DO ! d
END DO ! re


CALL TimerStop(Timer_Remesh_MakeLambdaArray)
CALL TimerStart(Timer_Remesh_FillTypeA)


! Work through Type A Coeffs
DO iU = 1,2
DO lm = 1,LM_Length
DO New_Node = 1,Num_R_Nodes

    Here_Old  = Map_To_FEM_Node(rn_in_re(New_Node),0)
    There_Old = Map_To_FEM_Node(rn_in_re(New_Node),Degree)

    dVA_Coeff_Vector(New_Node,lm,iU) = SUM( dVA_Coeff_Old(Here_Old:There_Old,lm,iU)     &
                                            * lm_at_rn(:,New_Node)                      )


END DO  ! New_Node
END DO  ! lm
END DO  ! iU


CALL TimerStop(Timer_Remesh_FillTypeA)


IF ( iPF_Core_Flags(iPF_Core_Method_Mode) .NE. iPF_Core_Method_Newtonian ) THEN
    CALL TimerStart(Timer_Remesh_FillX)


    ! Work through Type B Coeffs

    dVB_Coeff_Vector = 0.0_idp

    iVB = iVB_X
    DO re = 0,Num_R_Elements-1
    DO d = 0,Degree
    DO iU = iU_X1,iU_X3
    DO lm = 1,LM_Length

        Here_New = FP_Array_Map_TypeB(iU,iVB,re,d,lm)
        FEM_Node = Map_To_FEM_Node(re,d)
        re_old   = rn_in_re(FEM_Node)


        TMP = 0.0_idp
        DO d_old = 0,DEGREE
            Here_Old = FP_Array_Map_TypeB(iU,iVB,re_old,d_old,lm)

            TMP = TMP                            &
                + dVB_Coeff_Old(Here_Old,iVB)    &
                * lm_at_rn(d_old,FEM_Node)


        END DO ! d_Old

        dVB_Coeff_Vector(Here_New,iVB) = TMP

    END DO  ! lm
    END DO  ! iU
    END DO  ! d
    END DO  ! re



    CALL TimerStop(Timer_Remesh_FillX)
    CALL TimerStart(Timer_Remesh_FillS)

    iVB = iVB_S
    DO re = 0,Num_R_Elements-1
    DO d = 0,Degree
    DO iU = iU_S1,iU_S3
    DO lm = 1,LM_Length

        Here_New = FP_Array_Map_TypeB(iU,iVB,re,d,lm)
        FEM_Node = Map_To_FEM_Node(re,d)
        re_old   = rn_in_re(FEM_Node)


        TMP = 0.0_idp
        DO d_old = 0,DEGREE
            Here_Old = FP_Array_Map_TypeB(iU,iVB,re_old,d_old,lm)

            TMP = TMP                            &
                + dVB_Coeff_Old(Here_Old,iVB)    &
                * lm_at_rn(d_old,FEM_Node)


        END DO ! d_Old

        dVB_Coeff_Vector(Here_New,iVB) = TMP

    END DO  ! lm
    END DO  ! iU
    END DO  ! d
    END DO  ! re
END IF

CALL TimerStop(Timer_Remesh_FillS)

DEALLOCATE( lm_at_rn )
DEALLOCATE( rn_in_re )


CALL TimerStart(Timer_Remesh_FillTotal)

END SUBROUTINE Fill_Coeff_Vector_From_Copy







 !+201+####################################################!
!                                                           !
!          Make_Remesh_Copies                               !
!                                                           !
 !#########################################################!
SUBROUTINE Make_Remesh_Copies()


IF ( Verbose_Flag ) CALL Run_Message('Making Copies of Variables for Remesh.')
CALL TimerStart(Timer_Remesh_MakeCopies)


iVB_Prob_Dim_Old  = iVB_Prob_Dim
Num_R_Nodes_Old   = Num_R_Nodes
NE_Old(1)         = Num_R_Elements
NE_Old(2)         = Num_T_Elements
NE_Old(3)         = Num_P_Elements

ALLOCATE( rlocs_old(0:Num_R_Elements) )
ALLOCATE( dVA_Coeff_Old(1:Num_R_Nodes_Old,1:LM_Length,1:2) )

rlocs_Old = rlocs
dVA_Coeff_Old = dVA_Coeff_Vector



IF ( iPF_Core_Flags(iPF_Core_Method_Mode) .NE. iPF_Core_Method_Newtonian ) THEN
    ALLOCATE( dVB_Coeff_Old(1:iVB_Prob_Dim_Old,1:2) )
    dVB_Coeff_Old = dVB_Coeff_Vector
END IF


CALL TimerStop(Timer_Remesh_MakeCopies)

END SUBROUTINE Make_Remesh_Copies



 !+202+####################################################!
!                                                           !
!          Destroy_Remesh_Copies                            !
!                                                           !
 !#########################################################!
SUBROUTINE Destroy_Remesh_Copies()

IF ( Verbose_Flag ) CALL Run_Message('Destroying Remesh Variables.')
CALL TimerStart(Timer_Remesh_DestroyCopies)

DEALLOCATE( rlocs_old )
DEALLOCATE( dVA_Coeff_Old )

IF ( iPF_Core_Flags(iPF_Core_Method_Mode) .NE. iPF_Core_Method_Newtonian ) THEN
    DEALLOCATE( dVB_Coeff_Old )
END IF

CALL TimerStop(Timer_Remesh_DestroyCopies)

END SUBROUTINE Destroy_Remesh_Copies





END MODULE Poseidon_Remesh_Module
