   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Initialization_Derived                                                       !##!
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
            ONLY : idp

USE Poseidon_Parameters, &
            ONLY :  Domain_Dim,             &
                    Degree,                 &
                    L_Limit,                &
                    Verbose_Flag,           &
                    Num_Vars

USE Poseidon_Message_Routines_Module, &
            ONLY :  Init_Message

USE Variables_Mesh, &
            ONLY :  Num_R_Elements

USE Variables_Derived, &
            ONLY :  Num_R_Nodes,            &
                    Num_R_Nodesp1,          &
                    Var_Dim,                &
                    Elem_Var_Dim,           &
                    Block_Var_Dim,          &
                    LM_Length,              &
                    LM_Short_Length,        &
                    ULM_Length,             &
                    Prob_Dim,               &
                    Elem_Prob_Dim,          &
                    Elem_Prob_Dim_Sqr,      &
                    Block_Prob_Dim,         &
                    Num_Off_Diagonals,      &
                    iVB_Prob_Dim,           &
                    iVB_Elem_Prob_Dim

USE Variables_Matrices, &
            ONLY :  iMB_Bandwidth,         &
                    iMB_Diagonals,         &
                    Laplace_NNZ

USE Variables_AMReX_Core, &
            ONLY :  iNumLeafElements

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_MTGV_Flags,    &
                    iPF_Init_MTGV_Derived


IMPLICIT NONE

CONTAINS
 !+101+####################################################!
!                                                           !
!       Initialize_Derived                                  !
!                                                           !
 !#########################################################!
SUBROUTINE Initialize_Derived()

IF ( Verbose_Flag ) CALL Init_Message('Calculating Derived Variables, Native.')


IF ( DOMAIN_DIM == 1 ) THEN

    LM_Length = 1
    LM_Short_Length = 1
    
ELSE IF ( DOMAIN_DIM == 2 ) THEN

    LM_Length = L_Limit + 1
    LM_Short_Length = L_Limit + 1

ELSE IF ( DOMAIN_DIM == 3 ) THEN

    LM_Length = (L_Limit + 1)*(L_Limit + 1)
    LM_Short_Length = (L_Limit+1)*(L_Limit+2)/2
    
END IF



NUM_R_NODES             = DEGREE*NUM_R_ELEMENTS + 1
Num_R_Nodesp1           = Num_R_Nodes + 1



!
!   VAR_DIM - Length of vector required to hold coefficients for 1 variable.
!
VAR_DIM             = LM_LENGTH*NUM_R_NODES
ELEM_VAR_DIM        = LM_LENGTH*(DEGREE + 1)


!
!   PROB_DIM - Length of vector required to hold coefficients for all variables
!
ULM_LENGTH          = Num_Vars*LM_LENGTH
PROB_DIM            = Num_Vars*VAR_DIM
ELEM_PROB_DIM       = Num_Vars*ELEM_VAR_DIM
ELEM_PROB_DIM_SQR   = ELEM_PROB_DIM*ELEM_PROB_DIM


iVB_Prob_Dim        = 3*Var_Dim
iVB_Elem_Prob_Dim  = 3*Elem_Var_Dim


NUM_OFF_DIAGONALS = ULM_LENGTH*(DEGREE + 1) - 1

lPF_Init_MTGV_Flags(iPF_Init_MTGV_Derived) = .TRUE.

END SUBROUTINE Initialize_Derived










 !+201+####################################################!
!                                                           !
!       Initialize_Derived_AMReX                            !
!                                                           !
 !#########################################################!
SUBROUTINE Initialize_Derived_AMReX_Part1

IF ( Verbose_Flag ) CALL Init_Message('Calculating Derived Variables, AMReX - Part 1.')

IF ( DOMAIN_DIM == 1 ) THEN

    LM_Length = 1
    LM_Short_Length = 1
    
ELSE IF ( DOMAIN_DIM == 2 ) THEN

    LM_Length = L_Limit + 1
    LM_Short_Length = L_Limit + 1

ELSE IF ( DOMAIN_DIM == 3 ) THEN

    LM_Length = (L_Limit + 1)*(L_Limit + 1)
    LM_Short_Length = (L_Limit+1)*(L_Limit+2)/2
    
END IF


!LM_LENGTH = (L_LIMIT + 1)*(L_LIMIT + 1)
!LM_Short_Length = (L_Limit+1)*(L_Limit+2)/2

!
!   VAR_DIM - Length of vector required to hold coefficients for 1 variable.
!
ELEM_VAR_DIM        = LM_LENGTH*(DEGREE + 1)


!
!   PROB_DIM - Length of vector required to hold coefficients for all variables
!
ULM_LENGTH          = Num_Vars*LM_LENGTH
ELEM_PROB_DIM       = Num_Vars*ELEM_VAR_DIM
ELEM_PROB_DIM_SQR   = ELEM_PROB_DIM*ELEM_PROB_DIM


iVB_Elem_Prob_Dim  = 3*Elem_Var_Dim


NUM_OFF_DIAGONALS   = ULM_LENGTH*(DEGREE + 1) - 1


lPF_Init_MTGV_Flags(iPF_Init_MTGV_Derived) = .TRUE.


END SUBROUTINE Initialize_Derived_AMReX_Part1




 !+201+####################################################!
!                                                           !
!       Initialize_Derived_AMReX                            !
!                                                           !
 !#########################################################!
SUBROUTINE Initialize_Derived_AMReX_Part2

IF ( Verbose_Flag ) CALL Init_Message('Calculating Derived Variables, AMReX - Part 2.')


Num_R_Elements = iNumLeafElements



! Determine Derived Varaibles
Num_R_Nodes         = Degree*NUM_R_ELEMENTS + 1
Num_R_Nodesp1       = Num_R_Nodes + 1


VAR_DIM             = LM_Length*NUM_R_NODES
PROB_DIM            = NUM_VARS*VAR_DIM
iVB_Prob_Dim        = 3*Var_Dim

iMB_Diagonals       = iVB_Elem_Prob_Dim-1
iMB_Bandwidth       = 2*iMB_Diagonals+1

Laplace_NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1


END SUBROUTINE Initialize_Derived_AMReX_Part2







END MODULE Initialization_Derived


