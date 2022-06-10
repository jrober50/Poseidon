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
                    Num_CFA_Vars

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
                    ULM_Length,             &
                    Prob_Dim,               &
                    Elem_Prob_Dim,          &
                    Elem_Prob_Dim_Sqr,      &
                    Block_Prob_Dim,         &
                    Num_Off_Diagonals,      &
                    Beta_Prob_Dim,          &
                    Beta_Elem_Prob_Dim

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

IF ( Verbose_Flag ) CALL Init_Message('Calculating Derived Variables.')


IF ( DOMAIN_DIM == 1 ) THEN

    LM_LENGTH = 1

ELSE IF ( DOMAIN_DIM == 2 ) THEN

    LM_LENGTH = L_LIMIT + 1

ELSE IF ( DOMAIN_DIM == 3 ) THEN

    LM_LENGTH = (L_LIMIT + 1)*(L_LIMIT + 1)

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
ULM_LENGTH          = NUM_CFA_VARS*LM_LENGTH
PROB_DIM            = NUM_CFA_VARS*VAR_DIM
ELEM_PROB_DIM       = NUM_CFA_VARS*ELEM_VAR_DIM
ELEM_PROB_DIM_SQR   = ELEM_PROB_DIM*ELEM_PROB_DIM


Beta_Prob_Dim       = 3*Var_Dim
Beta_Elem_Prob_Dim  = 3*Elem_Var_Dim


NUM_OFF_DIAGONALS = ULM_LENGTH*(DEGREE + 1) - 1

lPF_Init_MTGV_Flags(iPF_Init_MTGV_Derived) = .TRUE.

END SUBROUTINE Initialize_Derived










 !+201+####################################################!
!                                                           !
!       Initialize_Derived_AMReX                            !
!                                                           !
 !#########################################################!
SUBROUTINE Initialize_Derived_AMReX

IF ( Verbose_Flag ) CALL Init_Message('Calculating Derived Variables.')


IF ( DOMAIN_DIM == 1 ) THEN

    LM_LENGTH = 1

ELSE IF ( DOMAIN_DIM == 2 ) THEN

    LM_LENGTH = L_LIMIT + 1

ELSE IF ( DOMAIN_DIM == 3 ) THEN

    LM_LENGTH = (L_LIMIT + 1)*(L_LIMIT + 1)

END IF





!
!   VAR_DIM - Length of vector required to hold coefficients for 1 variable.
!
ELEM_VAR_DIM        = LM_LENGTH*(DEGREE + 1)


!
!   PROB_DIM - Length of vector required to hold coefficients for all variables
!
ULM_LENGTH          = NUM_CFA_VARS*LM_LENGTH
ELEM_PROB_DIM       = NUM_CFA_VARS*ELEM_VAR_DIM
ELEM_PROB_DIM_SQR   = ELEM_PROB_DIM*ELEM_PROB_DIM


Beta_Elem_Prob_Dim  = 3*Elem_Var_Dim


NUM_OFF_DIAGONALS   = ULM_LENGTH*(DEGREE + 1) - 1


lPF_Init_MTGV_Flags(iPF_Init_MTGV_Derived) = .TRUE.


END SUBROUTINE Initialize_Derived_AMReX






END MODULE Initialization_Derived


