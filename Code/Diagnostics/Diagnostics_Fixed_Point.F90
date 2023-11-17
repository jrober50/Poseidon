   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Diagnostics_Fixed_Point                                                      !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
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
            ONLY :  L_Limit
                    
Use Variables_Derived, &
            ONLY :  Num_R_Nodes,            &
                    LM_Length
                    
USE Variables_Matrices,  &
            ONLY :  Factored_NNZ,               &
                    Laplace_Factored_VAL,       &
                    Laplace_Factored_ROW,       &
                    Laplace_Factored_COL
                    
USE Variables_FP,  &
            ONLY :  FP_Residual_Vector,     &
                    FP_Laplace_Vector,      &
                    FP_Iter_Load_Storage,   &
                    FP_Iter_Matrix_Storage, &
                    Resid_Norms,            &
                    Update_Norms
                    
USE Variables_Vectors,  &
            ONLY :  dVA_Coeff_Vector
            
USE Maps_Domain, &
            ONLY :  Map_To_lm
            
USE Functions_Matrix, &
            ONLY :  Matrix_CCS_MVMult,          &
                    Matrix_CCS_MtransVMult

IMPLICIT NONE


INTERFACE Calculate_Residual
    MODULE PROCEDURE Calculate_Residual_A
    MODULE PROCEDURE Calculate_Residual_B
END INTERFACE Calculate_Residual

CONTAINS
 !+202+################################################!
!                                                       !
!           Calculate_Residual_A                        !
!                                                       !
 !#####################################################!
SUBROUTINE Calculate_Residual_A( iter, iU )

INTEGER,        INTENT(IN)                      ::  iter
INTEGER,        INTENT(IN)                      ::  iU


INTEGER                                         ::  Convergence_Type = 3


INTEGER                                         ::  ui, l, m, lm_loc, map_loc, i

REAL(idp),      DIMENSION(1:3,1:LM_LENGTH,1:2)  ::  LOne_Norm
REAL(idp),      DIMENSION(1:3,1:LM_LENGTH,1:2)  ::  LTwo_Norm
REAL(idp),      DIMENSION(1:3,1:LM_LENGTH,1:2)  ::  LInf_Norm

REAL(idp),      DIMENSION(1:NUM_R_NODES)        ::  Work_Vec

LOne_Norm = 0.0_idp
LTwo_Norm = 0.0_idp
LInf_Norm = 0.0_idp


DO l = 0,L_LIMIT
DO m = -l,l

    lm_loc = Map_To_lm(l,m)
    
    CALL Matrix_CCS_MtransVMult(NUM_R_NODES,                    &
                                Factored_NNZ,                   &
                                Laplace_Factored_Row(:,l),      &
                                Laplace_Factored_COL(:,l),      &
                                FP_Iter_Matrix_Storage(:,l),    &
                                dVA_Coeff_Vector(:,lm_loc,iU),  &
                                Work_Vec                        )
           
!    PRINT*,"Work_Vec"
!    PRINT*,Work_Vec

    CALL Matrix_CCS_MVMult( NUM_R_NODES,                    &
                            Factored_NNz,                   &
                            Laplace_Factored_Row(:,l),      &
                            Laplace_Factored_COL(:,l),      &
                            FP_Iter_Matrix_Storage(:,l),    &
                            Work_Vec,                       &
                            FP_Laplace_Vector(:,lm_loc,iU)  )

!    PRINT*,"After MVMult"
!    print*,FP_Laplace_Vector(:,lm_loc,iU)
!    print*,"Load Storage"
!    Print*,FP_Iter_Load_Storage(:,lm_loc)
    
    IF ( ANY( abs(FP_Iter_Load_Storage(:,lm_loc)) == 0.0_idp ) ) THEN
        DO i = 1,Num_R_Nodes
            IF ( abs(FP_Iter_Load_Storage(i,lm_loc)) == 0.0_idp ) THEN
                FP_Residual_Vector(i,lm_loc,iU) = ( FP_Laplace_Vector(i,lm_loc,iU)  &
                                                - FP_Iter_Load_Storage(i,lm_loc) )
            ELSE
                FP_Residual_Vector(i,lm_loc,iU) = ( FP_Laplace_Vector(i,lm_loc,iU)  &
                                                - FP_Iter_Load_Storage(i,lm_loc) )  &
                                                / FP_Iter_Load_Storage(i,lm_loc)
            
            END IF
    
        END DO
    ELSE
        FP_Residual_Vector(:,lm_loc,iU) = ( FP_Laplace_Vector(:,lm_loc,iU)          &
                                        - FP_Iter_Load_Storage(:,lm_loc) )          &
                                        / FP_Iter_Load_Storage(:,lm_loc)
    END IF

!    PRINT*,"Residual Vector"
!    PRINT*,FP_Residual_Vector(:,lm_loc,iU)
!    CALL sleep(1)
    
!    PRINT*,"Laplace Vector"
!    PRINT*,FP_Laplace_Vector(:,lm_loc,iU)
!    CALL sleep(1)
!    PRINT*,"Load Vector"
!    PRINT*,FP_Iter_Load_Storage(:,lm_loc)


    Resid_Norms(1,lm_loc,iter,iU) = SUM(ABS(FP_Residual_Vector(:,lm_loc,iU) ) )

    Resid_Norms(2,lm_loc,iter,iU) = SQRT(REAL(DOT_PRODUCT(FP_Residual_Vector(:,lm_loc,iU),         &
                                    FP_Residual_Vector(:,lm_loc,iU) ), idp) )

    Resid_Norms(3,lm_loc,iter,iU) = MAXVAL( (/ Resid_Norms(3,lm_loc,iter,iU), ABS(FP_Residual_Vector(:,lm_loc,iU) ) /) )
!
!    PRINT*,"Residuals : ",Resid_Norms(1,lm_loc,iter,iU), &
!                          Resid_Norms(2,lm_loc,iter,iU), &
!                          Resid_Norms(3,lm_loc,iter,iU)

END DO ! m Loop
END DO ! l Loop





END SUBROUTINE Calculate_Residual_A




 !+202+################################################!
!                                                       !
!           Calculate_Residual_B                        !
!                                                       !
 !#####################################################!
SUBROUTINE Calculate_Residual_B(  Residual, iU, iVB )

REAL(idp),          INTENT(OUT)             ::  Residual
INTEGER,            INTENT(IN)              ::  iU
INTEGER,            INTENT(IN)              ::  iVB



END SUBROUTINE Calculate_Residual_B






 !+202+################################################!
!                                                       !
!          Clear_FP_Diagnostics_Tables                  !
!                                                       !
 !#####################################################!
SUBROUTINE Clear_FP_Diagnostics_Tables( Iter, iU )

INTEGER,        INTENT(IN)                      ::  iter
INTEGER,        INTENT(IN)                      ::  iU

Resid_Norms(:,:,iter,iU) = 0.0_idp
Update_Norms(:,iter,iU) = 0.0_idp

END SUBROUTINE Clear_FP_Diagnostics_Tables



END MODULE Diagnostics_Fixed_Point
