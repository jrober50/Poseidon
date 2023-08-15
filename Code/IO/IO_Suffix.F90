   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE IO_Suffix_Module                                                      !##!
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
                    L_Limit
            
USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    Num_T_Elements,         &
                    Num_P_Elements
            
USE Variables_FP, &
            ONLY :  FP_Anderson_M

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Max_Grid_Size,    &
                    AMReX_MaxLevel,        &
                    AMReX_Num_Levels
                    
                    
CONTAINS


 !+102+################################################!
!                                                       !
!       Create_Suffix                                   !
!                                                       !
 !#####################################################!
SUBROUTINE Create_Suffix( Suffix_Return,            &
                          Suffix_Flag_Option,       &
                          Frame_Option,             &
                          Suffix_Param_Type_Option, &
                          Suffix_Tail_Option      )

CHARACTER(LEN=40),  INTENT(OUT),    OPTIONAL    ::  Suffix_Return
CHARACTER(LEN=10),  INTENT(IN),     OPTIONAL    ::  Suffix_Flag_Option
INTEGER,            INTENT(IN),     OPTIONAL    ::  Frame_Option
INTEGER,            INTENT(IN),     OPTIONAL    ::  Suffix_Param_Type_Option

CHARACTER(LEN=4),   INTENT(IN),     OPTIONAL    ::  Suffix_Tail_Option

333 FORMAT('RE',I5.5,'_D',I2.2,'_L',I2.2)
444 FORMAT('RE',I5.5,'_TE',I3.3,'_D',I2.2,'_L',I2.2)
555 FORMAT('RE',I5.5,'_LVL',I3.3,'_D',I2.2,'_M',I2.2)

IF ( PRESENT(Suffix_Flag_Option) ) THEN

    !===========================!
    !                           !
    !           Params          !
    !                           !
    !===========================!
    IF ( Suffix_Flag_Option == "Params") THEN
    
    
        IF ( PRESENT(Suffix_Param_Type_Option) ) THEN
        
            IF ( Suffix_Param_Type_Option == 1 ) THEN
                WRITE(Suffix_Return,444)  Num_R_Elements,         &
                                        Num_T_Elements,         &
                                        Degree,                 &
                                        L_Limit
            ELSEIF ( Suffix_Param_Type_Option == 2 ) THEN
                WRITE(Suffix_Return,555)  Num_R_Elements,         &
                                        AMReX_MaxLevel,        &
                                        Degree,                 &
                                        FP_Anderson_M
            ELSE
                WRITE(Suffix_Return,333)  Num_R_Elements,         &
                                        Degree,                 &
                                        L_Limit
            END IF
        
        ELSE
            ! Default Suffix
            WRITE(Suffix_Return,333)  Num_R_Elements,         &
                                    Degree,                 &
                                    L_Limit
        END IF


    !===========================!
    !                           !
    !           Frame           !
    !                           !
    !===========================!
    ELSEIF ( SUffix_Flag_Option == "Frame") THEN
        
        IF ( PRESENT(Frame_Option) ) THEN
            WRITE(Suffix_Return,'(I5.5)') Frame_Option
        ELSE
            WRITE(Suffix_Return,'(I5.5)') 1
        END IF
    END IF
    
    
ELSE

    WRITE(Suffix_Return,'(I5.5)') 2
    
END IF





IF ( PRESENT(Suffix_Tail_Option) ) THEN
    WRITE(Suffix_Return,'(A,A,A)') TRIM(Suffix_Return),"_",TRIM(Suffix_Tail_Option)
END IF



END SUBROUTINE Create_Suffix






END MODULE IO_Suffix_Module
