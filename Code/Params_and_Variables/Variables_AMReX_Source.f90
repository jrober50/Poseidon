   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Variables_AMReX_Source                                                       !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

USE iso_c_binding

USE Poseidon_Kinds_Module,  &
            ONLY : idp


IMPLICIT NONE

INTEGER(c_int), PARAMETER                               ::  iTrunk = 0
INTEGER(c_int), PARAMETER                               ::  iLeaf  = 1

INTEGER(c_int), PARAMETER                               ::  iCovered    = 2
INTEGER(c_int), PARAMETER                               ::  iNotCovered = 3
INTEGER(c_int), PARAMETER                               ::  iOutside    = 4
INTEGER(c_int), PARAMETER                               ::  iInterior   = 5


REAL(idp), CONTIGUOUS, POINTER                          :: Source_PTR(:,:,:,:)
INTEGER,   CONTIGUOUS, POINTER                          :: Mask_PTR(:,:,:,:)
INTEGER,   CONTIGUOUS, POINTER                          :: Ghost_PTR(:,:,:,:)


END MODULE Variables_AMReX_Source





