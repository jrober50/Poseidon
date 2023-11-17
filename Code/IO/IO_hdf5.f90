   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE IO_hdf5_Module                                                        !##!
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
            ONLY : idp

USE Poseidon_Message_Routines_Module, &
            ONLY :  Warning_Message

USE Poseidon_Numbers_Module, &
            ONLY : pi

USE Poseidon_Units_Module, &
            ONLY :  Centimeter,     &
                    Shift_Units

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
            ONLY :  Degree,                 &
                    L_Limit,                &
                    Eq_Flags
                
USE Variables_MPI, &
            ONLY :  myID_Poseidon,          &
                    MasterID_Poseidon,      &
                    nProcs_Poseidon,        &
                    Poseidon_Comm_World

USE Variables_Derived, &
            ONLY :  Num_R_Nodes


USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    Num_T_Elements,         &
                    Num_P_Elements,         &
                    rlocs,                  &
                    tlocs,                  &
                    plocs,                  &
                    R_Inner,                &
                    R_Outer

USE Variables_IO, &
            ONLY :  Write_Flags,            &
                    Report_IDs,             &
                    Write_Results_R_Samps,  &
                    Write_Results_T_Samps,  &
                    Write_Results_P_Samps,  &
                    Total_Run_Iters,        &
                    Iter_Report_Num_Samples,&
                    Iter_Time_Table,        &
                    File_Suffix,            &
                    iWF_Source,             &
                    iWF_Results,            &
                    iRF_Run,                &
                    iRF_Frame,              &
                    iRF_Time,               &
                    iRF_Iter


USE Variables_External, &
            ONLY :  SelfSim_T

USE Return_Functions_FP,   &
            ONLY :  Calc_FP_Values_At_Location
                    
USE Variables_Interface, &
            ONLY :  Caller_NQ,                      &
                    Caller_Quad_DOF,                &
                    Caller_xL,                      &
                    Caller_RQ_xlocs,                &
                    Caller_TQ_xlocs,                &
                    Caller_PQ_xlocs
                
USE Variables_Tables, &
            ONLY :  Level_dx,                   &
                    Level_Ratios

USE Poseidon_Return_Routines_Module, &
            ONLY :  Calc_Var_at_Location

USE Return_Functions_FP , &
            ONLY :  Calc_Drv_At_Location_Type_B
            
USE Poseidon_Return_Routines_All, &
            ONLY :  Poseidon_Return_All_AMReX_Caller

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature_Locations,     &
                    Initialize_LGL_Quadrature_Locations


USE Functions_Mesh, &
            ONLY :  Create_Logarithmic_1D_Mesh,             &
                    Create_Uniform_1D_Mesh

USE Poseidon_IO_Parameters, &
            ONLY :  Poseidon_Output_Dir,        &
                    CFA_ShortVars,              &
                    Kij_ShortVars

USE Poseidon_File_Routines_Module, &
            ONLY :  Open_New_File,                  &
                    Open_Existing_File_Append

USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,               &
                    iPF_IO_Write_Results

USE Flags_Core_Module, &
            ONLY :  iPF_Core_Flags,         &
                    iPF_Core_Method_Mode,   &
                    iPF_Core_Method_Newtonian
                    
                    
#ifdef POSEIDON_AMREX_FLAG
use amrex_base_module

USE amrex_box_module,   &
            ONLY :  amrex_box

USE amrex_boxarray_module, &
            ONLY :  amrex_boxarray

use amrex_fort_module, &
            ONLY :  amrex_spacedim

USE amrex_multifab_module,  &
            ONLY :  amrex_multifab,         &
                    amrex_multifab_build,   &
                    amrex_imultifab_build,  &
                    amrex_imultifab_destroy

USE Variables_AMReX_Core, &
            ONLY :  MF_Source,              &
                    AMReX_Num_Levels,       &
                    Source_PTR,             &
                    Mask_PTR
                    
USE Parameters_AMReX, &
            ONLY :  iLeaf,                &
                    iTrunk

USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask

USE Poseidon_AMReX_Multilayer_Utilities_Module, &
            ONLY :  Find_Coarsest_Parent

#endif

IMPLICIT NONE



CONTAINS





 !+201+############################################################!
!                                                                   !
!        Calculate_Max_Difference                                   !
!                                                                   !
 !#################################################################!
SUBROUTINE IO_HDF5_Write_Mesh_Fields()

CHARACTER(LEN = 500)                    ::  Filename



END SUBROUTINE IO_HDF5_Write_Mesh_Fields









END MODULE IO_hdf5_Module
