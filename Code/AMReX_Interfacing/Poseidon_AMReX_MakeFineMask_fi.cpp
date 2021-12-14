
#ifdef POSEIDON_AMREX_FLAG
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Geometry.H>

using namespace amrex;


extern "C"
{

    void amrex_fi_makefinemask(iMultiFab*& Mask,
                               const BoxArray& Coarse_Box,
                               const DistributionMapping& Coarse_DM,
                               const BoxArray& Fine_Box,
                               int C_coarse, int C_fine )
    {
        Mask = new iMultiFab( makeFineMask(Coarse_Box, Coarse_DM, Fine_Box,
                                           amrex::IntVect(2), C_coarse, C_fine )            );
    }

}

#endif
