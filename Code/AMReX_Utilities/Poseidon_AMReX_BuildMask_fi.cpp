
#ifdef POSEIDON_AMREX_FLAG
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Geometry.H>
#include <AMReX_FabArray.H>

using namespace amrex;


extern "C"
{

    void amrex_fi_buildmask( iMultiFab*& mask,
                             Geometry const* geom,
                               int C_Covered,
                               int C_NotCovered,
                               int C_PhysBnd,
                               int C_Interior     )
    {
        mask->BuildMask( geom->Domain(),
                        geom->periodicity(),
                         C_Covered,
                         C_NotCovered,
                         C_PhysBnd,
                         C_Interior );
    }

}


#endif

