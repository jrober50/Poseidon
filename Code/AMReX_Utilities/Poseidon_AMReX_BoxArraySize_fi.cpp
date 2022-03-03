
#ifdef POSEIDON_AMREX_FLAG
#include <AMReX_BoxArray.H>
#include <iostream>

using namespace amrex;
using namespace std;

extern "C"
{

    Long amrex_fi_boxarraysize(const BoxArray* BA  )
    {
        return BA->size();
        
    }

}

#endif

