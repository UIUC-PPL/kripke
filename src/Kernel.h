/*--------------------------------------------------------------------------
 * Header file for the Grid_Data data structures
 *--------------------------------------------------------------------------*/

#ifndef KRIPKE_KERNEL_H__
#define KRIPKE_KERNEL_H__

#include "Kripke.h"

/* Forward declarations */
struct Grid_Data;
struct SubTVec;
struct Subdomain;

/**
 * This is the Kernel base-class and interface definition.
 * This abstracts the storage of Psi, Phi, L, L+ from the rest of the code,
 * providing data-layout specific routines.
 */
class Kernel {
  public:
    virtual Nesting_Order nestingPsi(void) const = 0;
    virtual Nesting_Order nestingPhi(void) const = 0;
    virtual Nesting_Order nestingSigt(void) const = 0;
    virtual Nesting_Order nestingEll(void) const = 0;
    virtual Nesting_Order nestingEllPlus(void) const = 0;
    virtual Nesting_Order nestingSigs(void) const = 0;

    // Computational Kernels
    virtual void LTimes(Grid_Data *grid_data) = 0;
    virtual void LPlusTimes(Grid_Data *grid_data) = 0;
    virtual void scattering(Grid_Data *grid_data) = 0;
    virtual void source(Grid_Data *grid_data) = 0;
    virtual void sweep(Subdomain *ga_set) = 0;
    
    // Pack-UnPack method
    virtual void pup(PUP::er &p) {};
};


// Factory to create correct kernel object
Kernel *createKernel(Nesting_Order, int num_dims);

#endif
