/*--------------------------------------------------------------------------
 * Header file for the Directions data structures
 *--------------------------------------------------------------------------*/

#ifndef KRIPKE_DIRECTIONS_H__
#define KRIPKE_DIRECTIONS_H__

#include <vector>
#include "charm++.h"

class Grid_Data;
struct Input_Variables;

/**
 * Contains information needed for one quadrature set direction.
 */
struct Directions{
public:
  Directions() {}
  ~Directions() {}

  double xcos;              /* Absolute value of the x-direction cosine. */
  double ycos;              /* Absolute value of the y-direction cosine. */
  double zcos;              /* Absolute value of the z-direction cosine. */
  double w;                 /* weight for the quadrature rule.*/
  int id;                   /* direction flag (= 1 if x-direction
                            cosine is positive; = -1 if not). */
  int jd;                   /* direction flag (= 1 if y-direction
                            cosine is positive; = -1 if not). */
  int kd;                   /* direction flag (= 1 if z-direction
                            cosine is positive; = -1 if not). */
  int octant;

  // Pack-UnPack method
  inline void pup(PUP::er &p) {
    p|xcos;
    p|ycos;
    p|zcos;
    p|w;
    p|id;
    p|jd;
    p|kd;
    p|octant;
  }
};

void InitDirections(Grid_Data *grid_data, Input_Variables *input_vars);

#endif
