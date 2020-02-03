/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(qe/cp,FixCP)

#else

#ifndef LMP_FIX_CP_H
#define LMP_FIX_CP_H

#include "fix_qe.h"

namespace LAMMPS_NS {

class FixCP : public FixQE {
 public:
  FixCP(class LAMMPS *, int, char **);
  ~FixCP() override;
  int setmask() override;
  // send new positions
  void post_integrate() override;
  void min_pre_force(int) override;
  // trigger calculation, receive forces
  void post_force(int) override;
  void min_post_force(int) override;
  void setup(int) override;
  void min_setup(int) override;
  // report internals
  double compute_scalar() override;
//  double memory_usage() override;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
