/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ONIOM_CLASS

OniomStyle(pwscf,OniomPW)

#else

#ifndef LMP_ONIOM_PW_H
#define LMP_ONIOM_PW_H

#include "oniom.h"

namespace LAMMPS_NS {

class OniomPW: public Oniom{
public:
    OniomPW(LAMMPS *l, int narg, char **args);
    ~OniomPW() override;
    int get_capabilities() override;
};

}

#endif // LMP_ONIOM_PW_H
#endif
