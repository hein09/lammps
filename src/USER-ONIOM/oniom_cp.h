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

OniomStyle(cp,OniomCP)

#else

#ifndef LMP_ONIOM_CP_H
#define LMP_ONIOM_CP_H

#include "oniom.h"

namespace LAMMPS_NS {

class OniomCP: public Oniom{
public:
    OniomCP(LAMMPS *l, int, char**);
    ~OniomCP() override;
    int get_capabilities() override;
};

}

#endif // LMP_ONIOM_CP_H
#endif
