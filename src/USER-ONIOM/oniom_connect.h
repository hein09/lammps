/* -*- c++ -*- -------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(oniom_connect,OniomConnect)

#else

#ifndef LMP_ONIOM_CONNECT_H
#define LMP_ONIOM_CONNECT_H

#include "pointers.h"

namespace LAMMPS_NS {

class OniomConnect : protected Pointers {
public:
    OniomConnect(LAMMPS *l);
    ~OniomConnect();
    void command(int, char **);
};

}

#endif // LMP_ONIOM_CONNECT_H
#endif

/* ERROR/WARNING messages:

TODO

*/
