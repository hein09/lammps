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

#ifdef COMMAND_CLASS

CommandStyle(oniom,Oniom)

#else

#ifndef LMP_ONIOM_H
#define LMP_ONIOM_H

#include "pointers.h"
#include "fix_oniom.h"

namespace LAMMPS_NS {

class Oniom : protected Pointers {
 public:
  Oniom(class LAMMPS *l): Pointers{l} {};
  void command(int, char **);
 private:
  void run(int, char **, FixONIOM*);
  void minimize(int, char **, FixONIOM*);
  FixONIOM* mkFix(int, char **);
};

}

#endif
#endif


/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid partition in oniom command

Make sure that all slave partitions are explicitely and correctly
assigned to ONIOM groups.

E: Could not find oniom ... group ID

The specified group name does not exist in this partition

E: Too many ... atoms for ONIOM calculation

Self-explanatory.

E: Could not establish connection to oniom slave/master

A slave partition may be not assigned to an ONIOM group.
May also be triggered by MPI problems.

E: Same partition used in multiple oniom regions

One slave partition has been assigned multiple times

E: Integrating fixes present in slave partition, ...

On slave partitions, integration will be performed by the oniom command.
Usually, no integrating fixes like NVE should be used.

E: Inconsistent number of ... atoms

The number of atoms differs between master and slave.
Ensure that the groups definition does not change in an unexpected manner.


*/
