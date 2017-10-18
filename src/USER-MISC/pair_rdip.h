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

#ifdef PAIR_CLASS

PairStyle(rdip,PairRDIP)

#else

#ifndef LMP_PAIR_RDIP_H
#define LMP_PAIR_RDIP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairRDIP : public Pair {
  public:
   PairRDIP(class LAMMPS *);
   virtual ~PairRDIP();

   virtual void compute(int, int);
   void settings(int,char**);
   void coeff(int,char**);
   void init_style();
   double init_one(int, int);
   void read_restart(FILE *);
   void write_restart(FILE *);
   void write_restart_settings(FILE *);
   void read_restart_settings(FILE *);

  private:
   int nneigh;            // # of nearest neighbors
   int *map;              // # map of internal atom types
   double cutoff;

   // parameters for potential
   double A, C, C0, C2, C4, delta, z0, lambda;

   void read_file(char *);
   void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

*/
