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
    double cutoff;   //cutoff for interaction
    int *tmap;       // map lammps to rdip types
    struct RDIPType{
      int nneigh;
      double cut_neigh, delta;
    };
    RDIPType *ntmap; // type information

    struct RDIPParam{
      double A, C, C0, C2, C4, z0, lambda;
    };
    RDIPParam **pmap; // pair interaction information

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
