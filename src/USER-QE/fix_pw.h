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

FixStyle(qe/pw,FixPW)

#else

#ifndef LMP_FIX_PW_H
#define LMP_FIX_PW_H

#include "fix.h"
#include <vector>

namespace LAMMPS_NS {

class FixPW : public Fix {
 public:
  FixPW(class LAMMPS *, int, char **);
  virtual ~FixPW() override;
  int setmask() override;
  void post_force(int) override;
  void post_integrate() override;
//  void min_post_force(int);
//  double compute_scalar() override;
//  double memory_usage() override;
 protected:
  int npool{1}, ntask{1}, nband{1}, ndiag{1}; // PWScf partition flags
  bigint nqm;    // number of QM atoms
  double fscale; // scale factor for forces
  std::string inp_file; // PWScf input file name
  std::string out_file; // PWScf output file name
  std::vector<tagint> tags; // save tags to correctly assign atoms
  std::map<tagint, int> hash; // assign tags to their qm-id
  struct commdata_t{
    int tag;
    double x,y,z;
  };
  std::vector<commdata_t> comm_buf; // buffer for internal communication
  double** double_buf{nullptr}; // buffer for communicating to pwscf
  int *recv_count_buf;
  int *displs_buf;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
