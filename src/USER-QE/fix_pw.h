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
#include <string>
#include <map>

namespace LAMMPS_NS {

class FixPW : public Fix {
 public:
  FixPW(class LAMMPS *, int, char **);
  virtual ~FixPW() override;
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
 protected:
  static bool initialized; // only one instance can run at a time

  int ec_group{-1}; // id of electrostatic coupling group
  bigint nqm{};    // number of QM atoms
  bigint nec{};    // number of EC atoms
  double escale{}, fscale{}; // scale factor for energies and forces
  std::string inp_file; // PWScf input file name
  std::string out_file; // PWScf output file name

  // keep track of atoms
  std::vector<tagint> qm_tags; // save tags to correctly assign atoms
  std::vector<tagint> ec_tags;
  std::map<tagint, int> qm_hash; // assign tags to their qm-id
  std::map<tagint, int> ec_hash;
  void collect_tags(int mask, bigint nat,
                    std::vector<tagint>& tags,
                    std::map<tagint, int>& hash);

  // link atoms
  struct linkgroup{
      tagint link_atom;
      tagint qm_atom;
      double ratio;
  };
  std::vector<linkgroup> linkatoms; // tags of link atoms

  // handle charges at link atoms
  enum class charge_t{Z1, Z2, Z3, RCD, CS};
  charge_t charge_scheme{charge_t::Z1};
  std::vector<std::pair<tagint, tagint>> chargeatoms; // tags of atoms with modified charges
  std::vector<tagint> get_bound_neighbors(tagint tag, int mask, int dist=0);

  // communication
  struct commdata_t{
    int tag;
    double x,y,z;
  };
  std::vector<commdata_t> comm_buf; // buffer for lammps-internal communication
  double** qm_buf{nullptr}; // buffer for communicating qm-atoms to pwscf
  double** ec_buf{nullptr}; // buffer for communicating ec-atoms to pscf
  int *recv_count_buf; // buffer for variable-size collectives
  int *displs_buf;
  void collect_positions(bigint nat, double** buffer, const std::map<tagint, int>& hash);
  void distribute_forces(bigint nat, double** buffer, const std::vector<tagint>& tags);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
