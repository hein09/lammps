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

#ifndef LMP_FIX_QE_H
#define LMP_FIX_QE_H

#include "fix_collect.h"

#include <vector>
#include <map>

namespace LAMMPS_NS {

class FixQE : public FixCollect {
 public:
  FixQE(LAMMPS *, int, char **);
 protected:
  // only one instance can run at a time per partition
  static bool initialized;

  // QE-MPI levels (npool, ntask, nband, ndiag)
  int mpi_partitions[4]{1,1,1,1};

  // input/output file name for respective QE-code
  std::string inp_file{}, out_file{};

  // map atoms between LAMMPS-tags and QE-index
  std::shared_ptr<collection_t> qm_coll{nullptr};
  std::shared_ptr<collection_t> ec_coll{nullptr};

  // communicate positions, forces and energies back and forth
  double escale{}, fscale{}; // scale factor for energies and forces
  void collect_positions();
  void distribute_forces();

  // link atoms
  struct linkgroup{
    tagint link_atom;
    tagint qm_atom;
    double ratio;
  };
  std::vector<linkgroup> linkatoms{};

  // charge-scheme for link-atoms
  enum class charge_t{Z1, Z2, Z3, RCD, CS};
  charge_t charge_scheme{charge_t::Z1};
  std::vector<std::pair<tagint, tagint>> chargeatoms; // tags of atoms with modified charges
  std::vector<tagint> get_bound_neighbors(tagint tag, int groupmask, int dist);
 private:
  void distribute_forces_impl(collection_t& coll);
};

}

#endif // LMP_FIX_QE_H

/* ERROR/WARNING messages:

*/
