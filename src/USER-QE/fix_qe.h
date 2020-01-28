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

#include "fix.h"

#include <vector>
#include <map>

namespace LAMMPS_NS {

class FixQE : public Fix {
 public:
  FixQE(LAMMPS *, int, char **);
  ~FixQE() override;
 protected:
  // only one instance can run at a time per partition
  static bool initialized;

  // id of electrostatic coupling group
  int ec_group{-1};

  // QE-MPI levels (npool, ntask, nband, ndiag)
  int mpi_partitions[4]{1,1,1,1};

  // input/output file name for respective QE-code
  std::string inp_file{}, out_file{};

  // map atoms between LAMMPS-tags and QE-index
  bigint nqm{}, nec{}; // number of QM and EC atoms
  std::vector<tagint> qm_idx2tag{};
  std::vector<tagint> ec_idx2tag{};
  std::map<tagint, int> qm_tag2idx{};
  std::map<tagint, int> ec_tag2idx{};
  void collect_tags(int groupmask, bigint nat,
                    std::vector<tagint>& idx2tag,
                    std::map<tagint, int>& tag2idx);

  // communicate positions, forces and energies back and forth
  double escale{}, fscale{}; // scale factor for energies and forces
  void collect_positions();
  void distribute_forces();

  // communication structures
  struct commdata_t{
    int tag;
    double x,y,z;
  };
  std::vector<commdata_t> comm_buf{}; // buffer for LAMMPS-internal communication
  double** qm_buf{nullptr}; // buffer for communicating qm-atoms to QE
  double** ec_buf{nullptr}; // buffer for communicating ec-atoms to QE
  int *recv_count_buf, *displs_buf; // buffer for variable-size MPI collectives

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
  std::vector<tagint> get_bound_neighbors(tagint tag, int groupmask, int dist=0);
 private:
  void collect_positions_impl(bigint nat, double** buffer, const std::map<tagint, int>& tag2idx);
  void distribute_forces_impl(bigint nat, double** buffer, const std::vector<tagint>& idx2tag);
};

}

#endif // LMP_FIX_QE_H

/* ERROR/WARNING messages:

*/
