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

#ifndef LMP_FIX_CONNECT
#define LMP_FIX_CONNECT

#include "fix.h"

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace LAMMPS_NS {

class FixCollect : public Fix {
 public:
  FixCollect(LAMMPS *, int, char **);

  double memory_usage() override;
  // data coupled to its atom's tag

  struct commdata_t{
      int idx;
      double x[3];
  };
 protected:
  /* container for a collection
   *
   * stores forward and backwards mapping between LAMMPS-tags
   * and indices local to this or shared-ownership fixes
   */
  struct collection_t{
      bigint nat{};
      int group_id{-1};
      std::vector<tagint> idx2tag{}; // fix to lammps
      std::map<tagint, int> tag2idx{}; // lammps to fix
      std::vector<commdata_t> buffer{}; // data-storage
  };

  // return a collection for the given group
  // will be shared between multiple callers
  std::shared_ptr<collection_t> get_collection(char *group);

  // collect data from lammps on partition's root
  void collect_lammps(collection_t& coll, double** source);

 private:
  std::vector<int> recv_count_buf, displs_buf; // buffer for variable-size MPI collective strides
  void collect_tags(int groupmask, collection_t& coll);
  // tracking-table for collections
  static std::map<std::string, std::weak_ptr<collection_t>> collections;
};

}

#endif // LMP_FIX_CONNECT

/* ERROR/WARNING messages:

*/
