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

FixStyle(ONIOM,FixONIOM)

#else

#ifndef LMP_FIX_ONIOM
#define LMP_FIX_ONIOM

#include "fix.h"
#include <vector>
#include <map>

namespace LAMMPS_NS {

class FixONIOM : public Fix {
 public:
  FixONIOM(class LAMMPS *, int, char **);
  ~FixONIOM() override;
  int setmask() override;
  void init() override;

  void post_integrate() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void min_pre_force(int) override;
  void min_post_force(int) override;
  void end_of_step() override;
  void post_run() override;

  double memory_usage() override;

 protected:
  void send_positions();        // send positions from master to slaves
  void receive_positions();     // receive positions from master
  void receive_forces();        // receive forces from slaves
  void send_forces();           // send forces from slave to master

 public:
  bool   master{false};         // identifies first == master partition
  bool   minimize{false};       // distinguish between MD and Minimization
  bool   run_once{false};       // make sure that Min and MD match communication
  int    verbose{0};            // print level (<= 0 means no output)
  enum {PLUS=0x0, MINUS=0x1};
  struct conn_t{
    conn_t(int mode=-1, int mc_group=-1, int mc_nat=-1, MPI_Comm comm=MPI_COMM_NULL,
           const std::vector<tagint>& tags={}, const std::map<tagint, int>& hash={})
        : mode{mode}, mc_group{mc_group}, mc_nat{mc_nat}, comm{comm},
          tags{tags}, hash{hash}
    {}
    int mode;
    int mc_group;
    int mc_nat;
    MPI_Comm comm;
    std::vector<tagint> tags;
    std::map<tagint,int> hash;
  };
  std::vector<conn_t> connections{};

 protected:
  struct commdata_t{
      tagint tag;
      double x,y,z;
      double q;
  };

  std::vector<commdata_t> comm_buf; // persistent storage for per-atom data
  int *recv_count_buf;
  int *displs_buf;
};

}

#endif
#endif
