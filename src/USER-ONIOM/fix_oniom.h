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

FixStyle(oniom,FixONIOM)

#else

#ifndef LMP_FIX_ONIOM
#define LMP_FIX_ONIOM

#include "fix.h"
#include <vector>

namespace LAMMPS_NS {

class FixONIOM : public Fix {
 public:
  FixONIOM(class LAMMPS *, int, char **);
  ~FixONIOM();
  int setmask();
  void init();

  // send up-to-date position information to QM and MM slave code
  void post_integrate();

  // receive and update forces
  void setup(int);
  void final_integrate();

  double memory_usage();

 protected:
  void exchange_positions();    // communicate positions to QM and MM slave
  void exchange_forces();       // collected forces from QM and MM slave

 protected:
  struct commdata_t{
      tagint tag;
      double x,y,z;
      double q;
  };
  struct conn_up_t{
      int partition;
      int mc_group;
      int mc_nat;
      int ec_group;
      int ec_nat;
      MPI_Comm comm;
      std::vector<tagint> tags;
      std::map<tagint,int> hash;
  };
  struct conn_down_t{
      int mc_group;
      int mc_nat;
      int ec_nat;
      std::vector<tagint> tags;
      std::map<tagint,int> hash;
      int part_low;
      int ec_low;
      MPI_Comm comm_low;
      int part_high;
      int ec_high;
      MPI_Comm comm_high;
  };

  conn_up_t *upwards{nullptr};
  std::map<std::string, conn_down_t> downwards{};

  double fscale;                    // scale factor for forces
  int    verbose{0};                // print level (<= 0 means no output)
  std::vector<commdata_t> comm_buf; // persistent storage for per-atom data
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not establish connection to oniom slave/master:

Mismatching partition IDs in master/slave 'fix oniom' calls.
Check the input script and your command line so all IDs match up correctly.

E: Could not find fix oniom master/slave group ID

Invalid group ID given for mechanical coupling (mandatory).

E: Could not find fix oniom master_ec/slave_ec group ID

Invalid group ID given for electrostatic coupling (optional environment atoms).

E: Fix oniom master cannot serve fix oniom slave_ec

When master only provides mechanical coupling,
slave must not request electrostatic coupling.

E: Electrostatic and mechanical group IDs must differ

Cannot use same group for both purposes

E: Fix oniom MC group must be static

Dynamic groups are only allowed for electrostatic coupling

W: Fix oniom should come before other integrating fixes.

If other integrating fixes are present, they shall come after fix oniom,
otherwise updated forces are not considered.

W: Using fix oniom slave with other integrating fixes.

On slave partitions, the effects of integrating fixes may be overwritten by
or interfering with fix oniom, at least if target groups overlap.

*/
