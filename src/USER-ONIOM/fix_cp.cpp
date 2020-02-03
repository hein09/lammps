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

/* ----------------------------------------------------------------------
   Contributing author: Sebastian Gsänger (FAU)
------------------------------------------------------------------------- */

#include "fix_cp.h"

#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

extern "C" {
void start_cp(MPI_Fint comm, int partitions[4],
              const char *input_file, const char *output_file,
              int* nat, double *x_qm,
              bigint nec);
void update_cp(double *x_qm, double *x_ec);
void calc_cp(double *f_qm, double *f_ec);
void end_cp();
//double energy_cp(void);
}

/* ---------------------------------------------------------------------- */

FixCP::FixCP(LAMMPS *l, int narg, char **arg):
    FixQE{l, narg, arg}
{
  // initially collect coordinates so CP will be initialized from LAMMPS instead of file
  auto nat = static_cast<int>(qm_coll->nat);
  collect_positions();

  // Boot up CP
//  start_cp(MPI_Comm_c2f(world), mpi_partitions,
//           inp_file.data(), out_file.data(),
//           &nat, qm_buf[0],
//           nec);
  initialized = true;

  // check if CP has been launched succesfully and with compatible input
  /* NOTE: using one, not sure if a stray CP process may run off.
   * may be if CP had been compiled against a different MPI version.
   * if so, may be resolved if building CP is done via LAMMPS' buildsystem.
   */
  if (nat<0){
    error->one(FLERR, "Error opening output file for fix qe/cp.");
  }else if(nat != qm_coll->nat){
    error->one(FLERR, "Mismatching number of atoms in fix qe/cp.");
  }
}

FixCP::~FixCP()
{
  if(initialized){
//    end_cp();
    initialized = false;
  }
}

int FixCP::setmask()
{
  return POST_FORCE
       | MIN_POST_FORCE
       | POST_INTEGRATE
       | MIN_PRE_FORCE
       | THERMO_ENERGY
  ;
}

double FixCP::compute_scalar()
{
// TODO
//  return energy_cp() * escale;
  return 0;
}

void FixCP::setup(int i)
{
  post_force(i);
}

void FixCP::min_setup(int i)
{
  post_force(i);
}

void FixCP::min_post_force(int i)
{
  post_force(i);
}

void FixCP::post_force(int)
{
  // execute a CP-step and receive forces
//  calc_cp(qm_buf[0], ec_buf[0]);

  // distribute forces across LAMMPS processes
  distribute_forces();
}

void FixCP::min_pre_force(int)
{
  post_integrate();
}

void FixCP::post_integrate()
{
  // collect positions of LAMMPS processes
  collect_positions();

  // transmit to CP
//  update_cp(qm_buf[0], ec_buf[0]);
}
