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

#include "fix_pw.h"

#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

extern "C" {
void start_pw(MPI_Fint comm, int partitions[4],
              const char *input_file, const char *output_file,
              int* nat, FixPW::commdata_t *dat);
//              bigint nec);
void update_pw(FixPW::commdata_t *x_qm);//, double *x_ec);
void calc_pw(FixPW::commdata_t *f_qm);//, double *f_ec);
void end_pw(int *exit_status);
double energy_pw(void);
}

/* ---------------------------------------------------------------------- */

FixPW::FixPW(LAMMPS *l, int narg, char **arg):
    FixQE{l, narg, arg}
{
  // initially collect coordinates so PW will be initialized from LAMMPS instead of file
  auto nat = static_cast<int>(qm_coll->nat);
  collect_positions();

  // Boot up PWScf
  start_pw(MPI_Comm_c2f(world), mpi_partitions,
           inp_file.data(), out_file.data(),
           &nat, qm_coll->buffer.data());
  initialized = true;

  // check if pw has been launched succesfully and with compatible input
  /* NOTE: using one, not sure if a stray pwscf process may run off.
   * may be if pw had been compiled against a different MPI version.
   * if so, may be resolved if building pw is done via lammps' buildsystem.
   */
  if (nat<0){
    error->one(FLERR, "Error opening output file for fix qe/pw.");
  }else if(nat != qm_coll->nat){
    char msg[50];
    sprintf(msg, "Mismatching number of atoms in fix qe/pw: %d(pw) vs %d(lmp)",
            nat, qm_coll->nat);
    error->one(FLERR, msg);
  }
}

FixPW::~FixPW()
{
  if(initialized){
    int result;
    end_pw(&result);
    initialized = false;
  }
}

int FixPW::setmask()
{
  return POST_FORCE
       | MIN_POST_FORCE
       | POST_INTEGRATE
       | MIN_PRE_FORCE
       | THERMO_ENERGY
  ;
}

double FixPW::compute_scalar()
{
  return energy_pw() * escale;
}

void FixPW::setup(int i)
{
  post_force(i);
}

void FixPW::min_setup(int i)
{
  post_force(i);
}

void FixPW::min_post_force(int i)
{
  post_force(i);
}

void FixPW::post_force(int)
{
  // run scf and receive forces
  calc_pw(qm_coll->buffer.data());

  // distribute forces across LAMMPS processes
  distribute_forces();
}

void FixPW::min_pre_force(int)
{
  post_integrate();
}

void FixPW::post_integrate()
{
  // collect positions of LAMMPS processes
  collect_positions();

  // transmit to PWScf
  update_pw(qm_coll->buffer.data());
}
