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

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "universe.h"
#include "update.h"

#include <algorithm>

using namespace LAMMPS_NS;
using namespace FixConst;

extern "C" {
void start_pw(MPI_Fint comm, int nimage, int npool, int ntaskgroup, int nband, int ndiag,
              const char *input_file, const char *output_file,
              int* nat);
//void step_pw(int *exit_status);
void update_pw(double *x, int root);
void calc_pw(double *f, int root);
void end_pw(int *exit_status);
}

/* ---------------------------------------------------------------------- */

FixPW::FixPW(LAMMPS *l, int narg, char **arg):
    Fix{l, narg, arg}
{
  if (narg < 4) error->all(FLERR, "Illegal fix qe/pw command");

  if (atom->tag_enable == 0)
    error->all(FLERR,"fix qe/pw requires atom IDs");

  // ensure that only one fix is active per partition
  if (modify->find_fix_by_style("qe/pw") != -1)
    error->all(FLERR, "Only one instance of fix qe/pw allowed at a time");

  // make sure we understand eachother
  if (strcmp(update->unit_style, "metal") == 0){
    fscale = 1.0/23.0609;
  }else if (strcmp(update->unit_style, "real") == 0){
    fscale = 1.0;
  }else error->all(FLERR, "fix qe/pw requires real or metal units");

  // save file-name for later
  inp_file = arg[3];

  // parse additional arguments
  int iarg = 4;
  while(iarg < narg-1){
    if(strcmp(arg[iarg], "npool") == 0){
      npool = force->inumeric(FLERR, arg[iarg+1]);
      iarg += 2;
    }else if(strcmp(arg[iarg], "ntask") == 0){
      ntask = force->inumeric(FLERR, arg[iarg+1]);
      iarg += 2;
    }else if(strcmp(arg[iarg], "nband") == 0){
      nband = force->inumeric(FLERR, arg[iarg+1]);
      iarg += 2;
    }else if(strcmp(arg[iarg], "ndiag") == 0){
      ndiag = force->inumeric(FLERR, arg[iarg+1]);
      iarg += 2;
    }
  }

  // redirect stdout when PWScf is running
  out_file = "log.pwscf";
  if(universe->existflag){
      out_file += '.' + std::to_string(universe->iworld);
  }

  // initialize PWScf, receive number of qm-atoms
  int nat{0};
  nqm = group->count(igroup);
  if(nqm > MAXSMALLINT){
      error->all(FLERR, "Too many QM atoms for fix qe/pw.");
  }

  start_pw(MPI_Comm_c2f(world), 1, npool, ntask, nband, ndiag, inp_file.data(), out_file.data(), &nat);

  // check if pw has been launched succesfully and with compatible input
  /* NOTE: using one, not sure if a stray pwscf process may run off.
   * may be if pw had been compiled against a different MPI version.
   * if so, may be resolved if building pw is done via lammps' buildsystem.
   */
  if (nat<0){
      error->one(FLERR, "Error opening output file for fix qe/pw.");
  }else if( nat != nqm){
      error->one(FLERR, "Mismatching number of atoms in fix qe/pw.");
  }

  // collect and save tags of main group
  recv_count_buf = new int[universe->nprocs];
  displs_buf = new int[universe->nprocs];

  // collect process-local tags
  decltype (tags) tmp{};
  tmp.reserve(static_cast<size_t>(nqm));
  for(int i=0; i<atom->nlocal; ++i){
    if(atom->mask[i] & groupbit)
      tmp.push_back(atom->tag[i]);
  }

  // gather number of tags
  int send_count = tmp.size() * sizeof(decltype (tmp)::value_type);
  MPI_Allgather(&send_count, 1, MPI_INT, recv_count_buf, 1, MPI_INT, world);

  // construct displacement
  displs_buf[0] = 0;
  for(int i=1; i<universe->nprocs; ++i){
      displs_buf[i] = displs_buf[i-1]+recv_count_buf[i-1];
  }

  // gather and sort tags
  tags.resize(static_cast<size_t>(nqm));
  MPI_Allgatherv(tmp.data(), send_count, MPI_BYTE,
                 tags.data(), recv_count_buf, displs_buf, MPI_BYTE, world);
  std::sort(tags.begin(), tags.end());

  // assign tags to qm-id
  for(int i=0; i<nqm; ++i){
    hash[tags[i]] = i;
  }
}

FixPW::~FixPW()
{
  int result;
  end_pw(&result);

  delete [] recv_count_buf;
  delete [] displs_buf;
  memory->destroy(double_buf);
}

int FixPW::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_INTEGRATE;
//  mask |= THERMO_ENERGY;
  return mask;
}

// TODO: collective communications may be optimized (away?)

void FixPW::post_force(int)
{
  // run scf and receive forces
  calc_pw(double_buf[0], static_cast<int>(comm->me == 0));
  // communicate forces to other lammps processes
  if(comm->me == 0){
    for(int i=0; i<nqm; ++i){
      comm_buf[i] = {i, double_buf[i][0],
                     double_buf[i][1], double_buf[i][2]};
    }
  }
  MPI_Bcast(comm_buf.data(), nqm*sizeof(decltype(comm_buf)::value_type), MPI_BYTE, 0, world);
  // add forces
  const tagint * const tag = atom->tag;
  double ** f = atom->f;
  for(int i=0; i<atom->nlocal; ++i){
    for(auto& dat: comm_buf){
      if(tag[i] == dat.tag){
        f[i][0] += dat.x*fscale;
        f[i][1] += dat.y*fscale;
        f[i][2] += dat.z*fscale;
      }
    }
  }
}

void FixPW::post_integrate()
{
  // collect coordinates
  const int *const mask = atom->mask;
  const tagint * const tag = atom->tag;
  const double * const * const x = atom->x;

  decltype(comm_buf) tmp;
  tmp.reserve(static_cast<size_t>(nqm));
  for(int i=0; i<atom->nlocal; ++i){
    if(mask[i] & groupbit){
      tmp.push_back({hash.at(tag[i]), x[i][0], x[i][1], x[i][2]});
    }
  }
  int send_count = tmp.size() * sizeof(decltype (tmp)::value_type);
  // gather number of local atoms
  MPI_Allgather(&send_count, 1, MPI_INT, recv_count_buf, 1, MPI_INT, world);
  // construct displacement
  displs_buf[0] = 0;
  for(int i=1; i<universe->nprocs; ++i){
    displs_buf[i] = displs_buf[i-1]+recv_count_buf[i-1];
  }
  // gather comm_buf
  comm_buf.resize(static_cast<size_t>(nqm));
  MPI_Allgatherv(tmp.data(), send_count, MPI_BYTE,
                 comm_buf.data(), recv_count_buf, displs_buf, MPI_BYTE, world);

  // extract data to contiguous memory
  memory->grow(double_buf, static_cast<int>(nqm), 3, "pw/qe:double_buf");
  for(auto& dat: comm_buf){
    double_buf[dat.tag][0] = dat.x;
    double_buf[dat.tag][1] = dat.y;
    double_buf[dat.tag][2] = dat.z;
  }

  MPI_Barrier(world);
  // transmit to PWScf
  update_pw(double_buf[0], static_cast<int>(comm->me == 0));
}
