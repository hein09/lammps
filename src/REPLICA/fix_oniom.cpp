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
   Contributing author:  Sebastian Gsänger (FAU), Axel Kohlmeyer (ICTP)
------------------------------------------------------------------------- */

#include <algorithm>
#include <unistd.h>

#include "fix_oniom.h"

#include "atom.h"
#include "comm.h"
#include "update.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "universe.h"
#include "modify.h"
#include "timer.h"

using namespace LAMMPS_NS;
using namespace FixConst;

FixONIOM::FixONIOM(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  master = universe->iworld == 0;

  if (narg < 3)
    error->all(FLERR,"Illegal fix oniom command");

  if (atom->tag_enable == 0)
    error->all(FLERR,"Oniom requires atom IDs");

  if (!atom->q_flag)
    error->all(FLERR,"Oniom requires atom attribute q");

  if (atom->tag_consecutive() == 0)
    error->all(FLERR,"Oniom requires consecutive atom IDs");

  recv_count_buf = new int[universe->nprocs];
  displs_buf = new int[universe->nprocs];
  // parse arguments
  int iarg=3;
  int iworld = universe->iworld;
  while(iarg < narg){

    /********************************************************
     * Needs a static group of atoms and two slave partitions
     ********************************************************/
    if(strcmp(arg[iarg], "group") == 0){
      if(iarg+3 > narg)
        error->universe_all(FLERR, "Illegal oniom command");

      // arg[iarg+1] is group
      int mc_group = group->find(arg[iarg+1]);
      if(mc_group == -1)
        error->all(FLERR,"Could not find oniom group ID");
      if(group->dynamic[mc_group])
        error->all(FLERR,"MC group must be static");
      bigint temp_mc = group->count(mc_group);
      if(temp_mc > MAXSMALLINT){
          error->all(FLERR,"Too many MC atoms for ONIOM calculation");
      }
      int mc_nat = static_cast<int>(temp_mc);

      // arg[iarg+2] is low-level slave partition
      int low = force->inumeric(FLERR, arg[iarg+2])-1;
      if((low < 1 ) || (low>=universe->nworlds)){
        error->universe_all(FLERR, "Invalid partition in oniom command");
      }

      // arg[iarg+3] is high-level slave partition
      int high = force->inumeric(FLERR, arg[iarg+3])-1;
      if((high < 1 ) || (high>=universe->nworlds)){
        error->universe_all(FLERR, "Invalid partition in oniom command");
      }

      // push information to fix and create intercommunicators
      if(!iworld){
        // master partition
        auto setupMPI = [&](int part, MPI_Comm& comm){
          int me, flag{0};
          MPI_Request req;
          MPI_Comm_rank(world, &me);
          if(me==0){
            MPI_Isend(&part, 1, MPI_INT,
                      universe->root_proc[part],
                      part, universe->uworld, &req);
            // wait for one minute to establish connection
            int j{0};
            for(; j<60; ++j){
                MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
                if(flag) break;
                usleep(1000000);
            }
            if(j >= 60) error->universe_one(FLERR, "Could not establish connection to oniom slave");
          }
          MPI_Intercomm_create(world, 0, universe->uworld,
                               universe->root_proc[part], part, &comm);
        };
        int lowmode = FixONIOM::MINUS;
        int highmode = FixONIOM::PLUS;
        connections.push_back({lowmode, mc_group, mc_nat});
        setupMPI(low, connections.back().comm);
        connections.push_back({highmode, mc_group, mc_nat});
        setupMPI(high, connections.back().comm);
      }else if((iworld == (low)) || (iworld == (high))){
        // slave partitions
        if(!connections.empty()){
          error->universe_one(FLERR, "Same partition used in multiple oniom regions");
        }
        connections.push_back({0, mc_group, mc_nat});
        int me, flag{0}, tag;
        MPI_Request req;
        MPI_Comm_rank(world, &me);
        if(me==0){
          MPI_Irecv(&tag, 1, MPI_INT,
                    universe->root_proc[0],
                    universe->iworld, universe->uworld, &req);
          // wait for one minute to establish connection
          int j;
          for(j=0; j<60; ++j){
            MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
            if(flag) break;
            usleep(1000000);
          }
          if(j>=60) error->universe_one(FLERR, "Could not establish connection to oniom master");
        }
        MPI_Intercomm_create(world, 0, universe->uworld,
                             universe->root_proc[0], iworld,
                             &connections.back().comm);
      }
      iarg += 4;

    /********************
     * Verbosity
     ********************/
    }else if(strcmp(arg[iarg], "verbose") == 0){
      verbose = force->inumeric(FLERR, arg[iarg+1]);
      iarg += 2;
    }else{
      error->universe_all(FLERR, "Illegal fix oniom command");
    }
  }
}

/*********************************
 * Clean up on deleting the fix. *
 *********************************/
FixONIOM::~FixONIOM()
{
    for(auto& con: connections){
        MPI_Comm_free(&con.comm);
    }
    delete [] recv_count_buf;
    delete [] displs_buf;
}

/* ---------------------------------------------------------------------- */
int FixONIOM::setmask()
{
  return POST_INTEGRATE
       | POST_FORCE
       | MIN_PRE_FORCE
       | MIN_POST_FORCE
       | POST_RUN
       | END_OF_STEP;
}

/* ----------------------------------------------------------------------
 * make sure that group-sizes match at beginning of run
 * ---------------------------------------------------------------------- */

void FixONIOM::init()
{
  if (!strstr(update->integrate_style,"verlet"))
    error->all(FLERR,"ONIOM currently only supports run_style verlet");

  if(!comm->me && !master && modify->n_final_integrate){
    error->warning(FLERR,
                   "ONIOM: Integrating fixes present in slave partition, "
                   "will probably not work as intended.");
  }

  int me = comm->me;
  MPI_Request req;

  /* - make sure number of atoms match-up
   * - create mapping between partition-local tags and
   *   oniom-specific atom-number
   * - initial matching of positions
   */

  if(!master){
    // broadcast nat to master for checking of consistency
    auto& con = connections.front();
    const int bitmask = group->bitmask[con.mc_group];
    auto root = (comm->me == 0) ? MPI_ROOT : MPI_PROC_NULL;
    MPI_Bcast(&con.mc_nat, 1, MPI_INT, root, con.comm);

    // collect process-local tags
    decltype (con.tags) tags{};
    tags.reserve(con.mc_nat);
    for(int i=0; i<atom->nlocal; ++i){
      if(atom->mask[i] & bitmask)
        tags.push_back(atom->tag[i]);
    }

    // gather number of tags
    int send_count = tags.size() * sizeof(decltype(con.tags)::value_type);
    MPI_Allgather(&send_count, 1, MPI_INT, recv_count_buf, 1, MPI_INT, world);

    // construct displacement
    displs_buf[0] = 0;
    for(int i=1; i<universe->nprocs; ++i){
      displs_buf[i] = displs_buf[i-1]+recv_count_buf[i-1];
    }

    // gather and sort tags
    con.tags.clear();
    con.tags.resize(con.mc_nat);
    MPI_Allgatherv(tags.data(), send_count, MPI_BYTE,
                   con.tags.data(), recv_count_buf, displs_buf, MPI_BYTE, world);
    std::sort(con.tags.begin(), con.tags.end());

    // create hashtable
    con.hash.clear();
    for(int i=0; i<con.mc_nat; ++i){
      con.hash[con.tags[i]] = i;
    }

    // report tags if requested
    if((comm->me == 0) && (verbose > 1)){
      const char fmt[] = "tags[%d]=" TAGINT_FORMAT
          "  hash[" TAGINT_FORMAT "]=" TAGINT_FORMAT "\n";
      if (screen) fputs("slave tags\n", screen);
      if (logfile) fputs("slave tags\n", logfile);
      for(int i=0; i<con.mc_nat; ++i){
        if (screen) fprintf(screen, fmt, i, con.tags[i],
                            con.tags[i], con.hash[con.tags[i]]);
        if (logfile) fprintf(logfile, fmt, i, con.tags[i],
                             con.tags[i], con.hash[con.tags[i]]);

      }
    }

    // finally, after all is set up, do a first position synchronization
    receive_positions();
  }else{
    for(auto& con: connections){
      const int bitmask = group->bitmask[con.mc_group];
      // receive nat from slaves
      int nat;
      MPI_Bcast(&nat, 1, MPI_INT, 0, con.comm);

      if(con.mc_nat != nat)
        error->all(FLERR, "Inconsistent number of MC atoms");

      // collect process-local tags
      decltype(con.tags) tags{};
      tags.reserve(con.mc_nat);
      for(int i=0; i<atom->nlocal; ++i){
        if(atom->mask[i] & bitmask)
          tags.push_back(atom->tag[i]);
      }

      // gather number of tags
      int send_count = tags.size() * sizeof(decltype(con.tags)::value_type);
      MPI_Allgather(&send_count, 1, MPI_INT, recv_count_buf, 1, MPI_INT, world);

      // construct displacement
      displs_buf[0] = 0;
      for(int i=1; i<universe->nprocs; ++i){
        displs_buf[i] = displs_buf[i-1]+recv_count_buf[i-1];
      }

      // gather and sort tags
      con.tags.clear();
      con.tags.resize(con.mc_nat);
      MPI_Allgatherv(tags.data(), send_count, MPI_BYTE,
                     con.tags.data(), recv_count_buf, displs_buf, MPI_BYTE, world);
      std::sort(con.tags.begin(), con.tags.end());

      // create hashtable
      con.hash.clear();
      for(int i=0; i<con.mc_nat; ++i){
        con.hash[con.tags[i]] = i;
      }
      // report tags
      if((comm->me == 0) && (verbose > 1)){
        if (screen) fprintf(screen, "master tags\n");
        if (logfile) fprintf(logfile, "master tags\n");
        const char fmt[] = "tags[%d]=" TAGINT_FORMAT
            "  hash[" TAGINT_FORMAT "]=" TAGINT_FORMAT "\n";
        for(int i=0; i<con.mc_nat; ++i){
          if (screen) fprintf(screen, fmt, i, con.tags[i],
                              con.tags[i], con.hash[tags[i]]);
          if (logfile) fprintf(logfile, fmt, i, con.tags[i],
                               con.tags[i], con.hash[tags[i]]);

        }
      }
    }
    // finally, after all is set up, do a first position synchronization
    send_positions();
  }
}

/* ---------------------------------------------------------------------- */

void FixONIOM::send_positions()
{
  if(!master) return;
  const int nlocal = atom->nlocal;
  const int * const mask = atom->mask;
  double * const * const x = atom->x;
  const double * const q = atom->q;
  const tagint * const tag = atom->tag;

  if((comm->me == 0) && (verbose > 0)){
    const char fmt[] = "ONIOM: sending positions\n";
    if (screen) fprintf(screen, fmt);
    if (logfile) fprintf(logfile, fmt);
  }

  // send updated positions to slaves
  for (const auto& con: connections) {
    const int bitmask = group->bitmask[con.mc_group];

    // collect data in world
    comm_buf.clear();
    comm_buf.reserve(con.mc_nat);
    for (int i=0; i<nlocal; ++i) {
      if (mask[i] & bitmask) {
        comm_buf.push_back({
          con.hash.at(tag[i]),
          x[i][0], x[i][1], x[i][2],
//          q[i]
        });
      }
    }
    int send_count = comm_buf.size() * sizeof(decltype (comm_buf)::value_type);
    // gather number of local atoms
    MPI_Gather(&send_count, 1, MPI_INT, recv_count_buf, 1, MPI_INT, 0, world);
    // construct displacement
    displs_buf[0] = 0;
    for(int i=1; i<universe->nprocs; ++i){
      displs_buf[i] = displs_buf[i-1]+recv_count_buf[i-1];
    }
    // gather comm_buf
    comm_buf.resize(con.mc_nat);
    auto gather_root = (comm->me == 0) ? MPI_IN_PLACE : comm_buf.data();
    MPI_Gatherv(gather_root, send_count, MPI_BYTE,
                comm_buf.data(), recv_count_buf, displs_buf, MPI_BYTE, 0, world);

    // transmit to slave-partition
    send_count = con.mc_nat * sizeof(decltype (comm_buf)::value_type);
    auto root = (comm->me == 0) ? MPI_ROOT : MPI_PROC_NULL;
    MPI_Bcast(comm_buf.data(), send_count, MPI_BYTE, root, con.comm);
  }
}

void FixONIOM::receive_positions()
{
  if(master) return;

  const auto& con = connections.front();
  const int nlocal = atom->nlocal;
  const int * const mask = atom->mask;
  const int bitmask = group->bitmask[con.mc_group];
  double * const * const x = atom->x;
  double * const q = atom->q;
  const tagint * const tag = atom->tag;

  // status report if requested
  if((comm->me == 0) && (verbose > 0)){
    const char fmt[] = "ONIOM: receiving positions from master\n";
    if (screen) fprintf(screen, fmt);
    if (logfile) fprintf(logfile, fmt);
  }

  // receive updated positions from master
  comm_buf.clear();
  comm_buf.resize(con.mc_nat);
  int recv_count = con.mc_nat * sizeof(decltype (comm_buf)::value_type);
  MPI_Bcast(comm_buf.data(), recv_count, MPI_BYTE, 0, con.comm);

  // update positions of relevant atoms
  for(int i=0; i<nlocal; ++i){
    if(mask[i] & bitmask){
      for(const auto& dat: comm_buf){
        if(con.tags[dat.tag] == tag[i]){
          x[i][0] = dat.x;
          x[i][1] = dat.y;
          x[i][2] = dat.z;
          q[i] = dat.q;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixONIOM::receive_forces()
{
  if(!master) return;
  const int nlocal = atom->nlocal;
  const int * const mask = atom->mask;
  double * const *const f = atom->f;
  const tagint * const tag = atom->tag;

  if((comm->me == 0) && (verbose > 0)){
    const char fmt[] = "ONIOM: collecting forces from slave\n";
    if (screen) fprintf(screen, fmt);
    if (logfile) fprintf(logfile, fmt);
  }

  for(auto& con: connections){
    const int bitmask = group->bitmask[con.mc_group];
    // collect new forces from slaves
    comm_buf.clear();
    comm_buf.resize(con.mc_nat);
    int recv_count = con.mc_nat * sizeof(decltype(comm_buf)::value_type);
    MPI_Bcast(comm_buf.data(), recv_count, MPI_BYTE, 0, con.comm);

    // update forces
    if(con.mode & MINUS){
      for(int i=0; i<nlocal; ++i){
        if(mask[i] & bitmask){
          for(int j=0; j<con.mc_nat; ++j){
            if(con.tags[comm_buf[j].tag] == tag[i]){
              f[i][0] -= comm_buf[j].x;
              f[i][1] -= comm_buf[j].y;
              f[i][2] -= comm_buf[j].z;
            }
          }
        }
      }
    }else{
      for(int i=0; i<nlocal; ++i){
        if(mask[i] & bitmask){
          for(int j=0; j<con.mc_nat; ++j){
            if(con.tags[comm_buf[j].tag] == tag[i]){
              f[i][0] += comm_buf[j].x;
              f[i][1] += comm_buf[j].y;
              f[i][2] += comm_buf[j].z;
            }
          }
        }
      }
    }
  }
}

void FixONIOM::send_forces()
{
  if(master) return;

  const int nlocal = atom->nlocal;
  const int * const mask = atom->mask;
  double * const *const f = atom->f;
  const tagint * const tag = atom->tag;
  const auto& con = connections.front();
  const int bitmask = group->bitmask[con.mc_group];

  // status report if requested
  if((comm->me == 0) && (verbose > 0)){
    const char fmt[] = "ONIOM: sending forces to master\n";
    if (screen) fprintf(screen, fmt);
    if (logfile) fprintf(logfile, fmt);
  }

  // collect data in world
  comm_buf.clear();
  comm_buf.reserve(con.mc_nat);
  for(int i=0; i<nlocal; ++i){
    if(mask[i] & bitmask){
      comm_buf.push_back({
        con.hash.at(tag[i]),
        f[i][0], f[i][1], f[i][2],
        -1
      });
    }
  }
  int send_count = comm_buf.size() * sizeof(decltype(comm_buf)::value_type);
  int recv_count[universe->nprocs];
  // gather number of local atoms
  MPI_Gather(&send_count, 1, MPI_INT, recv_count, 1, MPI_INT, 0, world);
  // construct displacement
  int displs[universe->nprocs];
  displs[0] = 0;
  for(int i=1; i<universe->nprocs; ++i){
    displs[i] = displs[i-1]+recv_count[i-1];
  }
  // gather comm_buf
  comm_buf.resize(con.mc_nat);
  auto gather_root = (comm->me == 0) ? MPI_IN_PLACE : comm_buf.data();
  MPI_Gatherv(gather_root, send_count, MPI_BYTE,
              comm_buf.data(), recv_count, displs, MPI_BYTE, 0, world);

  //transmit to master
  auto root = (comm->me == 0) ? MPI_ROOT : MPI_PROC_NULL;
  MPI_Bcast(comm_buf.data(), send_count, MPI_BYTE, root, con.comm);
}

/* ---------------------------------------------------------------------- */

void FixONIOM::post_integrate()
{
  if(master){
    // MD master
    send_positions();
  }else{
    // MD & MIN slave
    receive_positions();
  }
}

/* ---------------------------------------------------------------------- */

void FixONIOM::setup(int)
{
  run_once = false;
  if(master){
    // MD master
    receive_forces();
  }else{
    // MD & MIN slave
    send_forces();
  }
}

void FixONIOM::min_setup(int)
{
  // MIN master
  receive_forces();
}

/* ---------------------------------------------------------------------- */

void FixONIOM::post_force(int)
{
  if(master){
    // MD master
    receive_forces();
  }else{
    // MD & MIN slave
    send_forces();
  }
}

/* ---------------------------------------------------------------------- */

void FixONIOM::min_pre_force(int)
{
  // MIN master
  // tell slaves to continue faux MD after first step
  if(!run_once){
      run_once = true;
  }else{
    if((comm->me == 0) && (verbose > 0)){
      const char msg[] = "ONIOM: Minimization not converged, continuing slave-MDs\n";
      if (screen) fprintf(screen, msg);
      if (logfile) fprintf(logfile, msg);
    }
    for(auto& con: connections){
      auto root = (comm->me == 0) ? MPI_ROOT : MPI_PROC_NULL;
      uint8_t f{false};
      MPI_Bcast(&f, 1, MPI_BYTE, root, con.comm);
    }
  }
  send_positions();
}

/* ---------------------------------------------------------------------- */

void FixONIOM::min_post_force(int)
{
  // MIN master
  receive_forces();
}

/* ---------------------------------------------------------------------- */

void FixONIOM::end_of_step()
{
  // MIN slave
  if(minimize && !master){
    // get convergence state
    uint8_t converged{0};
    MPI_Bcast(&converged, 1, MPI_BYTE, 0, connections.front().comm);
    if(converged){
      // trigger timeout if converged
      if((comm->me == 0) && (verbose > 0)){
        const char msg[] = "ONIOM: Ending run early because master converged minimization\n";
        if (screen) fprintf(screen, msg);
        if (logfile) fprintf(logfile, msg);
      }
      timer->force_timeout();
    }else if((comm->me == 0) && (verbose > 0)){
      const char msg[] = "ONIOM: Continuing MD to accompany minimization\n";
      if (screen) fprintf(screen, msg);
      if (logfile) fprintf(logfile, msg);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixONIOM::post_run()
{
  // MIN master
  if(minimize && master){
    // tell slaves to end faux MD
    if((comm->me == 0) && (verbose > 0)){
      const char msg[] = "ONIOM: Minimization converged, cancelling slave-MDs\n";
      if (screen) fprintf(screen, msg);
      if (logfile) fprintf(logfile, msg);
    }
    for(auto& con: connections){
      auto root = (comm->me == 0) ? MPI_ROOT : MPI_PROC_NULL;
      uint8_t t{true};
      MPI_Bcast(&t, 1, MPI_BYTE, root, con.comm);
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixONIOM::memory_usage(void)
{
  double bytes;

  bytes = sizeof(FixONIOM);
  bytes += connections.capacity() * sizeof(conn_t);
  bytes += comm_buf.capacity() * sizeof(commdata_t);

  return bytes;
}

