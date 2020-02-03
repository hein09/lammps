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
   Contributing author:  Sebastian Gsänger (FAU)
------------------------------------------------------------------------- */

#include <algorithm>
#include <unistd.h>
#include <cstring>

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
  FixCollect(lmp, narg, arg)
{
  master = universe->iworld == 0;

  if (narg < 3)
    error->all(FLERR,"Illegal fix oniom command");

  if (atom->tag_enable == 0)
    error->all(FLERR,"Oniom requires atom IDs");

  if (atom->tag_consecutive() == 0)
    error->all(FLERR,"Oniom requires consecutive atom IDs");

  // parse arguments
  int iarg=3;
  int iworld = universe->iworld;
  while(iarg < narg){
    if(strcmp(arg[iarg], "plus") == 0){
      /******************************************************
       * Add forces of target partition
       ******************************************************/
      if(iarg+2 > narg)
        error->universe_all(FLERR, "Illegal oniom command");

      // arg[iarg+1] is slave partition
      int target = force->inumeric(FLERR, arg[iarg+1])-1;
      if((target < 1 ) || (target>=universe->nworlds)){
        char msg[40];
        sprintf(msg, "Invalid partition %i in oniom command", target);
        error->universe_all(FLERR, msg);
      }

      // if we're not involved in this exchange, skip creating a connection
      if(!master && iworld != target){
          iarg += 3;
          continue;
      }

      // arg[iarg+2] is group, only evaluated on involved partitions
      auto coll = get_collection(arg[iarg+2]);
      if(group->dynamic[coll->group_id])
        error->all(FLERR,"Oniom group must be static");

      // create connections
      if(master){
        connections.push_back(conn_t{FixONIOM::MINUS, target, connect_to_slave(target), coll});
      }else{
        if(!connections.empty()){
          error->universe_one(FLERR, "Same partition used in multiple oniom regions");
        }
        connections.push_back(conn_t{0, 0, connect_to_master(), coll});
      }
      iarg += 3;
    }else if(strcmp(arg[iarg], "minus") == 0){
      /******************************************************
       * Subtract forces of target partition
       ******************************************************/
      if(iarg+2 > narg)
        error->universe_all(FLERR, "Illegal oniom command");

      // arg[iarg+1] is slave partition
      int target = force->inumeric(FLERR, arg[iarg+1])-1;
      if((target < 1 ) || (target>=universe->nworlds)){
        char msg[40];
        sprintf(msg, "Invalid partition %i in oniom command", target);
        error->universe_all(FLERR, msg);
      }

      // if we're not involved in this exchange, skip creating a connection
      if(!master && iworld != target){
          iarg += 3;
          continue;
      }

      // arg[iarg+2] is group, only evaluated on involved partitions
      auto coll = get_collection(arg[iarg+2]);
      if(group->dynamic[coll->group_id])
        error->all(FLERR,"Oniom group must be static");

      // create connections
      if(master){
        connections.push_back(conn_t{FixONIOM::PLUS, target, connect_to_slave(target), coll});
      }else{
        if(!connections.empty()){
          error->universe_one(FLERR, "Same partition used in multiple oniom regions");
        }
        connections.push_back(conn_t{0, 0, connect_to_master(), coll});
      }
      iarg += 3;
    }else if(strcmp(arg[iarg], "group") == 0){
      /********************************************************
       * Classic oniom setup
       *
       * `group low high g` is equivalent to
       * `minus low g plus high g`
       ********************************************************/
      if(iarg+3 > narg)
        error->universe_all(FLERR, "Illegal oniom command");

      // arg[iarg+1] is low-level slave partition
      int low = force->inumeric(FLERR, arg[iarg+1])-1;
      if((low < 1 ) || (low>=universe->nworlds)){
        char msg[40];
        sprintf(msg, "Invalid partition %i in oniom command", low);
        error->universe_all(FLERR, msg);
      }

      // arg[iarg+2] is high-level slave partition
      int high = force->inumeric(FLERR, arg[iarg+2])-1;
      if((high < 1 ) || (high>=universe->nworlds)){
        char msg[40];
        sprintf(msg, "Invalid partition %i in oniom command", high);
        error->universe_all(FLERR, msg);
      }

      // if we're not involved in this exchange, skip creating a connection
      if(!master && iworld != low && iworld != high){
          iarg += 4;
          continue;
      }

      // arg[iarg+3] is group, only evaluated on involved partitions
      auto coll = get_collection(arg[iarg+3]);
      if(group->dynamic[coll->group_id])
        error->all(FLERR,"Oniom group must be static");

      // create connections
      if(master){
        connections.push_back(conn_t{FixONIOM::MINUS, low, connect_to_slave(low), coll});
        connections.push_back(conn_t{FixONIOM::PLUS, high, connect_to_slave(high), coll});
      }else{
        if(!connections.empty()){
          error->universe_one(FLERR, "Same partition used in multiple oniom regions");
        }
        connections.push_back(conn_t{0, 0, connect_to_master(), coll});
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
 * connect partitions
 *
 * this may fail spectacularly (unexpected input, strangely invalid MPI state)
 * trying to fail gracefully via nonblocking-comm and timing out
 * ---------------------------------------------------------------------- */

MPI_Comm FixONIOM::connect_to_slave(int target)
{
  if(!master){
    char msg[40];
    sprintf(msg, "fix oniom: trying to connect non-master"
                 " partition %i to partition %i", universe->iworld, target);
    error->all(FLERR, msg);
  }
  int flag{0}, tag;
  MPI_Comm out_comm;
  MPI_Request req;
  if(comm->me==0){
    MPI_Isend(&tag, 1, MPI_INT,
              universe->root_proc[target],
              target, universe->uworld, &req);
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
                       universe->root_proc[target],
                       target, &out_comm);
  return out_comm;
}

MPI_Comm FixONIOM::connect_to_master()
{
  if(master){
    error->all(FLERR, "fix oniom: trying to connect partition 0 to itself");
  }
  int flag{0}, tag;
  MPI_Comm out_comm;
  MPI_Request req;
  if(comm->me==0){
    MPI_Irecv(&tag, 1, MPI_INT,
              universe->root_proc[0],
              universe->iworld, universe->uworld, &req);
    // wait for one minute to establish connection
    int j{0};
    for(; j<60; ++j){
      MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
      if(flag) break;
      usleep(1000000);
    }
    if(j>=60) error->universe_one(FLERR, "Could not establish connection to oniom master");
  }
  MPI_Intercomm_create(world, 0, universe->uworld,
                       universe->root_proc[0],
                       universe->iworld,
                       &out_comm);
  return out_comm;
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

  if(!master){
    // broadcast nat to master for checking of consistency
    const auto& con = connections.front();
    const auto& coll = *con.coll;
    const auto root = (comm->me == 0) ? MPI_ROOT : MPI_PROC_NULL;
    MPI_Bcast(&con.coll->nat, 1, MPI_INT, root, con.comm);

    // report tags if requested
    if((comm->me == 0) && (verbose > 1)){
      if (screen) fputs("slave tags\n", screen);
      if (logfile) fputs("slave tags\n", logfile);
      const char fmt[] = "LAMMPS-tag[" TAGINT_FORMAT "] == "
                         "ONIOM-idx[%d]";
      for(int i=0; i<coll.nat; ++i){
        if (screen) fprintf(screen, fmt, coll.idx2tag[i],
                            coll.tag2idx.at(coll.idx2tag[i]));
        if (logfile) fprintf(logfile, fmt, coll.idx2tag[i],
                             coll.tag2idx.at(coll.idx2tag[i]));

      }
    }

    // {fter all is set up, do a first position synchronization
    receive_positions();
  }else{
    for(const auto& con: connections){
      const auto& coll = *con.coll;
      // receive nat from slaves
      int nat;
      MPI_Bcast(&nat, 1, MPI_INT, 0, con.comm);

      if(coll.nat != nat)
        error->all(FLERR, "Inconsistent number of MC atoms");

      // report tags if requested
      if((comm->me == 0) && (verbose > 1)){
        if (screen) fprintf(screen, "master tags (to %d)\n", con.target);
        if (logfile) fprintf(logfile, "master tags (to %d)\n", con.target);
        const char fmt[] = "LAMMPS-tag[" TAGINT_FORMAT "] == "
                           "ONIOM-idx[%d]";
        for(int i=0; i<coll.nat; ++i){
          if (screen) fprintf(screen, fmt, coll.idx2tag[i],
                              coll.tag2idx.at(coll.idx2tag[i]));
          if (logfile) fprintf(logfile, fmt, coll.idx2tag[i],
                               coll.tag2idx.at(coll.idx2tag[i]));
        }
      }
    }

    // after all is set up, do a first position synchronization
    send_positions();
  }
}

/* ---------------------------------------------------------------------- */

void FixONIOM::send_positions()
{
  if(!master) return;

  for (const auto& con: connections) {
    auto& coll = *con.coll;
    auto& buffer = coll.buffer;

    // status report if requested
    if((comm->me == 0) && (verbose > 0)){
      const char fmt[] = "ONIOM: sending positions to slave %d\n";
      if (screen) fprintf(screen, fmt, con.target);
      if (logfile) fprintf(logfile, fmt, con.target);
    }

    // collect updated positions in buffer on proc 0
    collect_lammps(coll, atom->x);

    // transmit to slave-partition
    MPI_Bcast(buffer.data(), coll.nat * sizeof(commdata_t), MPI_BYTE,
              (comm->me == 0) ? MPI_ROOT : MPI_PROC_NULL, con.comm);
  }
}

void FixONIOM::receive_positions()
{
  if(master) return;

  const auto& con = connections.front();
  auto& coll = *con.coll;
  auto& buffer = coll.buffer;

  // status report if requested
  if((comm->me == 0) && (verbose > 0)){
    const char fmt[] = "ONIOM: receiving positions from master\n";
    if (screen) fprintf(screen, fmt);
    if (logfile) fprintf(logfile, fmt);
  }

  // receive updated positions from master
  buffer.resize(coll.nat);
  MPI_Bcast(buffer.data(), coll.nat * sizeof(commdata_t),
            MPI_BYTE, 0, con.comm);

  // update positions of relevant atoms
  for(int i=0; i<atom->nlocal; ++i){
    for(const auto& dat: buffer){
      if(atom->tag[i] == dat.tag){
        atom->x[i][0] = dat.x[0];
        atom->x[i][1] = dat.x[1];
        atom->x[i][2] = dat.x[2];
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
    // receive forces from slaves
    auto& coll = *con.coll;
    auto& buffer = coll.buffer;
    buffer.resize(coll.nat);
    MPI_Bcast(buffer.data(), coll.nat * sizeof(commdata_t),
              MPI_BYTE, 0, con.comm);

    // update forces of relevant atoms
    double fscale = (con.mode & MINUS) ? -1 : 1;
    for(int i=0; i<atom->nlocal; ++i){
      for(const auto& dat: buffer){
        if(atom->tag[i] == dat.tag){
          atom->f[i][0] += dat.x[0] * fscale;
          atom->f[i][1] += dat.x[1] * fscale;
          atom->f[i][2] += dat.x[2] * fscale;
        }
      }
    }
  }
}

void FixONIOM::send_forces()
{
  if(master) return;

  const auto& con = connections.front();
  auto& coll = *con.coll;
  auto& buffer = coll.buffer;

  // status report if requested
  if((comm->me == 0) && (verbose > 0)){
    const char fmt[] = "ONIOM: sending forces to master\n";
    if (screen) fprintf(screen, fmt);
    if (logfile) fprintf(logfile, fmt);
  }

  collect_lammps(coll, atom->f);

  //transmit to master
  MPI_Bcast(buffer.data(), coll.nat * sizeof(commdata_t), MPI_BYTE,
            (comm->me == 0) ? MPI_ROOT : MPI_PROC_NULL, con.comm);
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
  ran_once = false;
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
  if(!ran_once){
      ran_once = true;
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
  bytes += FixCollect::memory_usage();

  return bytes;
}

