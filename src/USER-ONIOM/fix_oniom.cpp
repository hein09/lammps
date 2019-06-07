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
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include "fix_oniom.h"
#include "atom.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "universe.h"
#include "modify.h"

// message tags for ONIOM inter communicator communication
enum {ONIOM_TAG_OTHER=0, ONIOM_TAG_SIZE=1, ONIOM_TAG_COORD=2,
      ONIOM_TAG_FORCE=3, ONIOM_TAG_FORCE2=4, ONIOM_TAG_CELL=5,
      ONIOM_TAG_RADII=6, ONIOM_TAG_CHARGE=7, ONIOM_TAG_TYPE=8,
      ONIOM_TAG_MASS=9};

using namespace LAMMPS_NS;
using namespace FixConst;

/***************************************************************
 * create class and parse arguments in LAMMPS script. Syntax:
 * fix ID group-ID oniom (master[_ec] <low> <high> | slave[_ec] <master>) mc_group [ec_group]
 *
 * sets up connections between partitions
 ***************************************************************/
FixONIOM::FixONIOM(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{

  if (narg < 5)
    error->all(FLERR,"Illegal fix oniom command");

  if (strcmp(update->unit_style,"metal") == 0) {
    fscale = 1.0/23.0609;
  } else if (strcmp(update->unit_style,"real") == 0) {
    fscale = 1.0;
  } else error->all(FLERR,"Fix oniom requires real or metal units");

  if (atom->tag_enable == 0)
    error->all(FLERR,"Fix oniom requires atom IDs");

  if (!atom->q_flag)
    error->all(FLERR,"Fix oniom requires atom attribute q");

  if (atom->tag_consecutive() == 0)
    error->all(FLERR,"Fix oniom requires consecutive atom IDs");

  int iarg = 3;
  while(iarg < narg) {
    if (strstr(arg[iarg], "master")){
      // check for electrostatic coupling
      int ec_flag = (bool)strstr(arg[iarg], "_ec");
      if(ec_flag){
          error->all(FLERR, "EC not implemented yet");
        if (iarg+4 >= narg) error->all(FLERR, "Illegal fix oniom master_ec command");
      }else{
        if (iarg+3 >= narg) error->all(FLERR, "Illegal fix oniom master command");
      }

      // arg[iarg+1] is low-level slave partition
      int low = force->inumeric(FLERR, arg[iarg+1])-1;
      if((low < 0 ) || (low == universe->iworld) || (low>=universe->nworlds)){
        error->all(FLERR, "Invalid partition in fix oniom master command");
      }

      // arg[iarg+2] is high-level slave partition
      int high = force->inumeric(FLERR, arg[iarg+2])-1;
      if((high < 0 ) || (high == universe->iworld) || (high>=universe->nworlds)){
        error->all(FLERR, "Invalid partition in fix oniom master command");
      }

      // arg[iarg+3] is mc-group
      int mc_group = group->find(arg[iarg+3]);
      if(mc_group == -1)
        error->all(FLERR,"Could not find fix oniom master group ID");
      if(group->dynamic[mc_group])
        error->all(FLERR,"Fix oniom MC group must be static");

      // arg[iarg+4] is ec-group (if needed)
      int ec_group{-1};
      if(ec_flag){
        int ec_group = group->find(arg[iarg+4]);
        if(ec_group == -1)
          error->all(FLERR,"Could not find fix oniom master_ec group ID");
        if(ec_group == mc_group)
          error->all(FLERR, "Electrostatic and mechanical group IDs must differ");
      }

      // create (hopefully) unique tag from partition-id + master's mc-group name
      int tag{universe->iworld<<2};
      for(size_t i=0; i<strlen(arg[iarg+3]); ++i){
        tag += arg[iarg+3][i];
      }

      // establish connection to other partitions, if succesful create intercommunicator
      auto setupMPI = [&](int part, int& ec, MPI_Comm& comm){
        int me, flag{0}, child_ec{0};
        MPI_Request req;
        MPI_Comm_rank(world, &me);
        if(me==0){
          MPI_Isend(&tag, 1, MPI_INT,
                    universe->root_proc[part],
                    part, universe->uworld, &req);
          // wait for one minute to establish connection
          int j;
          for(j=0; j<60; ++j){
              MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
              if(flag) break;
              usleep(1000000);
          }
          if(j >= 60) error->one(FLERR, "Could not establish connection to oniom slave");
          MPI_Recv(&child_ec, 1, MPI_INT, universe->root_proc[part],
                   tag, universe->uworld, MPI_STATUS_IGNORE);
        }
        MPI_Bcast(&child_ec, 1, MPI_INT, 0, world);
        if(!ec_flag && child_ec){
            error->all(FLERR, "Fix oniom master cannot serve fix oniom slave_ec");
        }
        ec = child_ec ? ec_group : -1;
        MPI_Intercomm_create(world, 0, universe->uworld,
                             universe->root_proc[part], tag, &comm);
      };

      // setup connections
      conn_down_t conn{};
      conn.mc_group = mc_group;
      conn.part_low = low;
      conn.part_high = high;
      setupMPI(conn.part_low, conn.ec_low, conn.comm_low);
      setupMPI(conn.part_high, conn.ec_high, conn.comm_high);
      // save connections
      downwards.insert({arg[iarg+3], conn});
      if(ec_flag){
        iarg += 5;
      }else{
        iarg += 4;
      }

    }else if (strstr(arg[iarg], "slave")){
      if(upwards)
        error->all(FLERR,"Illegal fix oniom slave command");

      // check for electrostatic_coupling
      int ec_flag = (bool)strstr(arg[iarg], "_ec");
      if(ec_flag){
          error->all(FLERR, "EC not implemented yet");
        if (iarg+3 >= narg) error->all(FLERR, "Illegal fix oniom slave_ec command");
      }else{
        if (iarg+2 >= narg) error->all(FLERR, "Illegal fix oniom slave command");
      }

      // arg[iarg+1] is master partition
      int master = force->inumeric(FLERR, arg[iarg+1])-1;
      if((master < 0 ) || (master == universe->iworld) || (master>=universe->nworlds)){
        error->all(FLERR, "Invalid partition in fix oniom slave command");
      }

      // arg[iarg+2] is mc-group
      int mc_group = group->find(arg[iarg+2]);
      if(mc_group == -1)
        error->all(FLERR,"Could not find fix oniom slave group ID");
      if(group->dynamic[mc_group])
        error->all(FLERR,"Fix oniom MC group must be static");

      // arg[iarg+3] is ec-group (if needed)
      int ec_group{-1};
      if(ec_flag){
        ec_group = group->find(arg[iarg+3]);
        if(ec_group == -1)
          error->all(FLERR,"Could not find fix oniom slave_ec group ID");
      }

      // establish connection to other partition, if succesful create intercommunicator
      int me, flag{0}, tag;
      MPI_Request req;
      MPI_Comm_rank(world, &me);
      if(me==0){
        MPI_Irecv(&tag, 1, MPI_INT,
                  universe->root_proc[master],
                  universe->iworld, universe->uworld, &req);
        // wait for one minute to establish connection
        int j;
        for(j=0; j<60; ++j){
          MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
          if(flag) break;
          usleep(1000000);
        }
        if(j>=60) error->one(FLERR, "Could not establish connection to oniom master");
      }
      MPI_Bcast(&tag, 1, MPI_INT, 0, world);
      // request mc/ec
      MPI_Send(&ec_flag, 1, MPI_INT,
               universe->root_proc[master],
               tag, universe->uworld);
      MPI_Comm newcomm;
      MPI_Intercomm_create(world, 0, universe->uworld,
                           universe->root_proc[master], tag, &newcomm);
      // save connection
      upwards = new conn_up_t{master, mc_group, -1, ec_group, -1, newcomm};
      if(ec_flag){
        iarg += 4;
      }else{
        iarg += 3;
      }

    }else if (!strcmp(arg[iarg], "verbose")){
      verbose = force->inumeric(FLERR, arg[iarg+1]);
      iarg += 2;

    }else{
      error->all(FLERR, "Illegal fix oniom command");
    }
  }

}

/*********************************
 * Clean up on deleting the fix. *
 *********************************/
FixONIOM::~FixONIOM()
{
    // TODO: free communicators?
    if(upwards) delete upwards;
}

/* ---------------------------------------------------------------------- */
int FixONIOM::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ----------------------------------------------------------------------
 * make sure that group-sizes match at beginning of run
 * ---------------------------------------------------------------------- */

void FixONIOM::init()
{
  if (strstr(update->integrate_style,"respa"))
    error->all(FLERR,"Fix oniom does not currently support r-RESPA");

  bool before{true}, flag{false};
  for(int i=0; i<modify->nfix; ++i){
    if(strcmp(id, modify->fix[i]->id) == 0) before = false;
    else if (before && (modify->fmask[i] & FINAL_INTEGRATE)) flag = true;
  }
  if(flag){
    error->all(FLERR, "Fix oniom should come before other integrating fixes.");
  }

  if(!comm->me && upwards && modify->n_final_integrate){
    error->warning(FLERR,
      "Using fix oniom slave with other integrating fixes.");
  }

  int me = comm->me;
  MPI_Request req;

  /* exchanging numbers of atoms between master and slaves
   * and creating mapping between partition-local tags and
   * oniom-specific atom-number
   */

  if(upwards){
    // send nat_mc
    int nat[2];
    bigint tmp_mc = group->count(upwards->mc_group);
    if(tmp_mc > MAXSMALLINT){
      error->all(FLERR,"Too many MC atoms for fix oniom");
    }
    nat[0] = static_cast<int>(tmp_mc);

    // send nat_ec
    if(upwards->ec_group != -1){
      bigint tmp_ec = group->count(upwards->ec_group);
      if(tmp_ec > MAXSMALLINT){
        error->all(FLERR,"Too many EC atoms for fix oniom");
      }
      nat[1] = static_cast<int>(tmp_ec);
    }else{
      nat[1] = -1;
    }

    // broadcast nat to master for checking of consistency
    MPI_Isend(nat, 2, MPI_INT, 0, ONIOM_TAG_SIZE, upwards->comm, &req);

    // save nat in whole partition
    upwards->mc_nat = nat[0];
    upwards->ec_nat = nat[1];

    // setup hashtable
    decltype (upwards->tags) tags{};
    tags.reserve(nat[0]);
    for(int i=0; i<atom->nlocal; ++i){
      if(atom->mask[i] & group->bitmask[upwards->mc_group])
        tags.push_back(atom->tag[i]);
    }
    int send_count = tags.size() * sizeof(decltype(upwards->tags)::value_type);
    int recv_count[universe->nprocs];
    // gather number of tags
    MPI_Allgather(&send_count, 1, MPI_INT, recv_count, 1, MPI_INT, world);
    // construct displacement
    int displs[universe->nprocs];
    displs[0] = 0;
    for(int i=1; i<universe->nprocs; ++i){
      displs[i] = displs[i-1]+recv_count[i-1];
    }
    // gather and sort tags
    upwards->tags.clear();
    upwards->tags.resize(nat[0]);
    MPI_Allgatherv(tags.data(), send_count, MPI_BYTE,
                   upwards->tags.data(), recv_count, displs, MPI_BYTE, world);
    std::sort(upwards->tags.begin(), upwards->tags.end());
    // create hashtable
    upwards->hash.clear();
    for(int i=0; i<nat[0]; ++i){
      upwards->hash[upwards->tags[i]] = i;
    }
    // report tags
    if((comm->me == 0) && (verbose > 1)){
      const char fmt[] = "tags[%d]=" TAGINT_FORMAT
          "  hash[" TAGINT_FORMAT "]=" TAGINT_FORMAT "\n";
      if (screen) fputs("slave tags\n", screen);
      if (logfile) fputs("slave tags\n", logfile);
      for(int i=0; i<nat[0]; ++i){
        if (screen) fprintf(screen, fmt, i, upwards->tags[i],
                            upwards->tags[i], upwards->hash[upwards->tags[i]]);
        if (logfile) fprintf(logfile, fmt, i, upwards->tags[i],
                             upwards->tags[i], upwards->hash[upwards->tags[i]]);

      }
    }
  }

  for(auto& pair: downwards){
    auto& con = pair.second;

    // get nat_mc
    bigint tmp_mc = group->count(con.mc_group);
    if(tmp_mc > MAXSMALLINT){
      error->all(FLERR,"Too many MC atoms for fix oniom");
    }
    int nat_mc = static_cast<int>(tmp_mc);

    // get nat_ec
    int nat_ec = -1;
    auto ec_group = con.ec_low != -1 ? con.ec_low :
                    con.ec_high != -1 ? con.ec_high : -1 ;
    if(ec_group != -1){
      bigint tmp_ec = group->count(ec_group);
      if(tmp_ec > MAXSMALLINT){
        error->all(FLERR,"Too many EC atoms for fix oniom");
      }
      nat_ec = static_cast<int>(tmp_ec);
    }

    // receive nat from slaves
    int nat_low[2], nat_high[2];
    if(me == 0){
      MPI_Request req[2];
      MPI_Irecv(nat_low, 2, MPI_INT, 0, ONIOM_TAG_SIZE,
                con.comm_low, &req[0]);
      MPI_Irecv(nat_high, 2, MPI_INT, 0, ONIOM_TAG_SIZE,
                con.comm_high, &req[1]);
      MPI_Waitall(2, req, MPI_STATUS_IGNORE);
    }
    MPI_Bcast(&nat_low, 2, MPI_INT, 0, world);
    MPI_Bcast(&nat_high, 2, MPI_INT, 0, world);

    if((nat_mc != nat_low[0]) || (nat_mc != nat_high[0]))
      error->all(FLERR, "Inconsistent number of MC atoms");
    con.mc_nat = nat_mc;

    if(nat_ec != -1){
      if(((nat_low[1] != -1) && (nat_low[1] != nat_ec)) ||
         ((nat_high[1] != -1) && (nat_high[1] != nat_ec)))
        error->all(FLERR, "Inconsistent number of EC atoms");
      con.ec_nat = nat_ec;
    }

    // setup hashtable
    decltype(con.tags) tags{};
    tags.reserve(nat_mc);
    for(int i=0; i<atom->nlocal; ++i){
      if(atom->mask[i] & group->bitmask[con.mc_group])
        tags.push_back(atom->tag[i]);
    }
    int send_count = tags.size() * sizeof(decltype(con.tags)::value_type);
    int recv_count[universe->nprocs];
    // gather number of tags
    MPI_Allgather(&send_count, 1, MPI_INT, recv_count, 1, MPI_INT, world);
    // construct displacement
    int displs[universe->nprocs];
    displs[0] = 0;
    for(int i=1; i<universe->nprocs; ++i){
      displs[i] = displs[i-1]+recv_count[i-1];
    }
    // gather and sort tags
    con.tags.clear();
    con.tags.resize(nat_mc);
    MPI_Allgatherv(tags.data(), send_count, MPI_BYTE,
                   con.tags.data(), recv_count, displs, MPI_BYTE, world);
    std::sort(con.tags.begin(), con.tags.end());
    // create hashtable
    con.hash.clear();
    for(int i=0; i<nat_mc; ++i){
      con.hash[con.tags[i]] = i;
    }
    // report tags
    if((comm->me == 0) && (verbose > 1)){
      const char fmt[] = "tags[%d]=" TAGINT_FORMAT
          "  hash[" TAGINT_FORMAT "]=" TAGINT_FORMAT "\n";
      const char fmt2[] = "master group %s\n";
      if (screen) fprintf(screen, fmt2, pair.first.c_str());
      if (logfile) fprintf(logfile, fmt2, pair.first.c_str());
      for(int i=0; i<nat_mc; ++i){
        if (screen) fprintf(screen, fmt, i, con.tags[i],
                            con.tags[i], con.hash[tags[i]]);
        if (logfile) fprintf(logfile, fmt, i, con.tags[i],
                             con.tags[i], con.hash[tags[i]]);

      }
    }
  }

  //make sure that master received our message
  if(upwards) MPI_Wait(&req, MPI_STATUS_IGNORE);

  // finally, after all is set up, do a first position synchronization
  exchange_positions();
}

/* ---------------------------------------------------------------------- */

void FixONIOM::exchange_positions()
{
  //TODO: EC atoms
  const int nlocal = atom->nlocal;
  const int * const mask = atom->mask;
  double * const * const x = atom->x;
  const double * const q = atom->q;
  const tagint * const tag = atom->tag;

  // receive updated positions from master
  if (upwards) {
    if((comm->me == 0) && (verbose > 0)){
      const char fmt[] = "ONIOM: receiving positions from partition %d\n";
      if (screen) fprintf(screen, fmt, upwards->partition);
      if (logfile) fprintf(logfile, fmt, upwards->partition);
    }
    comm_buf.clear();
    comm_buf.resize(upwards->mc_nat);
    int recv_count = upwards->mc_nat * sizeof(decltype (comm_buf)::value_type);
    MPI_Request req;
    MPI_Ibcast(comm_buf.data(), recv_count, MPI_BYTE, 0, upwards->comm, &req);
    MPI_Wait(&req, MPI_STATUS_IGNORE);

    // update positions of relevant atoms
    for(int i=0; i<nlocal; ++i){
      if(mask[i] & group->bitmask[upwards->mc_group]){
        for(const auto& dat: comm_buf){
          if(upwards->tags[dat.tag] == tag[i]){
            x[i][0] = dat.x;
            x[i][1] = dat.y;
            x[i][2] = dat.z;
          }
        }
      }
    }
  }

  // send updated positions to slaves
  for (const auto& pair: downwards) {
    const auto& con = pair.second;
    if((comm->me == 0) && (verbose > 0)){
      const char fmt[] = "ONIOM: sending positions to partitions %d and %d\n";
      if (screen) fprintf(screen, fmt, con.part_low, con.part_high);
      if (logfile) fprintf(logfile, fmt, con.part_low, con.part_high);
    }

    // collect data in world
    comm_buf.clear();
    comm_buf.reserve(con.mc_nat);
    for (int i=0; i<nlocal; ++i) {
      if (mask[i] & group->bitmask[con.mc_group]) {
        comm_buf.push_back({
          con.hash.at(tag[i]),
          x[i][0], x[i][1], x[i][2],
          q[i]
        });
      }
    }
    int send_count = comm_buf.size() * sizeof(decltype (comm_buf)::value_type);
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

    // transmit to slave-partitions
    send_count = con.mc_nat * sizeof(decltype (comm_buf)::value_type);
    MPI_Request req[2];
    auto root = (comm->me == 0) ? MPI_ROOT : MPI_PROC_NULL;
    MPI_Ibcast(comm_buf.data(), send_count, MPI_BYTE, root, con.comm_low, &req[0]);
    MPI_Ibcast(comm_buf.data(), send_count, MPI_BYTE, root, con.comm_high, &req[1]);
    MPI_Waitall(2, req, MPI_STATUS_IGNORE);
  }
}

/* ---------------------------------------------------------------------- */

void FixONIOM::exchange_forces()
{
  //TODO: EC atoms
  const int nlocal = atom->nlocal;
  const int * const mask = atom->mask;
  double * const *const f = atom->f;
  const double * const q = atom->q;
  const tagint * const tag = atom->tag;

  if(upwards){
    // transmitting new forces to master
    if((comm->me == 0) && (verbose > 0)){
      const char fmt[] = "ONIOM: sending forces to partition %d\n";
      if (screen) fprintf(screen, fmt, upwards->partition);
      if (logfile) fprintf(logfile, fmt, upwards->partition);
    }

    // collect data in world
    comm_buf.clear();
    comm_buf.reserve(upwards->mc_nat);
    for(int i=0; i<nlocal; ++i){
      if(mask[i] & group->bitmask[upwards->mc_group]){
        comm_buf.push_back({
          upwards->hash.at(tag[i]),
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
    comm_buf.resize(upwards->mc_nat);
    auto gather_root = (comm->me == 0) ? MPI_IN_PLACE : comm_buf.data();
    MPI_Gatherv(gather_root, send_count, MPI_BYTE,
                comm_buf.data(), recv_count, displs, MPI_BYTE, 0, world);

    //transmit to master
    MPI_Request req;
    auto root = (comm->me == 0) ? MPI_ROOT : MPI_PROC_NULL;
    MPI_Ibcast(comm_buf.data(), send_count, MPI_BYTE, root, upwards->comm, &req);
    MPI_Wait(&req, MPI_STATUS_IGNORE);
  }

  for(auto& pair: downwards){
    // collect new forces from slaves
    const auto& con= pair.second;
    if((comm->me == 0) && (verbose > 0)){
      const char fmt[] = "ONIOM: collecting forces from partitions %d and %d\n";
      if (screen) fprintf(screen, fmt, con.part_low, con.part_high);
      if (logfile) fprintf(logfile, fmt, con.part_low, con.part_high);
    }
    comm_buf.clear();
    comm_buf.resize(2*con.mc_nat);
    int recv_count = con.mc_nat * sizeof(decltype(comm_buf)::value_type);
    MPI_Request req[2];
    MPI_Ibcast(comm_buf.data(), recv_count, MPI_BYTE, 0, con.comm_low, &req[0]);
    MPI_Ibcast(comm_buf.data()+con.mc_nat, recv_count, MPI_BYTE, 0, con.comm_high, &req[1]);
    MPI_Waitall(2, req, MPI_STATUS_IGNORE);

    /* update forces
     *
     * normal oniom-case: master-forces - low-forces + high-forces
     */
    for(int i=0; i<nlocal; ++i){
      if(mask[i] & group->bitmask[con.mc_group]){
        for(int j=0; j<con.mc_nat; ++j){
          if(con.tags[comm_buf[j].tag] == tag[i]){
            f[i][0] -= fscale * comm_buf[j].x;
            f[i][1] -= fscale * comm_buf[j].y;
            f[i][2] -= fscale * comm_buf[j].z;
          }
          int k=j+con.mc_nat;
          if(con.tags[comm_buf[k].tag] == tag[i]){
            f[i][0] += fscale * comm_buf[k].x;
            f[i][1] += fscale * comm_buf[k].y;
            f[i][2] += fscale * comm_buf[k].z;
          }
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixONIOM::post_integrate()
{
  exchange_positions();
}

/* ---------------------------------------------------------------------- */

void FixONIOM::setup(int)
{
  exchange_forces();
}

/* ---------------------------------------------------------------------- */

void FixONIOM::final_integrate()
{
  exchange_forces();
}

/* ---------------------------------------------------------------------- */
/* local memory usage. approximately. */
double FixONIOM::memory_usage(void)
{
  double bytes;

  bytes = sizeof(FixONIOM);
  if(upwards) bytes += sizeof(conn_up_t);
  bytes += downwards.size() * sizeof(decltype(downwards)::value_type);

  return bytes;
}

// TODO: make sure ec and mc groups don't overlap (ec dynamic?)
