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

#include "oniom.h"
#include "fix_oniom.h"

#include <unistd.h>

#include "run.h"

#include "error.h"
#include "universe.h"
#include "modify.h"
#include "input.h"
#include "force.h"
#include "group.h"

using namespace LAMMPS_NS;

/****************************************
 * usage:
 *
 * use one of:
 * > oniom run <runargs> <slave definition>
 * > oniom minimize <runargs> <slave definition>
 *
 * with <slave definition being one or multiple of:
 * > mc <pair of partitions> <mc_group>
 * > ec <pair of partitions> <mc_group> <ec_group>
 ****************************************/
void Oniom::command(int narg, char **arg)
{
  if(strcmp(arg[0], "run")==0){
    run(narg-1, &arg[1]);
  }else if(strcmp(arg[1], "minimize")==0){
    // minimization not implemented yet
    error->universe_all(FLERR, "Illegal oniom command");
    minimize(narg-1, &arg[1]);
  }else{
    error->universe_all(FLERR, "Illegal oniom command");
  }
}

void Oniom::run(int narg, char **arg)
{
  // extract run call
  int nrunarg = -1;
  char **runarg = arg;
  for(int i=0; i<narg; ++i){
    if((strcmp(arg[i], "mc")==0) ||
       (strcmp(arg[i], "ec")==0) ||
       (strcmp(arg[i], "verbose")==0)){
      nrunarg = i;
      break;
    }
  }
  if(nrunarg == -1){
    error->universe_all(FLERR, "Illegal oniom run command");
  }
  // create oniom fix
  mkFix(narg-nrunarg, &arg[nrunarg]);

  // perform calculation
  Run run{lmp};
  run.command(nrunarg, arg);

  // delete oniom fix
  modify->delete_fix("ONIOM");
}

void Oniom::minimize(int narg, char **arg)
{
  // extract minimize call
  int nminarg = 5;
  // create oniom fix
  mkFix(narg-5, &arg[5]);

  // delete oniom fix
  modify->delete_fix("ONIOM");
}

FixONIOM* Oniom::mkFix(int narg, char **arg)
{
  // create private oniom fix
  char **fixarg = new char*[3];
  fixarg[0] = (char *) "ONIOM";
  fixarg[1] = (char *) "all";
  fixarg[2] = (char *) "ONIOM";
  modify->add_fix(3, fixarg);
  delete[] fixarg;
  auto* fix = (FixONIOM*) modify->fix[modify->nfix-1];

  // parse arguments
  int iarg=0;
  int iworld = universe->iworld;
  while(iarg < narg){

    /********************************************************
     * Mechanical or electrostatic coupling
     *
     * Needs two slave partitions and a static group of atoms
     * Additional group needed for electrostatic coupling,
     * should not overlap with first, but may be dynamic
     ********************************************************/
    if((strcmp(arg[iarg], "mc") == 0) || (strcmp(arg[iarg], "ec") == 0)){
      bool ecflag = strcmp(arg[iarg], "ec")==0;
      int arg_count = ecflag ? 4 : 3;
      if(iarg+arg_count > narg)
        error->universe_all(FLERR, "Illegal oniom command");

      // arg[iarg+1] is low-level slave partition
      int low = force->inumeric(FLERR, arg[iarg+1])-1;
      if((low < 1 ) || (low>=universe->nworlds)){
        error->universe_all(FLERR, "Invalid partition in oniom command");
      }

      // arg[iarg+2] is high-level slave partition
      int high = force->inumeric(FLERR, arg[iarg+2])-1;
      if((high < 1 ) || (high>=universe->nworlds)){
        error->universe_all(FLERR, "Invalid partition in oniom command");
      }

      // arg[iarg+3] is mc-group
      int mc_group = group->find(arg[iarg+3]);
      if(mc_group == -1)
        error->all(FLERR,"Could not find oniom mc group ID");
      if(group->dynamic[mc_group])
        error->all(FLERR,"MC group must be static");
      bigint temp_mc = group->count(mc_group);
      if(temp_mc > MAXSMALLINT){
          error->all(FLERR,"Too many MC atoms for ONIOM calculation");
      }
      int mc_nat = static_cast<int>(temp_mc);

      // arg[iarg+4] is ec-group (does not need to exist on slaves)
      // TODO: check for overlap with mc_group
      int ec_group = -1;
      int ec_nat = -1;
      if(ecflag){
        ec_group = group->find(arg[iarg+4]);
        if(!iworld && (ec_group == -1)){
          error->all(FLERR,"Could not find oniom ec group ID");
        }else if(ec_group != -1){
          if(ec_group == mc_group)
            error->all(FLERR, "Electrostatic and mechanical groups overlap");
          // TODO: may remove this in the future
          bigint temp_ec = group->count(ec_group);
          if(temp_ec > MAXSMALLINT){
              error->all(FLERR,"Too many EC atoms for ONIOM calculation");
          }
          ec_nat = static_cast<int>(temp_mc);
        }
      }

      // push information to fix and create intercommunicators
      if(!iworld){
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
        int lowmode = FixONIOM::MINUS | (ecflag ? FixONIOM::EC : FixONIOM::MC);
        int highmode = FixONIOM::PLUS | (ecflag ? FixONIOM::EC : FixONIOM::MC);
        fix->connections.push_back({lowmode, mc_group, mc_nat, ec_group, ec_nat});
        setupMPI(low, fix->connections.back().comm);
        fix->connections.push_back({highmode, mc_group, mc_nat, ec_group, ec_nat});
        setupMPI(high, fix->connections.back().comm);
      }else if((iworld == (low)) || (iworld == (high))){
        if(!fix->connections.empty()){
          error->universe_one(FLERR, "Same partition used in multiple oniom regions");
        }
        fix->connections.push_back({0, mc_group, mc_nat, ec_group, ec_nat});
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
                             &fix->connections.back().comm);
      }
      iarg += arg_count+1;

    /********************
     * Verbosity
     ********************/
    }else if(strcmp(arg[iarg], "verbose") == 0){
      fix->verbose = force->inumeric(FLERR, arg[iarg+1]);
      iarg += 2;
    }else{
      error->universe_all(FLERR, "Illegal oniom command");
    }
  }
  return fix;
}
