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

#include <cstring>

#include "oniom.h"
#include "fix_oniom.h"

#include "run.h"
#include "minimize.h"

#include "timer.h"

#include "error.h"
#include "modify.h"

using namespace LAMMPS_NS;

/****************************************
 * usage:
 *
 * use one of:
 * > oniom run <runargs> <keyword args ...>
 * > oniom minimize <minargs> <keyword args ...>
 *
 * keyword args can be one or multiple of:
 * group <group-ID> <pair of partitons>
 * verbose 0/1/2
 ****************************************/
void Oniom::command(int narg, char **arg)
{
  // scan for first oniom command
  int ncmdarg = -1;
  for(int i=0; i<narg; ++i){
    if(!strcmp(arg[i], "group") ||
       !strcmp(arg[i], "plus") ||
       !strcmp(arg[i], "minus") ||
       !strcmp(arg[i], "verbose")){
        ncmdarg = i;
        break;
    }
  }
  if(ncmdarg == -1){
      error->universe_one(FLERR, "Illegal oniom command: oniom arguments missing");
      // could be universe_all, but would prob. break if user prefixes the oniom call with a partition cmd
  }

  // create oniom fix
  auto fix = mkFix(narg-ncmdarg, &arg[ncmdarg]);

  // execute command
  if(strcmp(arg[0], "run")==0){
    run(ncmdarg-1, &arg[1], fix);
  }else if(strcmp(arg[0], "minimize")==0){
    minimize(ncmdarg-1, &arg[1], fix);
  }else{
    char msg[50];
    sprintf(msg, "Illegal oniom command: command %s not supported", arg[0]);
    error->universe_one(FLERR, msg);
  }

  // delete oniom fix
  modify->delete_fix("ONIOM");
}

void Oniom::run(int narg, char **arg, FixONIOM *fix)
{
  fix->minimize = false;
  // perform calculation
  Run run{lmp};
  run.command(narg, arg);
}

void Oniom::minimize(int narg, char **arg, FixONIOM *fix)
{
  fix->minimize = true;

  // perform calculation
  if(fix->master){
      // master performs minimization
      Minimize min{lmp};
      min.command(narg, arg);
  }else{
      // slave performs faux md with N = maxeval
      Run run{lmp};
      run.command(1, &arg[3]);
      timer->reset_timeout();
  }
}

FixONIOM* Oniom::mkFix(int narg, char **arg)
{
  // create private oniom fix
  char **fixarg = new char*[narg+3];
  fixarg[0] = (char *) "ONIOM";
  fixarg[1] = (char *) "all";
  fixarg[2] = (char *) "ONIOM";
  for(int i=0; i<narg; ++i){
      fixarg[i+3] = arg[i];
  }
  modify->add_fix(narg+3, fixarg);
  delete[] fixarg;

  return static_cast<FixONIOM*>(modify->fix[modify->nfix-1]);
}
