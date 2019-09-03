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

#include "run.h"
#include "minimize.h"
#include "integrate.h"

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
  if(strcmp(arg[0], "run")==0){
    run(narg-1, &arg[1]);
  }else if(strcmp(arg[0], "minimize")==0){
    minimize(narg-1, &arg[1]);
  }else{
    error->universe_all(FLERR, "Illegal oniom command");
  }
}

void Oniom::run(int narg, char **arg)
{
  // extract run call
  int nrunarg = -1;
  for(int i=0; i<narg; ++i){
    if((strcmp(arg[i], "group")==0) ||
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
  int nminarg = -1;
  for(int i=0; i<narg; ++i){
    if((strcmp(arg[i], "group")==0) ||
       (strcmp(arg[i], "verbose")==0)){
      nminarg = i;
      break;
    }
  }
  if(nminarg != 4){
    error->universe_all(FLERR, "Illegal oniom minimize command");
  }
  // create oniom fix
  auto fix = mkFix(narg-nminarg, &arg[nminarg]);
  fix->minimize = true;

  // perform calculation
  if(fix->master){
      // master performs minimization
      Minimize min{lmp};
      min.command(nminarg, arg);
  }else{
      // slave performs faux md with N = maxeval
      Run run{lmp};
      run.command(1, &arg[3]);
      timer->reset_timeout();
  }

  // delete oniom fix
  modify->delete_fix("ONIOM");
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
