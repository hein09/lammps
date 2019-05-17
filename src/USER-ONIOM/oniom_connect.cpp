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

#include "force.h"
#include "update.h"
#include "universe.h"
#include "domain.h"
#include "error.h"

#include "verlet_oniom.h"
#include "oniom_connect.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

OniomConnect::OniomConnect(LAMMPS *l) : Pointers{l} {}

OniomConnect::~OniomConnect(){}

/* ----------------------------------------------------------------------
   create logical connections to other oniom-partitions
------------------------------------------------------------------------- */

void OniomConnect::command(int narg, char **arg)
{
  if (universe->nworlds == 1)
    error->all(FLERR, "Must have more than one partition to connect oniom-groups");
  if (domain->box_exist == 0)
    error->all(FLERR, "oniom_connect command before simulation box is defined");
  VerletOniom* VO = dynamic_cast<VerletOniom*>(update->integrate);
  if ((strcmp(update->integrate_style, "verlet/oniom") != 0) || (VO == nullptr))
    error->all(FLERR, "oniom_connect without oniom-capable run_style");

  // reset connections
  VO->outer = -1;
  VO->inner.clear();

  // parse new connections
  int iarg{0};
  while(iarg < narg){
    if(strcmp(arg[iarg], "couple") == 0){
      int out = force->inumeric(FLERR, arg[iarg+1]);
      if((out <= 0 ) || (out > universe->nworlds))
        error->all(FLERR, "Invalid partition in oniom_connect command");

      int in_a = force->inumeric(FLERR, arg[iarg+2]);
      if((in_a <= 0 ) || (in_a > universe->nworlds) || (in_a == out))
        error->all(FLERR, "Invalid partition in oniom_connect command");

      int in_b = force->inumeric(FLERR, arg[iarg+3]);
      if((in_b <= 0 ) || (in_b > universe->nworlds) || (in_b == in_a) || (in_b == out))
        error->all(FLERR, "Invalid partition in oniom_connect command");

      if(universe->iworld == out-1){
        // am outer layer, add inner pair
        VO->inner.push_back(std::pair<int,int>{in_a-1, in_b-1});
      }else if((universe->iworld == in_a-1) || (universe->iworld == in_b-1)){
        // am inner layer, add outer
        if(VO->outer != -1)
          // should also catch multiple assignments of one inner layer
          error->all(FLERR,
            "Invalid oniom_connect command: assigning two outer layers to an inner layer");
        VO->outer = out-1;
      }

      iarg+=4;
    }else{
      error->all(FLERR, "Illegal oniom_connect command");
    }
  }
  int me;
  MPI_Comm_rank(world, &me);
  // pretty-print topology
  if(me == 0){
    if (screen) fprintf(screen, "Local ONIOM Topology:\n");
    if (logfile) fprintf(logfile, "Local ONIOM Topology:\n");
    const auto& in = VO->inner;
    int off = in.size() * 6 - 1;
    if (VO->outer != -1){
      if (screen) fprintf(screen, "%*d\n%*c\n%*d(me)\n", off, VO->outer+1, off, '|', off, universe->iworld+1);
      if (logfile) fprintf(logfile, "%*d\n%*c\n%*d(me)\n", off, VO->outer+1, off, '|', off, universe->iworld+1);
    }else{
      if (screen) fprintf(screen, "%*d(me, outermost)\n", off, universe->iworld+1);
      if (logfile) fprintf(logfile, "%*d(me, outermost)\n", off, universe->iworld+1);
    }
    if(!VO->inner.empty()){
      for(size_t i=0; i<in.size(); ++i){
        if(i < (in.size()/2)){
          if (screen) fprintf(screen, "  /  /");
          if (logfile) fprintf(logfile, "  /  /");
        }else if((in.size() % 2) && (i == (in.size()/2))){
          if (screen) fprintf(screen, "  /  \\");
          if (logfile) fprintf(logfile, "  /  \\");
        }else{
          if (screen) fprintf(screen, "  \\  \\");
          if (logfile) fprintf(logfile, "  \\  \\");
        }
      }
      if (screen) fprintf(screen, "\n");
      if (logfile) fprintf(logfile, "\n");
      for(size_t i=0; i<in.size(); ++i){
        if (screen) fprintf(screen, "%3d%3d", in[i].first+1, in[i].second+1);
        if (logfile) fprintf(logfile, "%3d%3d", in[i].first+1, in[i].second+1);
      }
      if (screen) fprintf(screen, "\n");
      if (logfile) fprintf(logfile, "\n");
    }
  }

}
