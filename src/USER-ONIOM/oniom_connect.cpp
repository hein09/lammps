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

#include <algorithm>

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
  Oniom *oniom = update->oniom;

  int iarg{0};
  while(iarg < narg){
    // parse new connections
    if(!strcmp(arg[iarg], "couple")){
      int out = force->inumeric(FLERR, arg[iarg+2]);
      if((out <= 0 ) || (out > universe->nworlds))
        error->all(FLERR, "Invalid partition in oniom_connect command");

      int in_a = force->inumeric(FLERR, arg[iarg+3]);
      if((in_a <= 0 ) || (in_a > universe->nworlds) || (in_a == out))
        error->all(FLERR, "Invalid partition in oniom_connect command");

      int in_b = force->inumeric(FLERR, arg[iarg+4]);
      if((in_b <= 0 ) || (in_b > universe->nworlds) || (in_b == in_a) || (in_b == out))
        error->all(FLERR, "Invalid partition in oniom_connect command");

      oniom->couplings[arg[iarg+1]] = {out-1, in_a-1, in_b-1};
      if(universe->iworld == out-1){
        // am outer layer, add inner pair
        oniom->inner.push_back(&*oniom->couplings.find(arg[iarg+1]));
      }else if((universe->iworld == in_a-1) || (universe->iworld == in_b-1)){
        // am inner layer, add outer
        if(oniom->outer)
          // should also catch multiple assignments of one inner layer
          error->all(FLERR,
            "Invalid oniom_connect command: assigning two outer layers to an inner layer");
        oniom->outer = &*oniom->couplings.find(arg[iarg+1]);
      }

      iarg+=5;
    // reset connections
    }else if(!strcmp(arg[iarg], "clear")){
      if(!strcmp(arg[iarg+1], "all")){
        oniom->couplings.clear();
        oniom->outer = nullptr;
        oniom->inner.clear();
        iarg+=2;
      }else{
        // disable connection to master
        if(oniom->outer->first == arg[iarg+1]){
            oniom->outer = nullptr;
        }
        // disable connection to children
        auto pos = std::find_if(oniom->inner.begin(), oniom->inner.end(),
                   [&](Oniom::Coupling* i){return i->first == arg[iarg+1];});
        if(pos != oniom->inner.end()){
          oniom->inner.erase(pos);
        }
        // delete couplings
        oniom->couplings.erase(arg[iarg+1]);
      }
    }else{
      error->all(FLERR, "Illegal oniom_connect command");
    }
  }
  int me;
  MPI_Comm_rank(world, &me);
  // pretty-print topology
  if(me == 0){
    const auto& out = oniom->outer;
    const auto& in = oniom->inner;
    int off = in.empty() ? 3 : 1 + (int)in.size()*3;
    auto print_topo = [&](FILE*& f){
      if (f) fprintf(f, "Local ONIOM Topology:\n");
      if (oniom->outer){
        if (f) fprintf(f, "%*d (%s)\n%*c\n%*d (me)\n", off, out->second.master+1,
                            out->first.c_str(), off, '|', off, universe->iworld+1);
      }else{
        if (f) fprintf(f, "%*d(me, outermost)\n", off, universe->iworld+1);
      }
      if(!oniom->inner.empty()){
        for(size_t i=0; i<in.size(); ++i){
          if(i < (in.size()/2)){
            if (f) fprintf(f, "  /  /");
          }else if((in.size() % 2) && (i == (in.size()/2))){
            if (f) fprintf(f, "  /  \\");
          }else{
            if (f) fprintf(f, "  \\  \\");
          }
        }
        if (f) fprintf(f, "\n");
        for(size_t i=0; i<in.size(); ++i){
          if (f) fprintf(f, "%3d%3d", in[i]->second.low_slave+1, in[i]->second.high_slave+1);
        }
        if (f) fprintf(f, "\n");
        for(size_t i=0; i<in.size(); ++i){
          if (f) fprintf(f, "%6.6s", in[i]->first.c_str());
        }
        if (f) fprintf(f, "\n");
      }
    };
    if (screen) print_topo(screen);
    if (logfile) print_topo(logfile);
  }

}
