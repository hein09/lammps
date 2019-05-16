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
    VO->inner_a = -1;
    VO->inner_b = -1;

    // parse new connections
    int iarg{0};
    while(iarg < narg){
        if(strcmp(arg[iarg], "outer") == 0){
            int tmp = force->inumeric(FLERR, arg[iarg+1]);
            if((tmp <= 0 ) || (tmp > universe->nworlds))
                error->all(FLERR, "Invalid partition in oniom_connect command");
            VO->outer = tmp-1;
            iarg+=2;
        }else if(strcmp(arg[iarg], "inner") == 0){
            int tmp = force->inumeric(FLERR, arg[iarg+1]);
            if((tmp <= 0 ) || (tmp > universe->nworlds))
                error->all(FLERR, "Invalid partition in oniom_connect command");
            VO->inner_a = tmp-1;
            tmp = force->inumeric(FLERR, arg[iarg+2]);
            if((tmp <= 0 ) || (tmp > universe->nworlds))
                error->all(FLERR, "Invalid partition in oniom_connect command");
            VO->inner_b = tmp-1;
            iarg+=3;
        }else{
            error->all(FLERR, "Illegal oniom_connect command");
        }
    }
}
