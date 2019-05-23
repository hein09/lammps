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

#include <mpi.h>
#include "style_oniom.h"
#include "verlet_oniom.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

VerletOniom::VerletOniom(LAMMPS *l, int narg, char **arg) :
    Verlet(l, narg, arg)
{
}

VerletOniom::~VerletOniom()
{
}

/* ----------------------------------------------------------------------
   initialization before run
------------------------------------------------------------------------- */

void VerletOniom::init()
{
    if(strcmp(update->oniom_style, "lammps")){
        update->oniom->init();
    }else{
        Verlet::init();
    }
}

/* ----------------------------------------------------------------------
   setup before run
------------------------------------------------------------------------- */

void VerletOniom::setup(int flag)
{
    if(strcmp(update->oniom_style, "lammps")){
        update->oniom->setup(flag);
    }else{
        Verlet::setup(flag);
    }
}

/* ----------------------------------------------------------------------
   setup without output
   flag = 0 = just force calculation
   flag = 1 = reneighbor and force calculation
------------------------------------------------------------------------- */

void VerletOniom::setup_minimal(int flag)
{
    if(strcmp(update->oniom_style, "lammps")){
        update->oniom->setup_minimal(flag);
    }else{
        Verlet::setup_minimal(flag);
    }
}

/* ----------------------------------------------------------------------
   run for N steps
------------------------------------------------------------------------- */

void VerletOniom::run(int n)
{
    if(strcmp(update->oniom_style, "lammps")){
        update->oniom->run(n);
    }else{
        Verlet::run(n);
    }
}

/* ---------------------------------------------------------------------- */

void VerletOniom::cleanup()
{
    if(strcmp(update->oniom_style, "lammps")){
        update->oniom->cleanup();
    }else{
        Verlet::cleanup();
    }
}

void VerletOniom::reset_dt()
{
    if(strcmp(update->oniom_style, "lammps")){
        update->oniom->reset_dt();
    }else{
        Verlet::reset_dt();
    }
}

bigint VerletOniom::memory_usage()
{
    if(strcmp(update->oniom_style, "lammps")){
        update->oniom->memory_usage();
    }else{
        return Verlet::memory_usage();
    }
}
