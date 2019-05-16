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
#include "verlet_oniom.h"
#include "style_oniom.h"
#include "error.h"
#include "universe.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

VerletOniom::VerletOniom(LAMMPS *l, int narg, char **arg) :
    Verlet(l, narg, arg)
{
    if (narg < 1) error->all(FLERR, "Illegal run_style verlet/oniom command");

    OniomCreatorMap oniom_map{};
#define ONIOM_CLASS
#define OniomStyle(key,Class) \
    oniom_map[#key] = &oniom_creator<Class>;
#include "style_oniom.h"
#undef OniomStyle
#undef ONIOM_CLASS

    // "lammps" falls back to default verlet-integrator
    if(strcmp(arg[0], "lammps") == 0){
        return;
    }

    auto pos = oniom_map.find(arg[0]);
    if(pos != oniom_map.end()){
        oniom = pos->second(lmp, narg-1, &arg[1]);
    }else{
        error->all(FLERR, "Illegal run_style verlet/oniom command");
    }
}
template <typename T>
Oniom* VerletOniom::oniom_creator(LAMMPS *l, int narg, char**arg)
{
    return new T(l, narg, arg);
}

VerletOniom::~VerletOniom()
{
    if(oniom) delete oniom;
}

/* ----------------------------------------------------------------------
   initialization before run
------------------------------------------------------------------------- */

void VerletOniom::init()
{
    verify_connections();
    if(oniom){
        oniom->init();
    }else{
        Verlet::init();
    }
}

/* ----------------------------------------------------------------------
   setup before run
------------------------------------------------------------------------- */

void VerletOniom::setup(int flag)
{
    verify_connections();
    if(oniom){
        oniom->setup(flag);
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
    verify_connections();
    if(oniom){
        oniom->setup_minimal(flag);
    }else{
        Verlet::setup_minimal(flag);
    }
}

/* ----------------------------------------------------------------------
   run for N steps
------------------------------------------------------------------------- */

void VerletOniom::run(int n)
{
    if(oniom){
        oniom->run(n);
    }else{
        Verlet::run(n);
    }
}

/* ---------------------------------------------------------------------- */

void VerletOniom::cleanup()
{
    if(oniom){
        oniom->cleanup();
    }else{
        Verlet::cleanup();
    }
}

void VerletOniom::reset_dt()
{
    if(oniom){
        oniom->reset_dt();
    }else{
        Verlet::reset_dt();
    }
}

bigint VerletOniom::memory_usage()
{
    if(oniom){
        return oniom->memory_usage();
    }else{
        return Verlet::memory_usage();
    }
}

/* ----------------------------------------------------------------------
    Make sure partitions are connected correctly
   ---------------------------------------------------------------------- */

void VerletOniom::verify_connections()
{
    uint8_t tmp;
    MPI_Status status;
    // wait for inner area
    if(inner_a != -1){
        MPI_Recv(&tmp, 1, MPI_BYTE,
                 universe->root_proc[inner_a],
                 universe->iworld,
                 universe->uworld,
                 &status);
    }
    if(inner_b != -1){
        MPI_Recv(&tmp, 1, MPI_BYTE,
                 universe->root_proc[inner_b],
                 universe->iworld,
                 universe->uworld,
                 &status);
    }
    // check if master, else notify outer area
    if(outer == -1){
        master = true;
    }else{
        MPI_Send(&tmp, 1, MPI_BYTE,
                 universe->root_proc[outer],
                 outer,
                 universe->uworld);
    }
}
