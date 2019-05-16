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

using namespace LAMMPS_NS;

Oniom::Oniom(LAMMPS *l, int, char**)
    : Pointers{l}
{}

Oniom::~Oniom()
{}

void Oniom::reset_dt()
{}

bigint Oniom::memory_usage()
{
    return 0;
}

uint8_t Oniom::get_capabilities()
{
    return 0;
}
