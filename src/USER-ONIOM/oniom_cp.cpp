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

#include "oniom_cp.h"

using namespace LAMMPS_NS;

OniomCP::OniomCP(LAMMPS *l, int narg, char **arg)
    : Oniom{l, narg, arg}
{}

OniomCP::~OniomCP()
{}

int OniomCP::get_capabilities()
{
    return MD|Slave;
}
