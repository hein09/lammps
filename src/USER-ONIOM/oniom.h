/* -*- c++ -*- -------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_ONIOM_H
#define LMP_ONIOM_H

#include "pointers.h"

namespace LAMMPS_NS {

class Oniom : public Pointers {
public:
    Oniom(LAMMPS *l, int, char**);
    ~Oniom();
    virtual void init(){}
    virtual void setup(int){}
    virtual void setup_minimal(int){}
    virtual void run(int){}
    virtual void cleanup(){}
    virtual void reset_dt();
    virtual bigint memory_usage();

    enum Capabilities:uint8_t {MD=0x1, Min=0x2};
    virtual uint8_t get_capabilities();
};

}

#endif // LMP_ONIOM_H
