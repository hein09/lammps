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

#ifdef INTEGRATE_CLASS

IntegrateStyle(verlet/oniom,VerletOniom)

#else

#ifndef LMP_VERLET_ONIOM_H
#define LMP_VERLET_ONIOM_H

#include "oniom.h"
#include "verlet.h"

namespace LAMMPS_NS {

class VerletOniom : public Verlet {
    friend class OniomConnect;
public:
    VerletOniom(LAMMPS *, int, char **);
    ~VerletOniom();
    void init() override;
    void setup(int) override;
    void setup_minimal(int) override;
    void run(int) override;
    void cleanup() override;
    void reset_dt() override;
    bigint memory_usage() override;

    using OniomCreator = Oniom* (*)(LAMMPS *, int, char**);
    using OniomCreatorMap = std::map<std::string,OniomCreator>;

protected:
    Oniom* oniom{nullptr};
    bool master{true};
    int outer{-1};
    int inner_a{-1};
    int inner_b{-1};

    void verify_connections();

private:
    template <typename T> static Oniom *oniom_creator(LAMMPS *, int, char**);
};

}

#endif // LMP_VERLET_ONIOM_H
#endif

/* ERROR/WARNING messages:

  TODO

*/
