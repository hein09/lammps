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

#include "fix_qe.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "universe.h"
#include "update.h"

#include <algorithm>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

bool FixQE::initialized = false;

FixQE::FixQE(LAMMPS *l, int narg, char **arg):
    FixCollect{l, narg, arg}
{
  // this fix contributes to the systems energy
  thermo_energy = 1;
  scalar_flag = 1;

  /* convert energies and forces
   * energy is returned as eV
   * forces are returned as eV/Angstrom
   * NOTE: other formats also need additional scaling for positions
   */
  if (strcmp(update->unit_style, "metal") == 0){
    // need eV and eV / Angstrom
    escale = fscale = 1.0;
  }else if (strcmp(update->unit_style, "real") == 0){
    // need kcal/mol and kcal/mol / Angstrom
    escale = fscale = 23.060549;
  }else{
    char msg[50];
    sprintf(msg, "fix %s requires real or metal units", style);
    error->all(FLERR, msg);
  }

  // require tags for communication
  if (atom->tag_enable == 0) {
    char msg[30];
    sprintf(msg, "fix %s requires atom IDs", style);
    error->all(FLERR, msg);
  }

  // limit to one instance per partition, QE probably has too much global state
  if (initialized)
      error->all(FLERR, "Only one QuantumEspresso-based fix allowed at a time");

  // parse arguments
  if (narg < 4) {
    char msg[30];
    sprintf(msg, "Illegal fix %s command", style);
    error->all(FLERR, msg);
  }

  // request collection for main qm group
  qm_coll = get_collection(arg[1]);

  // require regular input-file to setup QE
  inp_file = arg[3];
  // redirect stdout of QE
  out_file = "log.qe_";
  out_file += &style[3];
  if(universe->existflag){
    out_file += '.' + std::to_string(universe->iworld);
  }

  // optional additional arguments
  int iarg = 4;
  while(iarg < narg-1){
    // QE Parallelization flags
    if(strcmp(arg[iarg], "npool") == 0){
      mpi_partitions[0] = force->inumeric(FLERR, arg[iarg+1]);
      iarg += 2;
    }else if(strcmp(arg[iarg], "ntask") == 0){
      mpi_partitions[1] = force->inumeric(FLERR, arg[iarg+1]);
      iarg += 2;
    }else if(strcmp(arg[iarg], "nband") == 0){
      mpi_partitions[2] = force->inumeric(FLERR, arg[iarg+1]);
      iarg += 2;
    }else if(strcmp(arg[iarg], "ndiag") == 0){
      mpi_partitions[3] = force->inumeric(FLERR, arg[iarg+1]);
      iarg += 2;
    // Electrostatic coupling
    }else if(strcmp(arg[iarg], "ec") == 0){
      if(comm->me == 0)
        error->warning(FLERR, "Electrostatic coupling is not efficient");
      ec_coll = get_collection(arg[iarg+1]);
      if(ec_coll->group_id == qm_coll->group_id){
        error->all(FLERR, "EC group must differ from main group of fix qe/pw");
      }
      iarg += 2;
    // Link atoms
    }else if(strcmp(arg[iarg], "link") == 0){
      int nlink = force->inumeric(FLERR, arg[iarg+1]);
      linkatoms.resize(nlink);
      for(int i=0; i<nlink; ++i){
        linkatoms[i] = {force->inumeric(FLERR, arg[iarg+2+2*i]),
                        -1,
                        force->numeric(FLERR, arg[iarg+3+2*i])};
      }
      auto last = std::unique(linkatoms.begin(), linkatoms.end(),
                              [](const linkgroup& l, const linkgroup& r)
                              {return l.link_atom == r.link_atom;});
      if(last != linkatoms.end()){
        error->all(FLERR, "Duplicate entries in list of link atoms");
      }
      iarg += 2+2*nlink;
    // Scheme for EC-handling of Link atoms
    }else if(strcmp(arg[iarg], "scheme") == 0){
      if(strcmp(arg[iarg+1], "Z1") == 0){
        charge_scheme = charge_t::Z1;
      }else if(strcmp(arg[iarg+1], "Z2") == 0){
        charge_scheme = charge_t::Z2;
      }else if(strcmp(arg[iarg+1], "Z3") == 0){
        charge_scheme = charge_t::Z3;
      }else if(strcmp(arg[iarg+1], "RCD") == 0){
        charge_scheme = charge_t::RCD;
      }else if(strcmp(arg[iarg+1], "CS") == 0){
        charge_scheme = charge_t::CS;
      }else{
        error->all(FLERR, "Invalid charge redistribution scheme");
      }
      iarg += 2;
    }else{
      char msg[30];
      sprintf(msg, "Illegal fix %s command", style);
      error->all(FLERR, msg);
    }
  }

  // calculate number of EC atoms and initialize_ec_buf
  if(ec_coll){
    auto nec = ec_coll->nat;
    if(comm->me == 0){
      switch(charge_scheme){
      case charge_t::Z1:
        break;
      case charge_t::Z2:
      case charge_t::Z3:
        nec -= chargeatoms.size();
        break;
      case charge_t::RCD:
        nec += chargeatoms.size();
        break;
      case charge_t::CS:
        nec += 2*chargeatoms.size();
        break;
      default:
        error->all(FLERR, "Invalid charge redistribution scheme");
      }
    }
    // TODO: preallocate space for additional atoms
  }

  // validate link and modified-charge atoms
  if (!linkatoms.empty()) {
    if(!(atom->molecular)){
      char msg[50];
      sprintf(msg, "Fix %s: Link atoms can only be used when bonds are declared", style);
      error->all(FLERR, msg);
    }
    // make sure map from global to local ids is initialized
    // TODO: ghost in atom->map could be problematic, could also be much slower than needed,
    // evaluate other solution
    bool mapflag{false};
    if(!atom->map_style){
      mapflag = true;
      atom->map_init();
      atom->map_set();
    }
    for(auto& link: linkatoms){
      // check if link atom is in QM group
      auto idx = atom->map(link.link_atom);
      if((idx != -1) && !(atom->mask[idx] & groupbit)){
          char msg[50];
          sprintf(msg, "Fix %s: Link atom %d is not included in QM group.",
                  style, link.link_atom);
          error->one(FLERR, msg);
      }
      auto neighbors = get_bound_neighbors(link.link_atom, qm_coll->idx2tag, 0);
      if(neighbors.empty()){
        char msg[50];
        sprintf(msg, "Fix %s: Link atom %d has no neighbors in QM group.",
                style, link.link_atom);
        error->all(FLERR, msg);
      }else if(neighbors.size()>1){
        char msg[50];
        sprintf(msg, "Fix %s: Link atom %d is bound to multiple QM atoms.",
                style, link.link_atom);
        error->all(FLERR, msg);
      }else{
        link.qm_atom = neighbors.front();
      }
    }
    // if needed, find atoms whose charges need to be modified
    if(ec_coll){
      const auto ecmask = group->bitmask[ec_coll->group_id];
      if(charge_scheme == charge_t::Z1){
        // only ignore charge of link-atoms,
        // which is assumed to be excluded from EC group
      }else if(charge_scheme == charge_t::Z2){
        // ignore charge of 1-2 neighbors
        for(const auto& link: linkatoms){
          auto neighbors = get_bound_neighbors(link.link_atom, ec_coll->idx2tag, 0);
          for(const auto& n: neighbors){
            chargeatoms.push_back({n, -1});
          }
        }
      }else if(charge_scheme == charge_t::Z3){
        // ignore charge of 1-2 and 1-3 neighbors
        for(const auto& link: linkatoms){
          auto neighbors = get_bound_neighbors(link.link_atom, ec_coll->idx2tag, 1);
          for(const auto& n: neighbors){
            chargeatoms.push_back({n, -1});
          }
        }
      }else if((charge_scheme == charge_t::RCD) || (charge_scheme == charge_t::CS)){
        // create virtual point charges depending on 1-2 neighbors
        for(const auto& link: linkatoms){
          auto neighbors = get_bound_neighbors(link.link_atom, ec_coll->idx2tag, 0);
          for(const auto& n: neighbors){
            chargeatoms.push_back({n, link.link_atom});
          }
        }
      }
    }
    // de-init map if we initialized it
    if(mapflag){
      atom->map_delete();
      atom->map_style = 0;
    }
  }
}

void FixQE::distribute_forces()
{
  // redistribute forces of link atoms
  auto& buffer = qm_coll->buffer;
  const auto& tag2idx = qm_coll->tag2idx;
  for(const auto& l: linkatoms){
    auto& forceL = buffer[tag2idx.at(l.link_atom)].x;
    auto& forceQ = buffer[tag2idx.at(l.qm_atom)].x;
    forceQ[0] += (1-l.ratio)*forceL[0];
    forceQ[1] += (1-l.ratio)*forceL[1];
    forceQ[2] += (1-l.ratio)*forceL[2];
    forceL[0] *= l.ratio;
    forceL[1] *= l.ratio;
    forceL[2] *= l.ratio;
  }

  // reassign forces to local lammps atoms and mark buffer contents
  distribute_forces_impl(*qm_coll);
  qm_coll->state = collection_t::State::Forces;
  if(ec_coll){
    distribute_forces_impl(*ec_coll);
    ec_coll->state = collection_t::State::Forces;
  }
}

void FixQE::distribute_forces_impl(collection_t& coll)
{
  auto& buffer = coll.buffer;
  const tagint * const tag = atom->tag;
  double ** f = atom->f;
  for(int i=0; i<atom->nlocal; ++i){
    for(auto& dat: buffer){
      if(tag[i] == coll.idx2tag[dat.idx]){
        f[i][0] += dat.x[0];
        f[i][1] += dat.x[1];
        f[i][2] += dat.x[2];
      }
    }
  }
}

void FixQE::collect_positions()
{
  // collect positions of lammps processes if needed
  if(qm_coll->state != collection_t::State::Positions){
    // TODO: benchmark this
//    gather_root(*qm_coll, atom->x);
//    auto &buffer = qm_coll->buffer;
//    MPI_Bcast(buffer.data(), buffer.size() * sizeof(commdata_t), MPI_BYTE, 0, world);
//    gather_all(*qm_coll, atom->x);
    gather_all_inplace(*qm_coll, atom->x);
    com_to_origin(*qm_coll);
    qm_coll->state = collection_t::State::Positions;
  }
  // modify positions of link atoms
  auto& buffer = qm_coll->buffer;
  const auto& tag2idx = qm_coll->tag2idx;
  for(const auto& l: linkatoms){
    auto& posL = buffer[tag2idx.at(l.link_atom)].x;
    auto& posQ = buffer[tag2idx.at(l.qm_atom)].x;
    posL[0] = posQ[0] + l.ratio * (posL[0]-posQ[0]);
    posL[1] = posQ[1] + l.ratio * (posL[1]-posQ[1]);
    posL[2] = posQ[2] + l.ratio * (posL[2]-posQ[2]);
  }

  // collect positions of EC atoms if needed
  if(ec_coll && ec_coll->state != collection_t::State::Positions){
    gather_all(*ec_coll, atom->x);
    ec_coll->state = collection_t::State::Positions;
  }
}

std::vector<tagint> FixQE::get_bound_neighbors(tagint tag, const std::vector<tagint> &grouptags, int dist)
{
  // TODO: this depends on atom->special being indexed by local id and only being available locally. True?
  std::vector<tagint> neighs;
  int idx = -1;
  const auto nlocal = atom->nlocal;
  const auto tags = atom->tag;
  for(int i=0; i < nlocal; ++i){
    if(tags[i] == tag){
      idx = i;
      break;
    }
  }
  int buf[2] = {0, 0}; // {source-process, nspecial}
  if(idx != -1){
    // we own the requested source atom, broadcast bound atoms
    buf[0] = comm->me;
    buf[1] = atom->nspecial[idx][dist];
    neighs.insert(neighs.begin(), atom->special[idx], atom->special[idx]+buf[1]);
  }
  MPI_Allreduce(MPI_IN_PLACE, &buf, 2, MPI_INT, MPI_SUM, world);
  neighs.resize(buf[1]);
  MPI_Bcast(neighs.data(), buf[1], MPI_LMP_TAGINT, buf[0], world);
  // check if bound atoms are part of group
  std::vector<tagint> groupneighs;
  std::sort(neighs.begin(), neighs.end());
  std::set_intersection(neighs.begin(), neighs.end(),
                        grouptags.begin(), grouptags.end(),
                        std::back_inserter(groupneighs));
  return groupneighs;
}
