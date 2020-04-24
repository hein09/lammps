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

#include "fix_collect.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "group.h"
#include "universe.h"

#include <algorithm>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

// static member
std::map<std::string, std::weak_ptr<FixCollect::collection_t>> FixCollect::collections;

FixCollect::FixCollect(LAMMPS *l, int narg, char **arg):
    Fix{l, narg, arg}
{
  // initialize buffers for variable-size collectives
  recv_count_buf.resize(universe->nprocs);
  displs_buf.resize(universe->nprocs);
}

double FixCollect::memory_usage()
{
  double bytes = sizeof(FixCollect);

  bytes += displs_buf.capacity() * sizeof(int);
  bytes += recv_count_buf.capacity() * sizeof(int);

  // NOTE: map-sizes not correct, but node_type is only available from C++17 onwards
  for(auto& entry: collections){
    bytes += sizeof(decltype(collections)::value_type);
    bytes += entry.first.size() * sizeof(std::string::value_type);
    auto coll = entry.second.lock();
    if(coll){
      bytes += sizeof(decltype(*coll));
      bytes += coll->idx2tag.capacity() * sizeof(tagint);
      bytes += coll->buffer.capacity() * sizeof(commdata_t);
      bytes += coll->tag2idx.size() * sizeof(decltype(coll->tag2idx)::value_type);
    }
  }
  return bytes;
}

std::shared_ptr<FixCollect::collection_t> FixCollect::get_collection(char *group_name)
{
  // if someone else already created a collection for this group, return this
  if(!collections[group_name].expired()){
    return collections[group_name].lock();
  }
  // else create a new collection
  int group_id = group->find(group_name);
  if(group_id < 0){
    char msg[50];
    sprintf(msg, "Invalid group %s in fix %s", group_name, style);
    error->all(FLERR, msg);
  }
  bigint nat = group->count(igroup);
  if(nat > MAXSMALLINT){
    char msg[70];
    sprintf(msg, "Too many atoms in group %s for fix %s", group_name, style);
    error->all(FLERR, msg);
  }
  auto coll = std::make_shared<collection_t>();
  coll->nat = nat;
  coll->group_id = group_id;
  collect_tags(group->bitmask[group_id], *coll);
  collections[group_name] = coll;
  return coll;
}

void FixCollect::collect_tags(int groupmask, collection_t& coll)
{
  // collect process-local tags
  std::vector<tagint> tmp{};
  tmp.reserve(static_cast<size_t>(coll.nat));
  for(int i=0; i<atom->nlocal; ++i){
    if(atom->mask[i] & groupmask)
      tmp.push_back(atom->tag[i]);
  }

  // gather number of tags
  int send_count = tmp.size() * sizeof(decltype (tmp)::value_type);
  MPI_Allgather(&send_count, 1, MPI_INT, recv_count_buf.data(), 1, MPI_INT, world);

  // construct displacement
  displs_buf[0] = 0;
  for(int i=1; i<universe->nprocs; ++i){
    displs_buf[i] = displs_buf[i-1]+recv_count_buf[i-1];
  }

  // gather and sort tags
  coll.idx2tag.resize(static_cast<size_t>(coll.nat));
  MPI_Allgatherv(tmp.data(), send_count, MPI_BYTE,
                 coll.idx2tag.data(), recv_count_buf.data(), displs_buf.data(), MPI_BYTE, world);
  std::sort(coll.idx2tag.begin(), coll.idx2tag.end());

  // assign tags to qm-id
  for(int i=0; i<coll.nat; ++i){
    coll.tag2idx[coll.idx2tag[i]] = i;
  }
}

void FixCollect::collect_lammps(collection_t &coll, double **source)
{
  const int *const mask = atom->mask;
  const tagint * const tag = atom->tag;
  auto& buffer = coll.buffer;

  // collect process-local atom positions
  buffer.clear();
  buffer.reserve(static_cast<size_t>(coll.nat));
  for(int i=0; i<atom->nlocal; ++i){
    if(mask[i] & groupbit){
      buffer.push_back({tag[i], {source[i][0], source[i][1], source[i][2]}});
    }
  }

  // collect each proc's number of atoms
  int send_count = buffer.size() * sizeof(commdata_t);
  MPI_Gather(&send_count, 1, MPI_INT, recv_count_buf.data(), 1, MPI_INT, 0, world);

  // construct displacement
  displs_buf[0] = 0;
  for(int i=1; i<universe->nprocs; ++i){
    displs_buf[i] = displs_buf[i-1]+recv_count_buf[i-1];
  }

  // collect in buffer on proc 0
  buffer.resize(static_cast<size_t>(coll.nat));
  MPI_Gatherv((comm->me == 0) ? MPI_IN_PLACE : buffer.data(),
              send_count, MPI_BYTE,
              buffer.data(), recv_count_buf.data(), displs_buf.data(), MPI_BYTE, 0, world);

}
