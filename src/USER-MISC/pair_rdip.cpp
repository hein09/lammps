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
   Contributing authors:
      Konstantin Weber 2012-11-22
      Florian Wullschlaeger, Sebastian Gsaenger 2017 (FAU Erlangen)
      e-mail: sebastian dot gsaenger at fau dot de
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include "pair_rdip.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include <set>
#include <map>
#include <algorithm>

#include "ctype.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define MAXNEIGH 6

/* -- Vector helper functions -- */

inline double vec_dot(const double a[3], const double b[3])
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline void vec_cross(const double a[3], const double b[3], double t[3])
{
  t[0] = a[1]*b[2] - a[2]*b[1];
  t[1] = a[2]*b[0] - a[0]*b[2];
  t[2] = a[0]*b[1] - a[1]*b[0];
}

/* ---------------------------------------------------------------------- */

PairRDIP::PairRDIP(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  one_coeff = 1;
  writedata = 1;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairRDIP::~PairRDIP()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(tmap);
    memory->destroy(ntmap);
    memory->destroy(pmap);
  }
}

/* ---------------------------------------------------------------------- */
void PairRDIP::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  const int newton_pair = force->newton_pair; // always true
  const int nlocal = atom->nlocal;
  const int *type  = atom->type;
  const double * const*x = atom->x;
  double **f       = atom->f;

  const int inum         = list->inum;
  const int *ilist       = list->ilist;
  const int *numneigh    = list->numneigh;
  int **firstneigh       = list->firstneigh;

  const double cut_sq = cutoff * cutoff;

  for (int ii=0; ii < inum; ++ii) {
    const int i      = ilist[ii];
    const int itype  = type[i];
    if (tmap[itype] <= 0) continue; // skip helper and ignored atoms

    const double *xi = x[i];
    const int *jlist = firstneigh[i];
    const int jnum   = numneigh[i];
    const int nneigh = ntmap[itype].nneigh;
    const double cut_neigh = ntmap[itype].cut_neigh;

    /* r_li:
     *
     * distances between reference atom i
     * and nearest neighbors l
     */
    double r_li_len[MAXNEIGH]; //squared, only needed for sorting
    int    i_li[MAXNEIGH];
    for (int l = 0; l < nneigh; ++l) {
      r_li_len[l] = (cut_neigh*cut_neigh) + l;
      i_li[l] = -1;
    }

    /* helper-atom or same type: determine next neighbors
     *
     * keep nearest neighbors for calculation
     * of surface normals
     */
    for (int jj = 0; jj < jnum; ++jj) {
      const int j      = jlist[jj] & NEIGHMASK;
      const int jtype  = type[j];
      if ((tmap[jtype] == -1) || (jtype == itype)) {
        double r_li[3] =
          {xi[0] - x[j][0],
           xi[1] - x[j][1],
           xi[2] - x[j][2]};
        const double rsq = r_li[0]*r_li[0] + r_li[1]*r_li[1] + r_li[2]*r_li[2];

        int pos = -1;
        // check whether l is one of nearest neighbors
        for (int l=0; l<nneigh; ++l) {
          if (rsq < r_li_len[l]) {
            pos = l;
            break;
          }
        }
        // save index and length if needed
        if (pos != -1) {
          for (int l = nneigh-1; l > pos; --l) {
            i_li[l] = i_li[l-1];
            r_li_len[l] = r_li_len[l-1];
          }
          i_li[pos] = j;
          r_li_len[pos] = rsq;
        }
      }
    }

    for (int l=0; l<nneigh; ++l) {
      if (i_li[l] == -1) {
        char str[128];
        sprintf(str, "Nearest neighbors could not be determined for atom %i", i);
        error->one(FLERR, str);
      }
    }

    // for nneigh > 3, need to sort the neighbors in pairs
    if (nneigh>3) {
      bool check[MAXNEIGH] = {};
      int new_li[MAXNEIGH] = {};
      int ja=0;
      for (int count=0; count < nneigh; ++count) {
        new_li[count] = i_li[ja];
        check[ja] = true;
        int mi = -1;
        double md = -1;
        for (int jb=0; jb<nneigh; ++jb) {
          if (check[jb]) continue;
          double r[3] =
            {x[i_li[ja]][0] - x[i_li[jb]][0],
             x[i_li[ja]][1] - x[i_li[jb]][1],
             x[i_li[ja]][2] - x[i_li[jb]][2]};
          double rsq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
          if ((md<0) || (rsq<md)) {
            md = rsq;
            mi = jb;
          }
        }
        ja = mi;
      }
      for (int l=0; l<nneigh; ++l) {
          i_li[l] = new_li[l];
      }
    }

    double r_li[MAXNEIGH][3];
    for (int l=0; l < nneigh; ++l) {
      r_li[l][0] = xi[0] - x[ i_li[l] ][0];
      r_li[l][1] = xi[1] - x[ i_li[l] ][1];
      r_li[l][2] = xi[2] - x[ i_li[l] ][2];
    }

    /* n_lm
     *
     * pairwise normal for nearest neighbors
     */
    double n_lm[MAXNEIGH][3], n_lm_norm[MAXNEIGH][3], n_lm_ilen[MAXNEIGH];
    for (int l=0; l < nneigh; ++l) {
      int m = (l+1) < nneigh ? l+1 : 0;
      vec_cross(r_li[l], r_li[m], n_lm[l]);
      n_lm_ilen[l] = 1.0 / sqrt(vec_dot(n_lm[l], n_lm[l]));
      n_lm_norm[l][0] = n_lm_ilen[l] * n_lm[l][0];
      n_lm_norm[l][1] = n_lm_ilen[l] * n_lm[l][1];
      n_lm_norm[l][2] = n_lm_ilen[l] * n_lm[l][2];
    }

    /* n_i
     *
     * local normals
     * normalized average over normalized pairwise normals
     */
    double n_i[3] = {0, 0, 0};
    for (int l=0; l < nneigh; ++l) {
      n_i[0] += n_lm_norm[l][0];
      n_i[1] += n_lm_norm[l][1];
      n_i[2] += n_lm_norm[l][2];
    }
    double n_i_ilen = 1.0 / sqrt(vec_dot(n_i, n_i));
    double n_i_norm[3] = {n_i[0]*n_i_ilen,
                          n_i[1]*n_i_ilen,
                          n_i[2]*n_i_ilen};

    /* non-helper-atom of same type:
     *
     * energy for i and j
     * forces on i, j and neighbors
     */
    for (int jj = 0; jj < jnum; ++jj) {
      const int j      = jlist[jj] & NEIGHMASK;
      const int jtype  = type[j];
      if ((tmap[jtype] > 0) && (jtype != itype)) {
        const double r_ij[3] = {x[j][0] - xi[0],
                                x[j][1] - xi[1],
                                x[j][2] - xi[2]};
        const double rsq = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
        const RDIPParam &param = pmap[itype][jtype];

        if (rsq < cut_sq) {
          // Helper variables in j
          const double delta_isq = 1.0/(param.delta * param.delta);
          const double z0_6 = param.z0*param.z0*param.z0*param.z0*param.z0*param.z0;
          const double r_ij_len   = sqrt(rsq);
          const double exp_V      = exp(-param.lambda * (r_ij_len - param.z0));
          const double ni_dot_rij = vec_dot(n_i_norm, r_ij);
          const double rho_ij_sq  = rsq - ni_dot_rij*ni_dot_rij;
          const double rho_ij     = sqrt(rho_ij_sq);
          const double exp_f      = exp(rho_ij_sq * -delta_isq);
          const double f_rho_ij   = exp_f *
            (param.C0 + param.C2 * rho_ij_sq * delta_isq +
             param.C4 * rho_ij_sq * rho_ij_sq * delta_isq * delta_isq);
          const double df_rho_ij  = -2.0 * delta_isq * f_rho_ij
            + exp_f * (2 * param.C2 * delta_isq +
            4 * param.C4 * rho_ij_sq * delta_isq * delta_isq);
          const double rn_ij[3] =
            {(r_ij[0] - n_i_norm[0] * ni_dot_rij) * n_i_ilen,
             (r_ij[1] - n_i_norm[1] * ni_dot_rij) * n_i_ilen,
             (r_ij[2] - n_i_norm[2] * ni_dot_rij) * n_i_ilen};

          // Helper variables in l
          double rn_ijlm[MAXNEIGH][3];
          for (int l=0; l < nneigh; ++l) {
            double nlm_dot_rnij = vec_dot(n_lm_norm[l], rn_ij);
            rn_ijlm[l][0] = (rn_ij[0] - n_lm_norm[l][0] * nlm_dot_rnij) * n_lm_ilen[l];
            rn_ijlm[l][1] = (rn_ij[1] - n_lm_norm[l][1] * nlm_dot_rnij) * n_lm_ilen[l];
            rn_ijlm[l][2] = (rn_ij[2] - n_lm_norm[l][2] * nlm_dot_rnij) * n_lm_ilen[l];
          }

          // FORCES
          /* i and j:
           *
           * tmp*r:
           *  0.5 * vdW-term                (added for ij and ji)
           *  0.5 * dExp*C                  (added for ij and ji)
           *  1.0 * dExp*f(rho_ij)          (added only in ij)
           *
           * 1.0 * Exp*(r-n*dot)*df(rho_ij) (added only in ij)
           *
           */
          double tmp = 3 * param.A * z0_6 / (rsq*rsq*rsq*rsq)
            - param.lambda * exp_V * (0.5 * param.C + f_rho_ij) / r_ij_len;
          double force_j[3] =
            {tmp * r_ij[0] + exp_V * df_rho_ij * (r_ij[0] - n_i_norm[0] * ni_dot_rij),
             tmp * r_ij[1] + exp_V * df_rho_ij * (r_ij[1] - n_i_norm[1] * ni_dot_rij),
             tmp * r_ij[2] + exp_V * df_rho_ij * (r_ij[2] - n_i_norm[2] * ni_dot_rij)};

          /* i and l:
           *
           * 1.0* (added only in ij)
           */
          double force_l[MAXNEIGH][3];
          for (int l=0; l < nneigh; ++l) {
            int m = ((l+1) < nneigh) ? (l+1) : 0;
            int n = (l == 0) ? (nneigh-1) : (l-1);
            double left[3], right[3];
            vec_cross(r_li[m], rn_ijlm[l], left);
            vec_cross(r_li[n], rn_ijlm[n], right);
            force_l[l][0] = exp_V * df_rho_ij * ni_dot_rij * (left[0] - right[0]);
            force_l[l][1] = exp_V * df_rho_ij * ni_dot_rij * (left[1] - right[1]);
            force_l[l][2] = exp_V * df_rho_ij * ni_dot_rij * (left[2] - right[2]);
          }

          for (int k=0; k < 3 ; ++k) {
            f[i][k] += force_j[k];
            f[j][k] -= force_j[k];
            for (int l=0; l<nneigh; ++l) {
              f[i_li[l]][k] -= force_l[l][k];
              f[i][k] += force_l[l][k];
            }
          }

          // ENERGY
          double evdwl;
          if (eflag) {
            /* C-term and A-term will be duplicated
             * in iterations ij and ji -> *0.5
             *
             * f_rho_ij will occur only here -> *1
             *
             * tally will halve all terms,
             * but they will be added on i and j,
             * countering the effect
             */
            evdwl =
                exp_V * (0.5 * param.C + f_rho_ij) -
                0.5 * param.A * z0_6 / (rsq * rsq * rsq);
          } else {
            evdwl = 0;
          }

          // Add energies and virial terms:
          if (evflag) {
            ev_tally_xyz(i, j, nlocal, newton_pair,
                         evdwl, 0,
                         force_j[0], force_j[1], force_j[2],
                         r_ij[0], r_ij[1], r_ij[2]);
          }
          // Add virial terms for i-l interaction only when needed
          if (vflag) {
            for (int l=0; l < nneigh; ++l) {
              ev_tally_xyz(i, i_li[l], nlocal, newton_pair,
                      0, 0,
                      force_l[l][0], force_l[l][1], force_l[l][2],
                      -r_li[l][0], -r_li[l][1], -r_li[l][2]);
            }
          }
        }
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairRDIP::allocate()
{
  allocated = 1;
  int n   = atom->ntypes;

  // setup setflag
  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++) {
    for (int j = i; j <= n; j++) {
      setflag[i][j] = 0;
    }
  }
  // needs to be created here
  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  // create mapping for the atom types
  memory->create(tmap,n+1,"pair:tmap");
  memory->create(ntmap,n+1,"pair:ntmap");
  memory->create(pmap,n+1,n+1,"pair:pmap");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairRDIP::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cutoff = force->numeric(FLERR,arg[0]);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairRDIP::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
      error->all(FLERR,"Incorrect args for pair coefficients");

  const int ntyp = atom->ntypes;
  if (narg != 3 + ntyp)
      error->all(FLERR,"Incorrect number of args for pair coefficients");

  // mapping for atom type number and C-atoms
  std::set<std::string> commset;
  std::string commvec[ntyp+1];
  for (int i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      tmap[i-2] = 0;
    } else if (strcmp(arg[i],"EDGE") == 0) {
      tmap[i-2] = -1;
    } else {
      tmap[i-2] = i-2;
      commvec[i-2] = arg[i];
      commset.insert(arg[i]);
    }
  }

  // read potential file
  if (comm->me == 0) {
    char s[MAXLINE];
    std::map<std::string,RDIPType> type_file_map;
    std::map<std::set<std::string>,RDIPParam> pair_file_map;
    std::set<std::string> fileset;
    if (comm->me == 0) {
      FILE *fp = fopen(arg[2], "r");
      if (fp == NULL) {
        char str[128];
        sprintf(str,"Cannot open RDIP potential file %s", arg[2]);
        error->one(FLERR, str);
      }
      while (1) {
        fgets(s, MAXLINE, fp);
        if (feof(fp)) break;
        if (s != NULL) {
          if (s[0] == '#') continue;
          int test;
          RDIPType t;
          char itype[12], jtype[12];
          //try to parse atom type
          test = sscanf(s,"%s %d %lg", &itype, &t.nneigh, &t.cut_neigh);
          if (test == 3) {
            type_file_map[itype] = t;
            fileset.insert(itype);
            continue;
          }
          //else, parse pair parameters
          test = sscanf(s,"%s %s", &itype, &jtype);
          if (test < 2) {
            char str[128];
            sprintf(str, "Wrong number of coefficients in line \"%s\"", s);
            error->one(FLERR, str);
          }
          std::set<std::string> tempset;
          tempset.insert(itype);
          tempset.insert(jtype);
          RDIPParam &p = pair_file_map[tempset];
          //read and check parameters
          test = sscanf(s,"%*s %*s %lg %lg %lg %lg %lg %lg %lg %lg",
                        &p.A, &p.C, &p.C0, &p.C2, &p.C4,
                        &p.delta, &p.z0, &p.lambda);
          if (test != 8) {
            char str[128];
            sprintf(str, "Wrong number of coefficients in line \"%s\"", s);
            error->one(FLERR, str);
          }
        }
      }
      fclose(fp);
    }

    if (!std::includes(fileset.begin(),fileset.end(),commset.begin(),commset.end()))
        error->one(FLERR, "RDIP potential file is missing needed parameters");

    // set setflag and parameters for all atom types
    int count = 0;
    for (int i = 1; i <= ntyp; ++i) {
      if(tmap[i] > 0) ntmap[i] = type_file_map[commvec[i]];
      for (int j = i; j <= ntyp; j++) {
        if (!tmap[i] || !tmap[j]){
          setflag[i][j] = 0;
        } else {
          setflag[i][j] = 1;
          count++;
          if ((tmap[i]<0) || (tmap[j]<0)) continue;
          //insert correct parameters into pmap
          std::set<std::string> tmpset;
          tmpset.insert(commvec[i]);
          tmpset.insert(commvec[j]);
          std::map<std::set<std::string>,RDIPParam>::iterator
            search = pair_file_map.find(tmpset);
          if(search != pair_file_map.end()){
            pmap[i][j] = search->second;
            pmap[j][i] = search->second;
          } else {
            char str[128];
            sprintf(str, "Parameter set missing for pair %s %s", commvec[i].c_str(), commvec[j].c_str());
            error->one(FLERR,str);
          }
        }
      }
    }
    if (count == 0) error->one(FLERR,"All interaction is turned off");
  }
  MPI_Bcast(tmap,ntyp+1,MPI_INT,0,world);
  MPI_Bcast(ntmap,sizeof(RDIPType)*(ntyp+1),MPI_BYTE,0,world);
  for(int i = 1; i <= ntyp; ++i) {
    MPI_Bcast(pmap[i],sizeof(RDIPParam)*(ntyp+1),MPI_BYTE,0,world);
    MPI_Bcast(setflag[i],ntyp+1,MPI_INT,0,world);
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
   ------------------------------------------------------------------------- */

void PairRDIP::init_style()
{
  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->half  = 0;
  neighbor->requests[irequest]->full  = 1;

  if (force->newton_pair == 0)
  error->all(FLERR,"Pair style RDIP requires newton pair on");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairRDIP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  return cutoff;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairRDIP::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  for (int i=1; i <= atom->ntypes; i++) {
    fwrite(&tmap[i], sizeof(int), 1, fp);
    fwrite(&ntmap[i],sizeof(RDIPType),1,fp);
    for (int j = i; j <= atom->ntypes; j++) {
      fwrite(&pmap[i][j], sizeof(RDIPParam), 1, fp);
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairRDIP::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  const int me = comm->me;
  const int ntyp = atom->ntypes;
  for (int i=1; i <= ntyp; i++) {
    if (me==0) {
      fread(&tmap[i],sizeof(int),1,fp);
      fread(&ntmap[i],sizeof(RDIPType),1,fp);
    }
    for (int j=i; j <= ntyp; j++) {
      if (me==0) {
        fread(&pmap[i][j],sizeof(RDIPParam),1,fp);
        fread(&setflag[i][j],sizeof(int),1,fp);
      }
    }
  }
  MPI_Bcast(tmap, ntyp+1, MPI_INT,0,world);
  MPI_Bcast(ntmap, sizeof(RDIPType)*(ntyp+1),MPI_BYTE,0,world);
  for (int i=1; i<=ntyp; ++i) {
    MPI_Bcast(pmap[i],sizeof(RDIPParam)*(ntyp+1),MPI_BYTE,0,world);
    MPI_Bcast(setflag[i],ntyp+1,MPI_INT,0,world);
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairRDIP::write_restart_settings(FILE *fp)
{
  fwrite(&cutoff, sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairRDIP::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cutoff,     sizeof(double),1,fp);
  }
  MPI_Bcast(&cutoff,1,MPI_DOUBLE,0,world);
}
