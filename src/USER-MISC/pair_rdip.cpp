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
      Florian Wullschlaeger, Sebastian Gsaenger 2017
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

#include "ctype.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define PGDELTA 1

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

    param.nneigh = 3;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairRDIP::~PairRDIP()
{
    if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(map);
  }
}

/* ---------------------------------------------------------------------- */
void PairRDIP::compute(int eflag, int vflag)
{
    if (eflag || vflag) ev_setup(eflag,vflag);
    else evflag = vflag_fdotr = 0;

    const int newton_pair = force->newton_pair; // always true
    const int nlocal = atom->nlocal; // not used
    const int *type  = atom->type;
    const double * const*x = atom->x;
    double **f       = atom->f;

    const int inum         = list->inum;
    const int *ilist       = list->ilist;
    const int *numneigh    = list->numneigh;
    int **firstneigh       = list->firstneigh;

    const double maximum = param.z0 * param.z0 * 2.25;
    const double cut_sq = cutoff * cutoff;
    const double delta_isq = 1.0/(param.delta * param.delta);
    const double z0_6 = param.z0*param.z0*param.z0*param.z0*param.z0*param.z0;

    for (int ii=0; ii < inum; ++ii) {
        const int i      = ilist[ii];
        const int itype  = type[i];
        if (map[itype] <= 0) continue; // skip helper and ignored atoms

        const double *xi = x[i];
        const int *jlist = firstneigh[i];
        const int jnum   = numneigh[i];

        /* r_li:
         *
         * distances between reference atom i
         * and nearest neighbors l
         */
        double r_li_len[param.nneigh]; //squared, only needed for sorting
        int    i_li[param.nneigh];
        for (int l = 0; l < param.nneigh; ++l) {
            r_li_len[l] = maximum+l;
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
            if (map[jtype] == -1 || (jtype == itype)) {
                double r_li[3] =
                    {xi[0] - x[j][0],
                     xi[1] - x[j][1],
                     xi[2] - x[j][2]};
                const double rsq = r_li[0]*r_li[0] + r_li[1]*r_li[1] + r_li[2]*r_li[2];

                int pos = -1;
                // check whether l is one of nearest neighbors
                for (int l=0; l<param.nneigh; ++l) {
                    if (rsq < r_li_len[l]) {
                        pos = l;
                        break;
                    }
                }
                // save index and length if needed
                if (pos != -1) {
                    for (int l = param.nneigh-1; l > pos; --l) {
                        i_li[l] = i_li[l-1];
                        r_li_len[l] = r_li_len[l-1];
                    }
                    i_li[pos] = j;
                    r_li_len[pos] = rsq;
                }
            }
        }

        for (int l=0; l<param.nneigh; ++l) {
            if (i_li[l] == -1) {
                char str[128];
                sprintf(str, "Nearest neighbors could not be determined for atom %i", i);
                error->one(FLERR, str);
            }
        }

        double r_li[param.nneigh][3];
        for (int l=0; l < param.nneigh; ++l) {
            r_li[l][0] = xi[0] - x[ i_li[l] ][0];
            r_li[l][1] = xi[1] - x[ i_li[l] ][1];
            r_li[l][2] = xi[2] - x[ i_li[l] ][2];
        }

        /* n_lm
         *
         * pairwise normal for nearest neighbors
         */
        double n_lm[param.nneigh][3], n_lm_norm[param.nneigh][3], n_lm_ilen[param.nneigh];
        for (int l=0; l < param.nneigh; ++l) {
            int m = (l+1) < param.nneigh ? l+1 : 0;
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
        for (int l=0; l < param.nneigh; ++l) {
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
            if ((map[jtype] > 0) && (jtype != itype)) {
                const double r_ij[3] = {x[j][0] - xi[0],
                                        x[j][1] - xi[1],
                                        x[j][2] - xi[2]};
                const double rsq = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];

                if (rsq < cut_sq) {
                    // Helper variables in j
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
                    double rn_ijlm[param.nneigh][3];
                    for (int l=0; l < param.nneigh; ++l) {
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
                    double force_l[param.nneigh][3];
                    for (int l=0; l < param.nneigh; ++l) {
                        int m = ((l+1) < param.nneigh) ? (l+1) : 0;
                        int n = (l == 0) ? (param.nneigh-1) : (l-1);
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
                        for (int l=0; l<param.nneigh; ++l) {
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
                        for (int l=0; l < param.nneigh; ++l) {
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
    int nat = atom->natoms;

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
    memory->create(map,n+1,"pair:map");
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
   read RDIP potential file
------------------------------------------------------------------------- */

void PairRDIP::read_file(char *filename)
{
  char s[MAXLINE];

  // read file on proc 0

  if (comm->me == 0) {
    FILE *fp = fopen(filename,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open RDIP potential file %s",filename);
      error->one(FLERR,str);
    }

    // skip initial comment lines

    while (1) {
      fgets(s,MAXLINE,fp);
      if (s[0] != '#') break;
    }

    // read parameters
    sscanf(s,"%*s %*s %d %lg %lg %lg %lg %lg %lg %lg %lg",
            &param.nneigh, &param.A, &param.C, &param.C0, &param.C2,
            &param.C4, &param.delta, &param.z0, &param.lambda);
  }

  MPI_Bcast(&param,sizeof(RDIPParam),MPI_BYTE,0,world);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairRDIP::coeff(int narg, char **arg)
{
    if (!allocated) allocate();

    if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
        error->all(FLERR,"Incorrect args for pair coefficients");

    if (narg != 3 + atom->ntypes)
        error->all(FLERR,"Incorrect number of args for pair coefficients");

    // mapping for atom type number and C-atoms
    for (int i = 3; i < narg; i++) {
        if (strcmp(arg[i],"NULL") == 0) {
            map[i-2] = 0;
        } else if (strcmp(arg[i],"EDGE") == 0) {
            map[i-2] = -1;
        } else if (strcmp(arg[i],"C") == 0) {
            map[i-2] = i-2;
        } else error->all(FLERR,"Incorrect args for pair coefficients");
    }

    // read potential file
    read_file(arg[2]);

    // set setflag for all atoms
    int count = 0;
    for (int i = 1; i <= atom->ntypes; i++) {
        for (int j = i; j <= atom->ntypes; j++) {
            if (!map[i]||!map[j]){
                setflag[i][j] = 0;
            } else {
                setflag[i][j] = 1;
                count++;
            }
        }
    }

    if (count == 0) error->all(FLERR,"All interaction is turned off");
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
        fwrite(&map[i], sizeof(int), 1, fp);
        for (int j = i; j <= atom->ntypes; j++) {
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

    int me = comm->me;
    for (int i=1; i <= atom->ntypes; i++) {
        if (me==0) fread(&map[i],sizeof(int),1,fp);
        MPI_Bcast(&map[i],1,MPI_INT,0,world);
        for (int j=i; j <= atom->ntypes; j++) {
            if (me==0) {
                fread(&setflag[i][j],sizeof(int),1,fp);
            }
            MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
        }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairRDIP::write_restart_settings(FILE *fp)
{
    fwrite(&cutoff,     sizeof(double),1,fp);
    fwrite(&param, sizeof(RDIPParam),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairRDIP::read_restart_settings(FILE *fp)
{
    if (comm->me == 0) {
       fread(&cutoff,     sizeof(double),1,fp);
       fread(&param,  sizeof(RDIPParam),1,fp);
    }
     MPI_Bcast(&cutoff,1,MPI_DOUBLE,0,world);
     MPI_Bcast(&param,sizeof(RDIPParam),MPI_BYTE,0,world);
}
