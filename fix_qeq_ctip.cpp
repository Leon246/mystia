
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
   Modified from fix_qeq_comb.cpp
   Xiaoyin Ji
-------------------------------------------------------------------------*/

#include "lmptype.h"
#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "pair_ctip.h"
#include "fix_qeq_ctip.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "group.h"
#include "respa.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include "iostream"
using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixQEQCtip::FixQEQCtip(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix qeq/ctip command");

  peratom_flag = 1;
  size_peratom_cols = 0;
  peratom_freq = 1;
  screentype = 0;

  nevery = force->inumeric(arg[3]);
  precision = force->numeric(arg[4]);
  screentype = force->inumeric(arg[5]);

  if (nevery <= 0 || precision <= 0.0)
    error->all(FLERR,"Illegal fix qeq/ctip command");

  MPI_Comm_rank(world,&me);

  // optional args

  fp = NULL;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qeq/ctip command");
      if (me == 0) {
        fp = fopen(arg[iarg+1],"w");
        if (fp == NULL) {
          char str[128];
          sprintf(str,"Cannot open fix qeq/ctip file %s",arg[iarg+1]);
          error->one(FLERR,str);
        }
      }
      iarg += 2;
    } else error->all(FLERR,"Illegal fix qeq/ctip command");
  }

  nmax = atom->nmax;
  memory->create(qf,nmax,"qeq:qf");
  memory->create(q1,nmax,"qeq:q1");
  memory->create(q2,nmax,"qeq:q2");
  vector_atom = qf;

  // zero the vector since dump may access it on timestep 0
  // zero the vector since a variable may access it before first run

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) qf[i] = 0.0;

  ctip = NULL;

  comm_forward = 1;
}

/* ---------------------------------------------------------------------- */

FixQEQCtip::~FixQEQCtip()
{
  if (me == 0 && fp) fclose(fp);
  memory->destroy(qf);
  memory->destroy(q1);
  memory->destroy(q2);
}

/* ---------------------------------------------------------------------- */

int FixQEQCtip::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixQEQCtip::init()
{
  if (!atom->q_flag)
    error->all(FLERR,"Fix qeq/ctip requires atom attribute q");

  ctip = (PairCtip *) force->pair_match("ctip",1);
  if (ctip == NULL)
    error->all(FLERR,"Must use pair_style ctip with fix qeq/ctip");

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  ngroup = group->count(igroup);
  if (ngroup == 0) error->all(FLERR,"Fix qeq/ctip group has no atoms");

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void FixQEQCtip::setup(int vflag)
{
  firstflag = 1;
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
  firstflag = 0;
}

/* ---------------------------------------------------------------------- */

void FixQEQCtip::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEQCtip::post_force(int vflag)
{
  int i,ii,iloop,loopmax,inum,*ilist;
  double heatpq,qmass,dtq,dtq2;
  double enegchkall,enegmaxall;

  if (update->ntimestep % nevery) return;

  // reallocate work arrays if necessary
  // qf = charge force
  // q1 = charge displacement
  // q2 = tmp storage of charge force for next iteration

  if (atom->nmax > nmax) {
    memory->destroy(qf);
    memory->destroy(q1);
    memory->destroy(q2);
    nmax = atom->nmax;
    memory->create(qf,nmax,"qeq:qf");
    memory->create(q1,nmax,"qeq:q1");
    memory->create(q2,nmax,"qeq:q2");
    vector_atom = qf;
  }

  // more loops for first-time charge equilibrium

  iloop = 0;
  if (firstflag) loopmax = 500;
  else loopmax = 200;

  // charge-equilibration loop

  if (me == 0 && fp)
    fprintf(fp,"Charge equilibration on step " BIGINT_FORMAT "\n",
            update->ntimestep);

  heatpq = 0.05;
  qmass  = 0.016;
  dtq    = 0.01;
  dtq2   = 0.5*dtq*dtq/qmass;

  double enegchk = 0.0;
  double enegtot = 0.0;
  double enegmax = 0.0;

  double *q = atom->q;
  int *mask = atom->mask;
  int *type = atom->type;

  inum = ctip->list->inum;
  ilist = ctip->list->ilist;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    q1[i] = q2[i] = qf[i] = 0.0;
  }
  
  // update natoms in group, dynamic deposition
  ngroup = group->count(igroup);

  for (iloop = 0; iloop < loopmax; iloop ++ ) {
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if ((mask[i] & groupbit) && type[i] != screentype) {
        q1[i] += qf[i]*dtq2 - heatpq*q1[i];
        q[i]  += q1[i];
      }
    }
    comm->forward_comm_fix(this);

    if(ctip) enegtot = ctip->yasu_char(qf,igroup);
    enegtot /= ngroup;
    enegchk = enegmax = 0.0;

    for (ii = 0; ii < inum ; ii++) {
      i = ilist[ii];
      if ((mask[i] & groupbit) && type[i] != screentype) {
        q2[i] = enegtot-qf[i];
        enegmax = MAX(enegmax,fabs(q2[i]));
        enegchk += fabs(q2[i]);
        qf[i] = q2[i];
      }
    }

    MPI_Allreduce(&enegchk,&enegchkall,1,MPI_DOUBLE,MPI_SUM,world);
    enegchk = enegchkall/ngroup;
    MPI_Allreduce(&enegmax,&enegmaxall,1,MPI_DOUBLE,MPI_MAX,world);
    enegmax = enegmaxall;

    if (enegchk <= precision && enegmax <= 100.0*precision) break;

    if (me == 0 && fp)
      fprintf(fp,"  iteration: %d, enegtot %.6g, "
              "enegmax %.6g, fq deviation: %.6g\n",
              iloop,enegtot,enegmax,enegchk);

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if ((mask[i] & groupbit) && type[i] != screentype)
        q1[i] += qf[i]*dtq2 - heatpq*q1[i];
    }
  }

  if (me == 0 && fp) {
    if (iloop == loopmax)
      fprintf(fp,"Charges did not converge in %d iterations\n",iloop);
    else
      fprintf(fp,"Charges converged in %d iterations to %.10f tolerance\n",
              iloop,enegchk);
  }
}

/* ---------------------------------------------------------------------- */

void FixQEQCtip::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixQEQCtip::memory_usage()
{
  double bytes = atom->nmax*3 * sizeof(double);
  return bytes;
}
/* ---------------------------------------------------------------------- */

int FixQEQCtip::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i ++) {
    j = list[i];
    buf[m++] = atom->q[j];
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

void FixQEQCtip::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n ;
  for (i = first; i < last; i++) atom->q[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */
