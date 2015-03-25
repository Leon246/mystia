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
   Contributing author: Xiaoyin Ji
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_ctip.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "group.h"
#include "update.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 1024
#define DELTA 4
#define PGDELTA 1
#define MAXNEIGH 24

/* ---------------------------------------------------------------------- */

PairCtip::PairCtip(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 1;
  restartinfo = 0;
  one_coeff = 1;

  nmax = 0;
  NCo = NULL;
  bbij = NULL;

  nelements = 0;
  elements = NULL;
  nparams = 0;
  maxparam = 0;
  params = NULL;
  //elem2param = NULL;

  intype = NULL;
  fafb = NULL;
  dfafb = NULL;
  ddfafb = NULL;
  afb = NULL;
  dafb = NULL;
  ddafb = NULL;
  phin = NULL;
  dphin = NULL;
  erpaw = NULL;

  maxpage = 0;
  pages = NULL;

  // set comm size needed by this Pair

  comm_forward = 1;
  comm_reverse = 1;
  
  alf = 0.2;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairCtip::~PairCtip()
{
  memory->destroy(NCo);

  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];

  delete [] elements;
  memory->sfree(params);
  //memory->destroy(elem2param);

  memory->destroy(intype);
  memory->destroy(fafb);
  memory->destroy(dfafb);
  memory->destroy(ddfafb);
  memory->destroy(afb);
  memory->destroy(dafb);
  memory->destroy(ddafb);
  memory->destroy(phin);
  memory->destroy(dphin);
  memory->destroy(erpaw);
  memory->destroy(bbij);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] esm;
  }
}

/* ---------------------------------------------------------------------- */

void PairCtip::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum;
  int itag,jtag,itype,jtype,iparam_ij;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int mr1,mr2,mr3;
  int inty;
  double iq,jq;
  double iz,jz;
  double yaself;
  double potal,fac11,fac11e;
  double vionij,fvionij,sr1,sr2,sr3,Eov,Fov;
  int nj;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = vflag_atom = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  yaself = vionij = fvionij = Eov = Fov = 0.0;

  // self energy correction term: potal

  potal_calc(potal,fac11,fac11e);

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    itype = type[i];
    if (itype > nelements) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    iq = q[i];
    iz = params[itype].Z;
    NCo[i] = 0;
    nj = 0;
    
    // self energy, only on i atom

    yaself = self(itype,iq,potal);

    if (evflag) ev_tally(i,i,nlocal,0,yaself,0.0,0.0,0.0,0.0,0.0);

    // two-body interactions

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtag = tag[j];
      jtype = type[j];
      if (jtype > nelements) continue;

      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < x[i][2]) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

      // Qj calculates 2-body Coulombic

      jq = q[j];
      jz = params[jtype].Z;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      // long range q-dependent

      //if (rsq > params[iparam_ij].lcutsq) continue;
      if (rsq > 144.0) continue;
      
      // polynomial three-point interpolation

      tri_point(rsq, mr1, mr2, mr3, sr1, sr2, sr3);

      // 1/r energy and forces

      direct(itype-1,jtype-1,mr1,mr2,mr3,rsq,sr1,sr2,sr3,iq,jq,iz,jz,
             potal,fac11,fac11e,vionij,fvionij);

      // field correction to self energy

      //field(&params[iparam_ij],rsq,iq,jq,vionij,fvionij);

      // polarization field
      // sums up long range forces

      f[i][0] += delx*fvionij;
      f[i][1] += dely*fvionij;
      f[i][2] += delz*fvionij;
      f[j][0] -= delx*fvionij;
      f[j][1] -= dely*fvionij;
      f[j][2] -= delz*fvionij;

      if (evflag)
        ev_tally(i,j,nlocal,newton_pair,0.0,vionij,fvionij,delx,dely,delz);

    }
  }

  cuo_flag = 0;

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairCtip::allocate()
{
 allocated = 1;
 int n = atom->ntypes;

 memory->create(setflag,n+1,n+1,"pair:setflag");
 memory->create(cutsq,n+1,n+1,"pair:cutsq");

 esm = new double[n];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCtip::settings(int narg, char **arg)
{
  if (narg > 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairCtip::coeff(int narg, char **arg)
{
  int i,j,n;
  
  n = atom->ntypes;

  if (!allocated) allocate();

  //if (narg != 3 + atom->ntypes)
  //  error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  nelements = atoi(arg[2]);
  if (nelements > atom->ntypes)
    error->all(FLERR,"CTIP atom types exceeded original atom types");

  // built in CTIP parameters
  // J.Phys Condens Matter 17 (2005) 3619
  
  nparams = nelements;
  params = (Param *) memory->srealloc(params,(4)*sizeof(Param),
                                          "pair:params");  
  // Al
  params[3].type = 3;
  params[3].q_min = 0;
  params[3].q_max = 3;
  params[3].chi = -1.47914;
  params[3].J = 9.07222;
  params[3].xi = 0.968;
  params[3].Z = 1.07514;
  // Ni
  params[2].type = 2;
  params[2].q_min = 0;
  params[2].q_max = 2;
  params[2].chi = -1.70804;
  params[2].J = 9.10954;
  params[2].xi = 1.087;
  params[2].Z = 1.44450;
  // O
  params[1].type = 1;
  params[1].q_min = -2;
  params[1].q_max = 0;
  params[1].chi = 2.000;
  params[1].J = 14.99523;
  params[1].xi = 2.144;
  params[1].Z = 0.000;
  
  cutmax = 12;
  cutmin = 3.05*3.05+0.2;
  
  setup();

  // generate streitz-mintmire direct 1/r energy look-up table

  if (comm->me == 0 && screen) fprintf(screen,"Pair CTIP:\n");
  if (comm->me == 0 && screen)
    fprintf(screen,"  generating Coulomb integral lookup table ...\n");
  sm_table();
  if (comm->me == 0 && screen)
    fprintf(screen,"  generating Coulomb integral lookup table ... done\n");
  
  // clear setflag since coeff() called once with I,J = * *

  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (i-1 >= 0 && j-1 >= 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairCtip::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style COMB requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style COMB requires newton pair on");
  if (!atom->q_flag)
    error->all(FLERR,"Pair style COMB requires atom attribute q");

  // ptr to QEQ fix

  //for (i = 0; i < modify->nfix; i++)
  //  if (strcmp(modify->fix[i]->style,"qeq") == 0) break;
  //if (i < modify->nfix) fixqeq = (FixQEQ *) modify->fix[i];
  //else fixqeq = NULL;

  // need a full neighbor list

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  pgsize = neighbor->pgsize;
  oneatom = neighbor->oneatom;
  if (maxpage == 0) add_pages();
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCtip::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairCtip::setup()
{
 
}

/* ---------------------------------------------------------------------- */

double PairCtip::self(int itype, double qi, double selfpot)
{
  double self_tmp, cmin, cmax, qmin, qmax;
  double s1=params[itype].chi, s2=0.5*params[itype].J;

  self_tmp = 0.0;
  qmin = params[itype].q_min;
  qmax = params[itype].q_max;
 
  // omega
  // Phys Rev B 69 035402
  cmin = cmax = 2.0*20.0;

  self_tmp = qi*(s1+qi*(s2+selfpot));

  if (qi < qmin) self_tmp += cmin*(qi-qmin)*(qi-qmin);
  if (qi > qmax) self_tmp += cmax*(qi-qmax)*(qi-qmax);

  return self_tmp;
}

/* ---------------------------------------------------------------------- */

void PairCtip::sm_table()
{
  int i,j,k,m,nntypes,ncoul;
  int inty, itype, jtype;
  double r,dra,drin,rc,z,zr,zrc,ea,eb,ea3,eb3;
  double exp2er,exp2ersh,fafash,dfafash,F1,dF1,ddF1,E1,E2,E3,E4;
  double exp2ear,exp2ebr,exp2earsh,exp2ebrsh,fafbsh,dfafbsh;
  double afbsh,dafbsh;

  int n = nelements;
  int nmax = atom->nmax;

  dra  = 0.001;  // lookup table step size
  drin = 0.1;    // starting distance of 1/r
  rc = cutmax;

  nntypes = int((n+1)*n/2); // interaction types
  ncoul = int((rc-drin)/dra)+1;

  // allocate arrays

  memory->create(intype,n,n,"pair:intype");
  memory->create(fafb,ncoul,nntypes,"pair:fafb");
  memory->create(dfafb,ncoul,nntypes,"pair:dfafb");
  memory->create(ddfafb,ncoul,nntypes,"pair:ddfafb");
  memory->create(phin,ncoul,nntypes,"pair:phin");
  memory->create(dphin,ncoul,nntypes,"pair:dphin");
  memory->create(erpaw,25000,2,"pair:erpaw");
  memory->create(NCo,nmax,"pair:NCo");
  memory->create(bbij,nmax,MAXNEIGH,"pair:bbij");
  
  memory->create(afb,ncoul,n,"pair:afb");
  memory->create(dafb,ncoul,n,"pair:dafb");
  memory->create(ddafb,ncoul,n,"pair:ddafb");

  // set interaction number: 0-0=0, 1-1=1, 0-1=1-0=2

  m = 0; k = n;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (j == i) {
        intype[i][j] = m;
        m += 1;
      } else if (j != i && j > i) {
        intype[i][j] = k;
        k += 1;
      } else if (j != i && j < i) {
        intype[i][j] = intype[j][i];
      }
    }
  }

  // default arrays to zero
  for (i = 0; i < ncoul; i ++) {
    for (j = 0; j < nntypes; j ++) {
      fafb[i][j] = 0.0;
      dfafb[i][j] = 0.0;
      ddfafb[i][j] = 0.0;
      phin[i][j] = 0.0;
      dphin[i][j] = 0.0;
    }
  }

  // direct 1/r energy with Slater 1S orbital overlap

  for (i = 0; i < n; i++) {
    r = drin;
    z = params[i+1].Z;
    for (j = 0; j < ncoul; j++) {
      exp2er = exp(-2.0 * z * r);
      phin[j][i] = 1.0 - exp2er * (1.0 + 2.0 * z * r * (1.0 + z * r));
      dphin[j][i] = (4.0 * exp2er * z * z * z * r * r);
      r += dra;
    }
  }
  
  // [i|fj] = f(j)
  
  for (j = 0; j < n; j++) {
    r = drin;
    eb = params[j+1].xi;
    exp2ebrsh = exp(-2.0*eb*rc);
    afbsh = -exp2ebrsh * (eb + 1.0/rc);
    dafbsh = exp2ebrsh * (2.0*eb*eb + 2.0*eb/rc + 1.0/(rc*rc));
    for (k = 0; k < ncoul; k++) {
      exp2ebr = exp(-2.0*eb*r);
      F1 = -exp2ebr*(eb + 1.0/r);
      dF1 = exp2ebr*(2.0*eb*eb + 2.0*eb/r + 1.0/(r*r));
      ddF1 = -2.0*exp2ebr*(1.0/(r*r*r) + 
              eb/r + eb*(1.0/(r*r) + 2*eb/r + 2*eb*eb));
      afb[k][j] = F1 - afbsh - (r-rc)*dafbsh;
      dafb[k][j] = dF1 - dafbsh;
      ddafb[k][j] = ddF1;
      r += dra;
    }
  }
  
  // [fi|fj]
  
  for (i = 0; i < n; i ++) {
    for (j = 0; j < n; j ++) {
      r = drin;
      if (j == i) {
        itype = params[i+1].type-1;
        inty = intype[itype][itype];
        z = params[i+1].xi;
        zrc = z * rc;
        exp2ersh = exp(-2.0 * zrc);
        fafash = -exp2ersh * (1.0 / rc +
                              z * (11.0/8.0 + 3.0/4.0*zrc + zrc*zrc/6.0));
        dfafash = exp2ersh * (1.0/(rc*rc) + 2.0*z/rc +
                              z*z*(2.0 + 7.0/6.0*zrc + zrc*zrc/3.0));
        for (k = 0; k < ncoul; k ++) {
          zr = z * r;
          exp2er = exp(-2.0*zr);
          F1 = -exp2er * (1.0 / r +
                          z * (11.0/8.0 + 3.0/4.0*zr + zr*zr/6.0));
          dF1 = exp2er * (1.0/(r*r) + 2.0*z/r +
                          z*z*(2.0 + 7.0/6.0*zr + zr*zr/3.0));
          ddF1 = -exp2er * (2.0/(r*r*r) + 4.0*z/(r*r) -
                            z*z*z/3.0*(17.0/2.0 + 5.0*zr + 2.0*zr*zr));
          fafb[k][inty] = F1-fafash-(r-rc)*dfafash;
          dfafb[k][inty] = (dF1 - dfafash);
          ddfafb[k][inty] = ddF1;
          r += dra;
        }
      } else if (j != i) {
        itype = params[i+1].type-1;
        jtype = params[j+1].type-1;
        inty = intype[itype][jtype];
        ea = params[i+1].xi;
        ea3 = ea*ea*ea;
        eb = params[j+1].xi;
        eb3 = eb*eb*eb;
        // hidden
        // E1 = ea*ea3*ea/((ea+eb)*(ea+eb)*(ea-eb)*(ea-eb));
        // E2 = eb*eb3*eb/((ea+eb)*(ea+eb)*(eb-ea)*(eb-ea));
        E1 = ea*eb3*eb/((ea+eb)*(ea+eb)*(ea-eb)*(ea-eb));
        E2 = eb*ea3*ea/((ea+eb)*(ea+eb)*(eb-ea)*(eb-ea));
        E3 = (3.0*ea*ea*eb3*eb-eb3*eb3) /
          ((ea+eb)*(ea+eb)*(ea+eb)*(ea-eb)*(ea-eb)*(ea-eb));
        E4 = (3.0*eb*eb*ea3*ea-ea3*ea3) /
          ((ea+eb)*(ea+eb)*(ea+eb)*(eb-ea)*(eb-ea)*(eb-ea));
        exp2earsh = exp(-2.0*ea*rc);
        exp2ebrsh = exp(-2.0*eb*rc);
        fafbsh = -exp2earsh*(E1 + E3/rc)-exp2ebrsh*(E2 + E4/rc);
        dfafbsh =
          exp2earsh*(2.0*ea*(E1+E3/rc)+E3/(rc*rc)) +
          exp2ebrsh*(2.0*eb*(E2+E4/rc)+E4/(rc*rc));
        for (k = 0; k < ncoul; k ++) {
          exp2ear = exp(-2.0*ea*r);
          exp2ebr = exp(-2.0*eb*r);
          fafb[k][inty] =
            - exp2ear*(E1+E3/r) - exp2ebr*(E2+E4/r)
            - fafbsh - (r-rc) * dfafbsh;
          dfafb[k][inty] = (exp2ear*(2.0*ea*(E1+E3/r) + E3/(r*r))
                            + exp2ebr*(2.0*eb*(E2+E4/r) + E4/(r*r))- dfafbsh);
          ddfafb[k][inty] = (- exp2ear*(E3/(r*r)*(1.0/r+2.0*ea/r+2.0/(r*r))
                                        + 2.0*ea*(E1+E3/r))-
                             exp2ebr*(E4/(r*r)
                                      *(1.0/r+2.0*eb/r+2.0/(r*r)) +
                                      2.0*eb*(E2+E4/r)));
          r += dra;
        }
      }
    }
  }

  for (i = 0; i < 25000; i ++) {
    r = dra * i + drin;
    erpaw[i][0] = erfc(r*alf);
    erpaw[i][1] = exp(-r*r*alf*alf);
  }
   
}

/* ---------------------------------------------------------------------- */

void PairCtip::potal_calc(double &calc1, double &calc2, double &calc3)
{
  double rcoul,esucon;
  int m;

  rcoul = 12.0;
  //for (m = 0; m < nparams; m++)
  //  if (params[m].lcut > rcoul) rcoul = params[m].lcut;

  esucon = force->qqr2e;

  calc2 = (erfc(rcoul*alf)/rcoul/rcoul+2.0*alf/MY_PIS*
           exp(-alf*alf*rcoul*rcoul)/rcoul)*esucon/rcoul;
  calc3 = (erfc(rcoul*alf)/rcoul)*esucon;
  calc1 = -(alf/MY_PIS*esucon+calc3*0.5);
}

/* ---------------------------------------------------------------------- */

void PairCtip::tri_point(double rsq, int &mr1, int &mr2,
                         int &mr3, double &sr1, double &sr2,
                         double &sr3)
{
 double r, rin, dr, dd, rr1, rridr, rridr2;

 rin = 0.10; dr = 0.0010;
 r = sqrt(rsq);
 if (r < rin + 2.0*dr) r = rin + 2.0*dr;
 if (r > cutmax - 2.0*dr) r = cutmax - 2.0*dr;
 rridr = (r-rin)/dr;

 mr1 = int(rridr)-1;
 dd = rridr - float(mr1);
 if (dd > 0.5) mr1 += 1;
 mr2 = mr1 + 1;
 mr3 = mr2 + 1;

 rr1 = float(mr1)*dr;
 rridr = (r - rin - rr1)/dr;
 rridr2 = rridr * rridr;

 sr1 = (rridr2 - rridr) * 0.50;
 sr2 = 1.0 - rridr2;
 sr3 = (rridr2 + rridr) * 0.50;
}

/* ---------------------------------------------------------------------- */

void PairCtip::direct(int itype, int jtype, int mr1, int mr2, int mr3, double rsq,
                      double sr1, double sr2, double sr3,
                      double iq, double jq,
                      double iz, double jz,
                      double potal, double fac11, double fac11e,
                      double &pot_tmp, double &pot_d)
{
  int inty;
  double r,erfcc,fafbn1,potij,sme2,esucon;
  double r3,erfcd,dfafbn1,smf2,dvdrr,alfdpi;
  double afbn1,bfan1,dafbn1,dbfan1;
  
  inty = intype[itype][jtype];

  r = sqrt(rsq);
  r3 = r * rsq;
  alfdpi = 2.0*alf/MY_PIS;
  esucon = force->qqr2e;
  pot_tmp = 0.0;
  pot_d = 0.0;

  // 1/r energy

  erfcc = sr1*erpaw[mr1][0] + sr2*erpaw[mr2][0] + sr3*erpaw[mr3][0];
  fafbn1= sr1*fafb[mr1][inty] + sr2*fafb[mr2][inty] + sr3*fafb[mr3][inty];
  potij = (erfcc/r * esucon - fac11e);
  sme2 = potij + fafbn1 * esucon;

  afbn1 = sr1*afb[mr1][jtype] + sr2*afb[mr2][jtype] + sr3*afb[mr3][jtype];
  bfan1 = sr1*afb[mr1][itype] + sr2*afb[mr2][itype] + sr3*afb[mr3][itype];

  // Vij(rij,qi,qj)
  // kc*qi*qj*[fi|fj] + kc*qi*Zj([j|fi]-[fi|fj]) + kc*qj*Zi([i|fj]-[fi|fj])
  // + kc*Zi*Zj([fi|fj]-[i|fj]-[j|fi]+1/rij) 
  pot_tmp = 1.0 * iq * jq * sme2;
  pot_tmp += 1.0 * iq * jz * (bfan1 - fafbn1) * esucon;
  pot_tmp += 1.0 * jq * iz * (afbn1 - fafbn1) * esucon;
  //pot_tmp += 1.0 * iz * jz * (fafbn1 - afbn1 - bfan1) * esucon;
    
  // 1/r force (wrt r)

  erfcd = sr1*erpaw[mr1][1] + sr2*erpaw[mr2][1] + sr3*erpaw[mr3][1];
  dfafbn1= sr1*dfafb[mr1][inty] + sr2*dfafb[mr2][inty] + sr3*dfafb[mr3][inty];
  dvdrr = (erfcc/r3+alfdpi*erfcd/rsq)*esucon-fac11;
  smf2 = dvdrr - dfafbn1 * esucon/r;
  pot_d =  1.0 * iq * jq * smf2;

  dafbn1= sr1*dafb[mr1][jtype] + sr2*dafb[mr2][jtype] + sr3*dafb[mr3][jtype];
  dbfan1= sr1*dafb[mr1][itype] + sr2*dafb[mr2][itype] + sr3*dafb[mr3][itype];

  pot_d -= 1.0 * iq * jz * (dbfan1 - dfafbn1) * esucon/r;
  pot_d -= 1.0 * jq * iz * (dafbn1 - dfafbn1) * esucon/r;
  //pot_d -= 1.0 * iz * jz * (dfafbn1 - dafbn1 - dbfan1) * esucon/r;

}

/* ---------------------------------------------------------------------- */

double PairCtip::yasu_char(double *qf_fix, int &igroup)
{
  int i,j,ii,jj,jnum,itag,jtag;
  int itype,jtype;
  double xtmp,ytmp,ztmp;
  double rsq1,delr1[3];
  int *ilist,*jlist,*numneigh,**firstneigh;
  double iq,jq,iz,jz,qfj,fqi,fij,fqij,fqjj;
  double potal,fac11,fac11e,sr1,sr2,sr3;
  int mr1,mr2,mr3,inty,nj;


  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  int *tag = atom->tag;

  int inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  int *mask = atom->mask;
  int groupbit = group->bitmask[igroup];

  qf = qf_fix;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit)
      qf[i] = 0.0;
  }

  // communicating charge force to all nodes, first forward then reverse

  comm->forward_comm_pair(this);

  // self energy correction term: potal

  potal_calc(potal,fac11,fac11e);

  // loop over full neighbor list of my atoms

  fqi = qfj = fij = 0.0;

  for (ii = 0; ii < inum; ii ++) {
    i = ilist[ii];
    itype = type[i];
    if (itype > nelements) continue;
    itag = tag[i];
    nj = 0;
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      iq = q[i];
      iz = params[itype].Z;

      // charge force from self energy

      fqi = qfo_self(itype,iq,potal);

      // two-body interactions

      jlist = firstneigh[i];
      jnum = numneigh[i];
      
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        jtype = type[j];
        if (jtype > nelements) continue;
        jtag = tag[j];
        
        if (itag > jtag) {
          if ((itag+jtag) % 2 == 0) continue;
        } else if (itag < jtag) {
          if ((itag+jtag) % 2 == 1) continue;
        } else {
          if (x[j][2] < x[i][2]) continue;
          if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
          if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
        }
        
        jq = q[j];
        jz = params[jtype].Z;
        
        delr1[0] = x[j][0] - xtmp;
        delr1[1] = x[j][1] - ytmp;
        delr1[2] = x[j][2] - ztmp;
        rsq1 = vec3_dot(delr1,delr1);

        // long range q-dependent

        if (rsq1 > 144.0) continue;

        // polynomial three-point interpolation

        tri_point(rsq1,mr1,mr2,mr3,sr1,sr2,sr3);

        // 1/r charge forces

        qfo_direct(itype-1,jtype-1,mr1,mr2,mr3,rsq1,sr1,sr2,sr3,iq,jq,iz,jz,fac11e,fqij,fqjj);
        fqi += fqij;  qf[j] += fqjj;
      }
      qf[i] += fqi; 
    }
  }

  comm->reverse_comm_pair(this);

  // sum charge force on each node and return it

  double eneg = 0.0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit)
      eneg += qf[i];
  }
  double enegtot;
  MPI_Allreduce(&eneg,&enegtot,1,MPI_DOUBLE,MPI_SUM,world);
  return enegtot;
}

/* ---------------------------------------------------------------------- */

double PairCtip::qfo_self(int itype, double qi, double selfpot)
{
 double self_d, cmin, cmax, qmin, qmax;
 double s1=params[itype].chi, s2=0.5*params[itype].J;

 self_d = 0.0;
 qmin = params[itype].q_min;
 qmax = params[itype].q_max;
 
 // omega
 // Phys Rev B 69 035402
 cmin = cmax = 2.0*20.0;

 self_d = s1+2.0*qi*(s2+selfpot);

 if (qi < qmin) self_d += 2.0*cmin*(qi-qmin);
 if (qi > qmax) self_d += 2.0*cmax*(qi-qmax);

 return self_d;
}

/* ---------------------------------------------------------------------- */

void PairCtip::qfo_direct(int itype, int jtype, int mr1, int mr2, int mr3,
                          double rsq, double sr1, double sr2, double sr3,
                          double iq, double jq, double iz, double jz,
                          double fac11e, double &fqij, double &fqjj)
{
  int inty = intype[itype][jtype];;
  double r, erfcc, fafbn1, vm, esucon, fij;
  double afbn1, bfan1;

  r = sqrt(rsq);
  esucon=force->qqr2e;

  // 1/r force (wrt q)

  erfcc = sr1*erpaw[mr1][0]   + sr2*erpaw[mr2][0]   + sr3*erpaw[mr3][0];
  fafbn1= sr1*fafb[mr1][inty] + sr2*fafb[mr2][inty] + sr3*fafb[mr3][inty];
  vm = (erfcc/r * esucon - fac11e);
  fij = 1.0 * (vm+esucon*fafbn1);
  
  fqij = jq*fij;  fqjj = iq*fij;
  
  afbn1 = sr1*afb[mr1][jtype] + sr2*afb[mr2][jtype] + sr3*afb[mr3][jtype];
  bfan1 = sr1*afb[mr1][itype] + sr2*afb[mr2][itype] + sr3*afb[mr3][itype];

  // Vij(rij,qi,qj)
  // kc*qi*qj*[fi|fj] + kc*qi*Zj([j|fi]-[fi|fj]) + kc*qj*Zi([i|fj]-[fi|fj])
  // + kc*Zi*Zj([fi|fj]-[i|fj]-[j|fi]+1/rij) 
  fqij += 1.0 * jz * (bfan1 - fafbn1) * esucon;
  fqjj += 1.0 * iz * (afbn1 - fafbn1) * esucon;
}

/* ---------------------------------------------------------------------- */

int PairCtip::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i ++) {
    j = list[i];
    buf[m++] = qf[j];
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

void PairCtip::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n ;
  for (i = first; i < last; i++) qf[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

int PairCtip::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = qf[i];
  return 1;
}

/* ---------------------------------------------------------------------- */

void PairCtip::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    qf[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairCtip::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);
  bytes += nmax * sizeof(int);
  bytes += MAXNEIGH * nmax * sizeof(int);
  return bytes;
}

/* ---------------------------------------------------------------------- */

void PairCtip::add_pages(int howmany)
{
  int toppage = maxpage;
  maxpage += howmany*PGDELTA;

  pages = (int **)
    memory->srealloc(pages,maxpage*sizeof(int *),"pair:pages");
  for (int i = toppage; i < maxpage; i++)
    memory->create(pages[i],pgsize,"pair:pages[i]");
}
