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
   Modified from FixGCMC
   Contributing author: Xiaoyin Ji (NCSU)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_mcswap.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "comm.h"
#include "group.h"
#include "domain.h"
#include "random_park.h"
#include "force.h"
#include "pair.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pointers.h"
#include "array.h"
#include <iostream>

using namespace std;

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define MAX_NLOOP 4096

/* ---------------------------------------------------------------------- 
 Usage
 fix fix-ID group-ID mcswap type1 type2 N-every N-MCMD seed temp displace
 swap N-swap nmaxattempt N-maxattempt
 ------------------------------------------------------------------------ */

FixMCSWAP::FixMCSWAP(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 8) error->all(FLERR,"Illegal fix MCSWAP command");

  rank = comm->me;
  size = comm->nprocs;
  cutoff = neighbor->cutneighmax;

  vector_flag = 1;
  size_vector = 6;
  global_freq = 1;
  extvector = 0;
  restart_global = 1;
  time_depend = 1;
  itype1 = 1;
  itype2 = 2;
  nevery = 0;
  nswap = 0;
  npostswap = 0;

  if (domain->triclinic == 0) {
    mcxlo = domain->boxlo[0];
    mcxhi = domain->boxhi[0];
  } else {
    mcxlo = domain->boxlo_bound[0];
    mcxhi = domain->boxhi_bound[0];
  }

  nevery = atoi(arg[3]);
  seed = atoi(arg[4]);
  reservoir_temperature = atof(arg[5]);
  displace = atof(arg[6]);
  cutoff = atof(arg[7]);

  if (seed <= 0) error->all(FLERR,"Illegal fix MCSWAP command");
  if (reservoir_temperature < 0.0) error->all(FLERR,"Illegal fix MCSWAP command");
  if (displace < 0.0) error->all(FLERR,"Illegal fix MCSWAP command");
  if (reservoir_temperature == 0) {
	beta = 1.0e30;
  } else {
    beta = 1.0/(force->boltz*reservoir_temperature);
  }

  // read options from end of input line
  
  int p = 8;
  options(narg-p,&arg[p]);

  // random number generator, same for all procs

  random_equal = new RanPark(lmp,seed);

  // random number generator, not the same for all procs

  random_unequal = new RanPark(lmp,seed);
  
  // compute the number of MC cycles that occur nevery timesteps
  
  ncycles = nmcmoves;
  
  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;

  // zero out counters
  
  nmove_attempts = 0.0;
  nmove_successes = 0.0;
  nswap_attempts = 0.0;
  nswap_successes = 0.0;

  nmax = 0;
  local_gas_list = NULL;
  local_type1_list = NULL;
  local_type2_list = NULL;
  scn = 0;
}

/* ---------------------------------------------------------------------- */

FixMCSWAP::~FixMCSWAP()
{
  delete random_equal;
  delete random_unequal;
  memory->sfree(local_gas_list);
  memory->sfree(local_type1_list);
  memory->sfree(local_type2_list);

  delete(ngas_local_all);
  delete(ntype1_local_all);
  delete(ntype2_local_all);
  delete(balance_all);
  delete(balance_type1);
  delete(balance_type2);
  delete(sc1tag);
  delete(sc2tag);
  delete(tag1);
  delete(tag2);
  delete(dtag1);
  delete(dtag2);
  delete(dtag1_all);
  delete(dtag2_all);
  delete(taglist);
  delete(sc1x);
  delete(sc1y);
  delete(sc1z);
  delete(sc2x);
  delete(sc2y);
  delete(sc2z);
  delete(svc1x);
  delete(svc1y);
  delete(svc1z);
  delete(svc2x);
  delete(svc2y);
  delete(svc2z);
  delete(velx1);
  delete(vely1);
  delete(velz1);
  delete(velx2);
  delete(vely2);
  delete(velz2);
  delete(coordx1);
  delete(coordy1);
  delete(coordz1);
  delete(coordx2);
  delete(coordy2);
  delete(coordz2);
  if (rank == 0) {
    delete(denergy1_all);
    delete(denergy2_all);
    delete(tag1list);
    delete(tag2list);
  }
}

/* ---------------------------------------------------------------------- */

int FixMCSWAP::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMCSWAP::init()
{
  // check that no deletable atoms are in atom->firstgroup
  // deleting such an atom would not leave firstgroup atoms first

  if (atom->firstgroup >= 0) {
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int firstgroupbit = group->bitmask[atom->firstgroup];

    int flag = 0;
    for (int i = 0; i < nlocal; i++)
      if ((mask[i] == groupbit) && (mask[i] && firstgroupbit)) flag = 1;

    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);

    if (flagall)
      error->all(FLERR,"Cannot do GCMC on atoms in atom_modify first group");
  }

  if (force->pair->single_enable == 0) 
    error->all(FLERR,"Fix GCMC incompatible with given pair_style");
  
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;  
  
  int size = comm->nprocs;
  
  create(balance_all,size,"Fix_MCSWAP::balance_all");
  create(balance_type1,size,"Fix_MCSWAP::balance_type1");
  create(balance_type2,size,"Fix_MCSWAP::balance_type2");
  create(ngas_local_all,size,"Fix_MCSWAP::ngas_local_all");
  create(ntype1_local_all,size,"Fix_MCSWAP::ntype1_local_all");
  create(ntype2_local_all,size,"Fix_MCSWAP::ntype2_local_all");

  create(sc1x,size,"Fix_MCSWAP::sc1x");
  create(sc1y,size,"Fix_MCSWAP::sc1y");
  create(sc1z,size,"Fix_MCSWAP::sc1z");
  create(sc2x,size,"Fix_MCSWAP::sc2x");
  create(sc2y,size,"Fix_MCSWAP::sc2y");
  create(sc2z,size,"Fix_MCSWAP::sc2z");
  create(svc1x,size,"Fix_MCSWAP::sc1x");
  create(svc1y,size,"Fix_MCSWAP::sc1y");
  create(svc1z,size,"Fix_MCSWAP::sc1z");
  create(svc2x,size,"Fix_MCSWAP::sc2x");
  create(svc2y,size,"Fix_MCSWAP::sc2y");
  create(svc2z,size,"Fix_MCSWAP::sc2z");
  create(sc1tag,size,"Fix_MCSWAP::sc1tag");
  create(sc2tag,size,"Fix_MCSWAP::sc2tag");
  create(tag1,size,"Fix_MCSWAP::tag1");
  create(tag2,size,"Fix_MCSWAP::tag2");
  create(dtag1,size,"Fix_MCSWAP::dtag1");
  create(dtag2,size,"Fix_MCSWAP::dtag2");
  create(dtag1_all,size,"Fix_MCSWAP::dtag1_all");
  create(dtag2_all,size,"Fix_MCSWAP::dtag2_all");
  create(taglist,2*size,"Fix_MCSWAP::taglist");
  create(coordx1,size,"Fix_MCSWAP::coordx1");
  create(coordy1,size,"Fix_MCSWAP::coordy1");
  create(coordz1,size,"Fix_MCSWAP::coordz1");
  create(coordx2,size,"Fix_MCSWAP::coordx2");
  create(coordy2,size,"Fix_MCSWAP::coordy2");
  create(coordz2,size,"Fix_MCSWAP::coordz2");
  create(velx1,size,"Fix_MCSWAP::coordx1");
  create(vely1,size,"Fix_MCSWAP::coordy1");
  create(velz1,size,"Fix_MCSWAP::coordz1");
  create(velx2,size,"Fix_MCSWAP::coordx2");
  create(vely2,size,"Fix_MCSWAP::coordy2");
  create(velz2,size,"Fix_MCSWAP::coordz2");

  if (rank == 0) {
    create(tag1list,size,"Fix_MCSWAP::tag1list");
    create(tag2list,size,"Fix_MCSWAP::tag2list");
    create(denergy1_all,size,"Fix_MCSWAP::denergy1_all");
    create(denergy2_all,size,"Fix_MCSWAP::denergy2_all");
  }
}

/* ---------------------------------------------------------------------- */

void FixMCSWAP::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixMCSWAP::update_atom_list()
{
  int i;

  // grow local_gas_list array if necessary

  if (atom->nlocal > nmax) {
    memory->sfree(local_type1_list);
    memory->sfree(local_type2_list);
    memory->sfree(local_gas_list);
    nmax = atom->nmax;
    local_type1_list = (int *) memory->smalloc(nmax*sizeof(int),
                                             "GCMC:local_type1_list");
    local_type2_list = (int *) memory->smalloc(nmax*sizeof(int),
                                             "GCMC:local_type2_list");
    local_gas_list = (int *) memory->smalloc(nmax*sizeof(int),
                                             "GCMC:local_gas_list");
  }

  int *type = atom->type;
  double **x = atom->x;
  ntype1_local = ntype2_local = 0;
  for (i = 0; i < atom->nlocal; i++) {
    if (x[i][0] > mcxhi || x[i][0] < mcxlo) continue;
    if (type[i] == itype1) {
      local_type1_list[ntype1_local] = i;
      ntype1_local++;
    } else if (type[i] == itype2) {
      local_type2_list[ntype2_local] = i;
      ntype2_local++;
    }
  }

  MPI_Allreduce(&ntype1_local,&ntype1,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&ntype2_local,&ntype2,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&ntype1_local,&ntype1_before,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&ntype2_local,&ntype2_before,1,MPI_INT,MPI_SUM,world);
  ntype1_before -= ntype1_local;
  ntype2_before -= ntype2_local;
  ngas_local = 0;

  for (i = 0; i < atom->nlocal; i++) {
    if (x[i][0] > mcxhi || x[i][0] < mcxlo) continue;
    if (type[i] == itype1 || type[i] == itype2) {
      local_gas_list[ngas_local] = i;
      ngas_local++;
    }
  }
  MPI_Allreduce(&ngas_local,&ngas,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&ngas_local,&ngas_before,1,MPI_INT,MPI_SUM,world);
  ngas_before -= ngas_local;

  MPI_Gather(&ngas_local,1,MPI_INT,ngas_local_all,1,MPI_INT,0,world);
  MPI_Gather(&ntype1_local,1,MPI_INT,ntype1_local_all,1,MPI_INT,0,world);
  MPI_Gather(&ntype2_local,1,MPI_INT,ntype2_local_all,1,MPI_INT,0,world);

  if (rank == 0) {
	int ngas_local_max,ntype1_local_max,ntype2_local_max;
	ngas_local_max = ntype1_local_max = ntype2_local_max = 0;
    for (i = 0; i < size; i++) {
      ngas_local_max  = ngas_local_max < ngas_local_all[i] ? ngas_local_all[i] : ngas_local_max;
      ntype1_local_max = ntype1_local_max < ntype1_local_all[i] ? ntype1_local_all[i] : ntype1_local_max;
      ntype2_local_max = ntype2_local_max < ntype2_local_all[i] ? ntype2_local_all[i] : ntype2_local_max;
    }

    for (i = 0; i < size; i++) {
      balance_all[i] = (double)ngas_local_all[i]/(double)ngas_local_max;
      balance_type1[i] = (double)ntype1_local_all[i]/(double)ntype1_local_max;
      balance_type2[i] = (double)ntype2_local_all[i]/(double)ntype2_local_max;
      //printf("%d %f %f %f \n",i,balance_all[i],balance_type1[i],balance_type2[i]);
    }
  }

  MPI_Barrier(world);
  MPI_Bcast(balance_all,size,MPI_DOUBLE,0,world);
  MPI_Bcast(balance_type1,size,MPI_DOUBLE,0,world);
  MPI_Bcast(balance_type2,size,MPI_DOUBLE,0,world);

}

/* ----------------------------------------------------------------------
   attempt particle insertions and deletions
   done before exchange, borders, reneighbor
   so that ghost atoms and neighbor lists will be correct
------------------------------------------------------------------------- */

void FixMCSWAP::pre_exchange()
{
  // just return if should not be called on this timestep
  
  if (next_reneighbor != update->ntimestep) return;

  if (domain->triclinic == 0) {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    xlo = domain->boxlo_bound[0];
    xhi = domain->boxhi_bound[0];
    ylo = domain->boxlo_bound[1];
    yhi = domain->boxhi_bound[1];
    zlo = domain->boxlo_bound[2];
    zhi = domain->boxhi_bound[2];
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  } 
  
  volume = domain->xprd * domain->yprd * domain->zprd;
  
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();
  update_atom_list();
      
  // perform ncycles MC cycles
  
  if (ncycles > 0) {
    attempt_move();
  }

  if (nswap > 0) {
	domain->pbc();
	comm->exchange();
	atom->nghost = 0;
	comm->borders();
	update_atom_list();
    attempt_swap();
  }
  
  next_reneighbor = update->ntimestep + nevery;
} 

/* ----------------------------------------------------------------------
   choose particle randomly across all procs and attempt displacement
 * Known problems: partial energy incorrect in serial
------------------------------------------------------------------------- */

void FixMCSWAP::attempt_move()
{
  int loop;
  int success,success_all;
  int i,iwhichlocal;
  int j,k;
  double ddenergy;
  double rx,ry,rz;
  double delx,dely,delz;
  double rsq;
  double cutoffsq = cutoff+2.0;
  cutoffsq = cutoffsq*cutoffsq*4.0;
  double coord[3];
  double **x = atom->x;
  int *type = atom->type;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  
  loop = 0;
  while (loop < ncycles) {
    
    success = 0;
    for (j = 0; j < size; j++) {
      dtag1[j] = dtag1_all[j] = 0;
      dtag2[j] = dtag2_all[j] = 0;
      coordx1[j] = coordx2[j] = 0.0;
      coordy1[j] = coordy2[j] = 0.0;
      coordz1[j] = coordz2[j] = 0.0;
    }
    
    // select atom
    iwhichlocal = static_cast<int> (ngas_local*random_equal->uniform());
    i = local_gas_list[iwhichlocal];
    dtag1[rank] = tag[i];
    // atom does not exist?
    if (dtag1[rank] == 0) {
      dtag2[rank]++;
    }

    MPI_Allreduce(dtag1,dtag1_all,size,MPI_INT,MPI_SUM,world);
    
    // 1st screen: probability balance
    if (random_equal->uniform() > balance_all[rank])
      dtag2[rank]++;

    // 2nd screen: distance check

    k = nlocal;
    while (k < nall && dtag2[rank] == 0) {
      for (j = 0; j < size; j++) {
        if (tag[k] == dtag1_all[j]) {
          delx = x[i][0] - x[k][0];
          dely = x[i][1] - x[k][1];
          delz = x[i][2] - x[k][2];
          rsq = delx*delx+dely*dely+delz*delz;
          if (rsq < cutoffsq) {
            dtag2[rank]++;
            break;
          }
        }
      }
      k++;
    }

    MPI_Allreduce(dtag2,dtag2_all,size,MPI_INT,MPI_SUM,world);

    // atom not kicked proceed to MCMD
    if (dtag2_all[rank] == 0) {
      
      for (j = 0; j < size; j++) {
        dtag2[j] = dtag2_all[j];
      }     
      
      // MC movement direction
      double rsq = 1.1;
      while (rsq > 1.0) {
        rx = 2*random_unequal->uniform() - 1.0;
        ry = 2*random_unequal->uniform() - 1.0;
        rz = 2*random_unequal->uniform() - 1.0;
        rsq = rx*rx + ry*ry + rz*rz;
      }
      coord[0] = x[i][0] + displace*rx;
      coord[1] = x[i][1] + displace*ry;
      coord[2] = x[i][2] + displace*rz;
      
      // compute energy change
      ddenergy = denergy(i,type[i],coord);      
      
      if (random_unequal->uniform() < exp(-beta*ddenergy)) {
        x[i][0] = coordx1[rank] = coord[0];
        x[i][1] = coordy1[rank] = coord[1];
        x[i][2] = coordz1[rank] = coord[2];
        success++;
      } else {
        // MCMD not triggered
        dtag2[rank]++;
      }
    }
    
    loop += size;
    nmove_attempts += size;
    
    success_all = 0;
    MPI_Allreduce(&success,&success_all,1,MPI_INT,MPI_SUM,world);
    
    MPI_Allreduce(dtag2,dtag2_all,size,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(coordx1,coordx2,size,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(coordy1,coordy2,size,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(coordz1,coordz2,size,MPI_DOUBLE,MPI_SUM,world);

    if (success_all) {
      nmove_successes += success_all;
      // update ghost atoms
      for (k = nlocal; k < nall; k++) {
        for (j = 0; j < size; j++) {
          if (tag[k] == tag[j] && dtag2_all[j] == 0) {
            x[k][0] = coordx2[j];
            x[k][1] = coordy2[j];
            x[k][2] = coordz2[j];
            break;
          }
        }
      }
    }

  }

}

/* ----------------------------------------------------------------------
   choose atoms randomly with different atom types globally
   and attempt type exchange
------------------------------------------------------------------------- */

void FixMCSWAP::attempt_swap()
{
  int i,j,k;
  int i1,i2;
  int iwhichlocal1;
  int iwhichlocal2;
  double coord1[3],coord2[3];
  double vel1[3],vel2[3];
  double delx,dely,delz;
  double **x = atom->x;
  double **v = atom->v;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int *type = atom->type;
  int *tag = atom->tag;
  
  int size = comm->nprocs;
  int rank = comm->me;
  int element1,element2,ntag1list,ntag2list;
  
  double cutoffsq = cutoff*cutoff*4.0;
  double rsq;

  int loop;
  
  double ddenergy,denergy1,denergy2;
  
  scn = 0;
  
  loop = 0;
  while (loop < nswap) {
    
    iswap = 0;
    scn = 0;
    MPI_Barrier(world);
    
    // try to select two atoms of different types from each procs
    
    for (i = 0; i < size; i++) {
      if (npostswap > 0) {
        sc1x[i] = sc1y[i] = sc1z[i] = 0.0;
        sc2x[i] = sc2y[i] = sc2z[i] = 0.0;
        svc1x[i] = svc1y[i] = svc1z[i] = 0.0;
        svc2x[i] = svc2y[i] = svc2z[i] = 0.0;
        coordx1[i] = coordx2[i] = 0.0;
        coordy1[i] = coordy2[i] = 0.0;
        coordz1[i] = coordz2[i] = 0.0;
        velx1[i] = velx2[i] = 0.0;
        vely1[i] = vely2[i] = 0.0;
        velz1[i] = velz2[i] = 0.0;
      }
      sc1tag[i] = sc2tag[i] = 0;
      tag1[i] = tag2[i] = 0;
      dtag1[i] = dtag1_all[i] = 0;
      dtag2[i] = dtag2_all[i] = 0;
      if (rank == 0) {
        denergy1_all[i] = 0.0;
        denergy2_all[i] = 0.0;
      }
      taglist[i] = 0;
      taglist[i+size] = 0;
    }
    
    i1 = -1;
    iwhichlocal1 = static_cast<int>(ntype1_local*random_equal->uniform());
    i1 = local_type1_list[iwhichlocal1];
    if (i1 > -1) {
      coord1[0] = x[i1][0];
      coord1[1] = x[i1][1];
      coord1[2] = x[i1][2];
      vel1[0] = v[i1][0];
      vel1[1] = v[i1][1];
      vel1[2] = v[i1][2];
      tag1[rank] = tag[i1];
    } else {
      dtag1[rank]++;
      tag1[rank] = 0;
    }

    // 1st screen: probability balance
    if (random_equal->uniform() > balance_type1[rank])
      dtag1[rank]++;

    i2 = -1;
    iwhichlocal2 = static_cast<int>(ntype2_local*random_equal->uniform());
    i2 = local_type2_list[iwhichlocal2];
    if (i2 > -1) {
      coord2[0] = x[i2][0];
      coord2[1] = x[i2][1];
      coord2[2] = x[i2][2];
      vel2[0] = v[i2][0];
      vel2[1] = v[i2][1];
      vel2[2] = v[i2][2];
      tag2[rank] = tag[i2];
    } else {
      dtag2[rank]++;
      tag2[rank] = 0;
    }

    // 1st screen: probability balance
    if (random_equal->uniform() > balance_type2[rank])
      dtag2[rank]++;

    if (i1 > -1 && i2 > -1) {
      delx = coord1[0] - coord2[0];
      dely = coord1[1] - coord2[1];
      delz = coord1[2] - coord2[2];
      rsq = delx*delx+dely*dely+delz*delz;
      if (rsq < cutoffsq) {
        dtag1[rank]++;
        dtag2[rank]++;
      }
    }
    
    MPI_Allreduce(tag1,sc1tag,size,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(tag2,sc2tag,size,MPI_INT,MPI_SUM,world);
    
    // check ghost atoms for possible close target atoms
    i = nlocal;
    while (i < nall && dtag1[rank] == 0) {
      for (j = 0; j < size; j++) {
        if (tag[i] == sc1tag[j] || tag[i] == sc2tag[j]) {
          delx = coord1[0] - x[i][0];
          dely = coord1[1] - x[i][1];
          delz = coord1[2] - x[i][2];
          rsq = delx*delx+dely*dely+delz*delz;
          if (rsq < cutoffsq) {
            dtag1[rank]++;
            break;
          }
        }
      }
      i++;
    }

    i = nlocal;
    while (i < nall && dtag2[rank] == 0) {
      for (j = 0; j < size; j++) {
        if (tag[i] == sc1tag[j] || tag[i] == sc2tag[j]) {
          delx = coord2[0] - x[i][0];
          dely = coord2[1] - x[i][1];
          delz = coord2[2] - x[i][2];
          rsq = delx*delx+dely*dely+delz*delz;
          if (rsq < cutoffsq) {
            dtag2[rank]++;
            break;
          }
        }
      }
      i++;
    }

    MPI_Allreduce(dtag1,dtag1_all,size,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(dtag2,dtag2_all,size,MPI_INT,MPI_SUM,world);

    // calculate partial PE difference
    if (dtag1_all[rank] == 0) {
      denergy1 = denergy(i1,itype2,coord1);
    } else {
      denergy1 = 1.0e30;
    }
    if (dtag2_all[rank] == 0) {
      denergy2 = denergy(i2,itype1,coord2);
    } else {
      denergy2 = 1.0e30;
    }
    
    MPI_Gather(&denergy1,1,MPI_DOUBLE,denergy1_all,1,MPI_DOUBLE,0,world);
    MPI_Gather(&denergy2,1,MPI_DOUBLE,denergy2_all,1,MPI_DOUBLE,0,world);
    MPI_Gather(&vel1[0],1,MPI_DOUBLE,velx1,1,MPI_DOUBLE,0,world);
    MPI_Gather(&vel1[1],1,MPI_DOUBLE,vely1,1,MPI_DOUBLE,0,world);
    MPI_Gather(&vel1[2],1,MPI_DOUBLE,velz1,1,MPI_DOUBLE,0,world);
    MPI_Gather(&vel2[0],1,MPI_DOUBLE,velx2,1,MPI_DOUBLE,0,world);
    MPI_Gather(&vel2[1],1,MPI_DOUBLE,vely2,1,MPI_DOUBLE,0,world);
    MPI_Gather(&vel2[2],1,MPI_DOUBLE,velz2,1,MPI_DOUBLE,0,world);
    if (npostswap > 0) {
      MPI_Gather(&coord1[0],1,MPI_DOUBLE,coordx1,1,MPI_DOUBLE,0,world);
      MPI_Gather(&coord1[1],1,MPI_DOUBLE,coordy1,1,MPI_DOUBLE,0,world);
      MPI_Gather(&coord1[2],1,MPI_DOUBLE,coordz1,1,MPI_DOUBLE,0,world);
      MPI_Gather(&coord2[0],1,MPI_DOUBLE,coordx2,1,MPI_DOUBLE,0,world);
      MPI_Gather(&coord2[1],1,MPI_DOUBLE,coordy2,1,MPI_DOUBLE,0,world);
      MPI_Gather(&coord2[2],1,MPI_DOUBLE,coordz2,1,MPI_DOUBLE,0,world);
    }
    
    // random atom pair and swap
    if (rank == 0) {
      
      k = 0;
      int maxloop = size*size;
      while (k < maxloop) {
        element1 = static_cast<int> (size*random_equal->uniform());
        element2 = static_cast<int> (size*random_equal->uniform());
        k++;
        if (element1 == element2) continue;
        if (dtag1_all[element1] > 0) continue;
        if (dtag2_all[element2] > 0) continue;
        ddenergy = denergy1_all[element1]+denergy2_all[element2];

        //printf("%4d%4d%10d%10d%20.10f\n",element1,element2,sc1tag[element1],sc2tag[element2],ddenergy);
        if (random_unequal->uniform() < exp(-beta*ddenergy)) {

          // found a pair
          taglist[scn] = sc1tag[element1];
          taglist[scn+size] = sc2tag[element2];
          svc1x[scn] = velx1[element1];
          svc1y[scn] = vely1[element1];
          svc1z[scn] = velz1[element1];
          svc2x[scn] = velx2[element2];
          svc2y[scn] = vely2[element2];
          svc2z[scn] = velz2[element2];
          if (npostswap > 0) {
            sc1x[scn] = coordx1[element1];
            sc1y[scn] = coordy1[element1];
            sc1z[scn] = coordz1[element1];
            sc2x[scn] = coordx2[element2];
            sc2y[scn] = coordy2[element2];
            sc2z[scn] = coordz2[element2];
          }
          //printf("%4d%4d%10d%10d%20.10f\n",element1,element2,sc1tag[element1],sc2tag[element2],ddenergy);
          scn++;
          iswap++;
          // accepted atoms cannot be used again
          dtag1_all[element1]++;
          dtag2_all[element2]++;
        }
      }

      /*
      k = 0;
      for (j = 0; j < size; j++) {
        if (dtag1_all[j] == 0) {
          tag1list[k] = j;
          k++;
        }

        // initialize taglist
        taglist[j] = -1;
        taglist[j+size] = -1;
      }
      ntag1list = k;
      
      // think about best useage of low energy.....
      while (ntag1list > 0) {
        
        // rebuild tag2list
        k = 0;
        for (j = 0; j < size; j++) {
          if (dtag2_all[j] == 0) {
            tag2list[k] = j;
            k++;
          }
        }
        ntag2list = k;
        
        // randomly pick one atom from list tag1list
        element1 = static_cast<int> (ntag1list*random_equal->uniform());
        
        // search in tag2list for a potential lower energy pair
        while (ntag2list > 0) {

          // randomly pick on atom from tag2list
          element2 = static_cast<int> (ntag2list*random_equal->uniform());
          // get total energy diff
          ddenergy = denergy1_all[tag1list[element1]] + 
                     denergy2_all[tag2list[element2]];
          printf("%4d%4d%10d%10d%20.10f\n",element1,element2,sc1tag[tag1list[element1]],sc2tag[tag2list[element2]],ddenergy);

          if (random_unequal->uniform() < exp(-beta*ddenergy)) {

            // found a pair
            taglist[scn] = sc1tag[tag1list[element1]];
            taglist[scn+size] = sc2tag[tag2list[element2]];
            svc1x[scn] = velx1[tag1list[element1]];
            svc1y[scn] = vely1[tag1list[element1]];
            svc1z[scn] = velz1[tag1list[element1]];
            svc2x[scn] = velx2[tag2list[element2]];
            svc2y[scn] = vely2[tag2list[element2]];
            svc2z[scn] = velz2[tag2list[element2]];
            if (npostswap > 0) {
              sc1x[scn] = coordx1[tag1list[element1]];
              sc1y[scn] = coordy1[tag1list[element1]];
              sc1z[scn] = coordz1[tag1list[element1]];
              sc2x[scn] = coordx2[tag2list[element2]];
              sc2y[scn] = coordy2[tag2list[element2]];
              sc2z[scn] = coordz2[tag2list[element2]];
            }
            //printf("%4d%4d%10d%10d%20.10f\n",element1,element2,taglist[scn],taglist[scn+size],ddenergy);
            //printf("%4d%10d <=> %4d%10d: denergy: %20.10f\n",element1,taglist[scn],element2,taglist[scn+size],ddenergy);
            scn++;
            iswap++;
            // accepted atom from tag2list cannot be used again
            dtag2_all[tag2list[element2]]++;
            // exit this loop
            ntag2list = 0;

          } else {

            // pop element2 from tag2list
            k = 0;
            for (j = 0; j < ntag2list; j++) {
              if (j != element2) {
                tag2list[k] = tag2list[j];
                k++;
              }
            }
            ntag2list = k;

          }
        }
        
        // pop element1 from tag1list
        k = 0;
        for (j = 0; j < ntag1list; j++) {
          if (j != element1) {
            tag1list[k] = tag1list[j];
            k++;
          }
        }
        ntag1list = k;
      }
      */
    }

    MPI_Barrier(world);
    MPI_Bcast(&scn,1,MPI_INT,0,world);
    MPI_Bcast(&iswap,1,MPI_INT,0,world);
    loop += size;
    nswap_attempts += size;
    nswap_successes += iswap;

    if (iswap > 0) {
      // bcast sc tags to all procs
      MPI_Bcast(taglist,2*size,MPI_INT,0,world);
	  // bcast sc velocites to all procs
      MPI_Bcast(svc1x,size,MPI_DOUBLE,0,world);
      MPI_Bcast(svc1y,size,MPI_DOUBLE,0,world);
      MPI_Bcast(svc1z,size,MPI_DOUBLE,0,world);
      MPI_Bcast(svc2x,size,MPI_DOUBLE,0,world);
      MPI_Bcast(svc2y,size,MPI_DOUBLE,0,world);
      MPI_Bcast(svc2z,size,MPI_DOUBLE,0,world);
      if (npostswap > 0) {
        MPI_Bcast(sc1x,size,MPI_DOUBLE,0,world);
        MPI_Bcast(sc1y,size,MPI_DOUBLE,0,world);
        MPI_Bcast(sc1z,size,MPI_DOUBLE,0,world);
        MPI_Bcast(sc2x,size,MPI_DOUBLE,0,world);
        MPI_Bcast(sc2y,size,MPI_DOUBLE,0,world);
        MPI_Bcast(sc2z,size,MPI_DOUBLE,0,world);
      }
      
      // find target atoms on diff procs
      // both local and ghost atoms updated!
      for (i = 0; i < nall; i++) {
        for (j = 0; j < scn; j++) {
          if (tag[i] == taglist[j]) {
            type[i] = itype2;
			v[i][0] = svc2x[j];
			v[i][1] = svc2y[j];
			v[i][2] = svc2z[j];
          } else if (tag[i] == taglist[j+size]) {
            type[i] = itype1;
			v[i][0] = svc1x[j];
			v[i][1] = svc1y[j];
			v[i][2] = svc1z[j];
          }
        }
      }
      
     if (npostswap > 0) {

        domain->pbc();
        comm->exchange();
        atom->nghost = 0;
        comm->borders();
        update_atom_list();

    	attempt_postswap();

        domain->pbc();
        comm->exchange();
        atom->nghost = 0;
        comm->borders();
        update_atom_list();

      } else {

        // update local_type1_list local_type2_list
        ntype1_local = ntype2_local = 0;
        for (i = 0; i < nlocal; i++) {
          if (x[i][0] > mcxhi || x[i][0] < mcxlo) continue;
          if (type[i] == itype1) {
            local_type1_list[ntype1_local++] = i;
          } else if (type[i] == itype2) {
            local_type2_list[ntype2_local++] = i;
          }
        }

        MPI_Allreduce(&ntype1_local,&ntype1,1,MPI_INT,MPI_SUM,world);
        MPI_Allreduce(&ntype2_local,&ntype2,1,MPI_INT,MPI_SUM,world);
        MPI_Scan(&ntype1_local,&ntype1_before,1,MPI_INT,MPI_SUM,world);
        MPI_Scan(&ntype2_local,&ntype2_before,1,MPI_INT,MPI_SUM,world);
        ntype1_before -= ntype1_local;
        ntype2_before -= ntype2_local;

        MPI_Gather(&ntype1_local,1,MPI_INT,ntype1_local_all,1,MPI_INT,0,world);
        MPI_Gather(&ntype2_local,1,MPI_INT,ntype2_local_all,1,MPI_INT,0,world);

        if (rank == 0) {
          int ntype1_local_max,ntype2_local_max;
          ntype1_local_max = ntype2_local_max = 0;
          for (i = 0; i < size; i++) {
            ntype1_local_max < ntype1_local_all[i] ? ntype1_local_all[i] : ntype1_local_max;
            ntype2_local_max < ntype2_local_all[i] ? ntype2_local_all[i] : ntype2_local_max;
          }
          for (i = 0; i < size; i++) {
            balance_type1[i] = (double)ntype1_local_all[i]/(double)ntype1_local_max;
            balance_type2[i] = (double)ntype2_local_all[i]/(double)ntype2_local_max;
          }
        }

        MPI_Barrier(world);
        MPI_Bcast(balance_type1,size,MPI_DOUBLE,0,world);
        MPI_Bcast(balance_type2,size,MPI_DOUBLE,0,world);
      }
    }
  }
  
}

/* ----------------------------------------------------------------------
   neighborlist MCMD
------------------------------------------------------------------------- */

void FixMCSWAP::attempt_postswap()
{
  int loop;
  int success,success_all;
  int i,iwhichlocal,iwhichglobal;
  int j,k;
  double ddenergy;
  double rx,ry,rz;
  double delx,dely,delz;
  double rsq;
  double cutoffsq = cutoff+2.0;
  cutoffsq = cutoffsq*cutoffsq;
  double coord[3];
  double xtmp,ytmp,ztmp;
  double **x = atom->x;
  int *type = atom->type;
  int *tag = atom->tag;
  int rank = comm->me;
  int size = comm->nprocs;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  
  loop = 0;
  success = 0;

  while (loop <= npostswap) {

    for (j = 0; j < size; j++) {
      dtag1[j] = dtag1_all[j] = 0;
      dtag2[j] = dtag2_all[j] = 0;
      coordx1[j] = coordx2[j] = 0.0;
      coordy1[j] = coordy2[j] = 0.0;
      coordz1[j] = coordz2[j] = 0.0;
    }

    // select atom
    iwhichlocal = static_cast<int> (ngas_local*random_equal->uniform());
    i = local_gas_list[iwhichlocal];
    dtag1[rank] = tag[i];
    // atom does not exist?
    if (dtag1[rank] == 0) {
      dtag2[rank]++;
    }

    MPI_Allreduce(dtag1,dtag1_all,size,MPI_INT,MPI_SUM,world);

    if (dtag2[rank] == 0) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];

      // check if atom belongs to any sc
      for (j = 0; j < scn; j++) {
        delx = xtmp - sc1x[j];
        dely = ytmp - sc1y[j];
        delz = ztmp - sc1z[j];
        rsq = delx*delx+dely*dely+delz*delz;
        if (rsq < cutoffsq)
          dtag2[rank]++;
        delx = xtmp - sc2x[j];
        dely = ytmp - sc2y[j];
        delz = ztmp - sc2z[j];
        rsq = delx*delx+dely*dely+delz*delz;
        if (rsq < cutoffsq)
          dtag2[rank]++;
      }
    }

    // inverse logic
    if (dtag2[rank] == 0) {
      dtag2[rank] = 1;
    } else {
      dtag2[rank] = 0;
    }

    // 1st screen: probability balance
    if (random_equal->uniform() > balance_all[rank])
      dtag2[rank]++;

    // check atom distance across different procs via ghost atoms

    k = nlocal;
    while (k < nall && dtag2[rank] == 0) {
      for (j = 0; j < size; j++) {
        if (tag[k] == dtag1_all[j]) {
          delx = xtmp - x[k][0];
          dely = ytmp - x[k][1];
          delz = ztmp - x[k][2];
          rsq = delx*delx+dely*dely+delz*delz;
          if (rsq < cutoffsq) {
            dtag2[rank]++;
            break;
          }
        }
      }
      k++;
    }

    MPI_Allreduce(dtag2,dtag2_all,size,MPI_INT,MPI_SUM,world);

    // atom not kicked proceed to MCMD
    if (dtag2_all[rank] == 0) {
      
      for (j = 0; j < size; j++) {
        dtag2[j] = dtag2_all[j];
      }     
      
      // MC movement direction
      double rsq = 1.1;
      while (rsq > 1.0) {
        rx = 2*random_unequal->uniform() - 1.0;
        ry = 2*random_unequal->uniform() - 1.0;
        rz = 2*random_unequal->uniform() - 1.0;
        rsq = rx*rx + ry*ry + rz*rz;
      }
      coord[0] = xtmp + displace*rx;
      coord[1] = ytmp + displace*ry;
      coord[2] = ztmp + displace*rz;
      
      // compute energy change
      ddenergy = denergy(i,type[i],coord);

      if (random_unequal->uniform() < exp(-beta*ddenergy)) {
        x[i][0] = coordx1[rank] = coord[0];
        x[i][1] = coordy1[rank] = coord[1];
        x[i][2] = coordz1[rank] = coord[2];
        success++;
      } else {
        // MCMD not triggered
        dtag2[j]++;
      }
    }
    
    loop += size;
    npostswap_attempts += size;
    
    success_all = 0;
    MPI_Allreduce(&success,&success_all,1,MPI_INT,MPI_SUM,world);
    
    MPI_Allreduce(dtag2,dtag2_all,size,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(coordx1,coordx2,size,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(coordy1,coordy2,size,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(coordz1,coordz2,size,MPI_DOUBLE,MPI_SUM,world);

    if (success_all) {
      npostswap_successes += success_all;
      
      // update ghost atoms
      for (k = nlocal; k < nall; k++) {
        for (j = 0; j < size; j++) {
          if (tag[k] == tag[j] && dtag2_all[j] == 0) {
            x[k][0] = coordx2[j];
            x[k][1] = coordy2[j];
            x[k][2] = coordz2[j];
            break;
          }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute particle's interaction energy with the rest of the system
------------------------------------------------------------------------- */

double FixMCSWAP::denergy(int i, int itype, double *coord) 
{
  double delx,dely,delz,rsq;
  double coord0[3];
  double **x = atom->x;
  int *type = atom->type;
  int nall = atom->nlocal + atom->nghost;
  pair = force->pair;
  cutsq = force->pair->cutsq;
  double cutoff = neighbor->cutneighmax;
  double cutoffsq1 = cutoff+2.0;
  cutoffsq1 = cutoffsq1*cutoffsq1;
  double cutoffsq2 = cutoff*cutoff;
  
  int itype0 = type[i];
  int tflag = itype0 == itype ? 0 : 1;
  coord0[0] = x[i][0];
  coord0[1] = x[i][1];
  coord0[2] = x[i][2];

  double fpair = 0.0;
  double factor_coul = 1.0;
  double factor_lj = 1.0;
  
  double total_energy = 0.0;
  double dpair,dembedded;
  
  dpair = dembedded = 0.0;
  
  // Move
  if (!tflag) {
    for (int j = 0; j < nall; j++) {

      if (i == j) continue;

      delx = coord0[0] - x[j][0];
      dely = coord0[1] - x[j][1];
      delz = coord0[2] - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      int jtype = type[j];

      if (rsq < cutoffsq1) {
        dpair -= 
      pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);
        delx = coord[0] - x[j][0];
        dely = coord[1] - x[j][1];
        delz = coord[2] - x[j][2];
        rsq = delx*delx+dely*dely+delz*delz;
        dpair += 
      pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);
      }
    }
  // type switch
  } else {
    for (int j = 0; j < nall; j++) {
      
      if (i == j) continue;

      delx = coord0[0] - x[j][0];
      dely = coord0[1] - x[j][1];
      delz = coord0[2] - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      int jtype = type[j];

      if (rsq < cutoffsq2) {
        dpair -= 
      pair->single(i,j,itype0,jtype,rsq,factor_coul,factor_lj,fpair);
      }
      if (rsq < cutoffsq2) {
        dpair += 
      pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);
      }
    }
  }
  
  dembedded = pair->singleembedding(i, itype, coord);
  
  total_energy = dpair + dembedded;
  
  return total_energy;
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixMCSWAP::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal fix MCSWAP command");

  int iarg = 0;
  while (iarg < narg) {
	if (strcmp(arg[iarg],"type") == 0) {
      itype1 = atoi(arg[++iarg]);
      itype2 = atoi(arg[++iarg]);
      iarg++;
	} else if (strcmp(arg[iarg],"mcmd") == 0) {
	  nmcmoves = atoi(arg[++iarg]);
	  iarg++;
	} else if (strcmp(arg[iarg],"mcswap") == 0) {
      nswap = atoi(arg[++iarg]);
      npostswap = atoi(arg[++iarg]);
      iarg++;
    } else if (strcmp(arg[iarg],"mcxlh") == 0) {
      mcxlo = atof(arg[++iarg]);
      mcxhi = atof(arg[++iarg]);
      iarg++;
    } else error->all(FLERR,"Illegal fix MCSWAP command");
  }
}

/* ----------------------------------------------------------------------
  return acceptance ratios
------------------------------------------------------------------------- */

double FixMCSWAP::compute_vector(int n)
{
  double value;
  switch (n) {
    case 0:
      value = nmove_attempts;
      nmove_attempts = 0.0;
      break;
    case 1:
      value = nmove_successes;
      nmove_successes = 0.0;
      break;
    case 2:
      value = npostswap_attempts;
      npostswap_attempts = 0.0;
      break;
    case 3:
      value = npostswap_successes;
      npostswap_successes = 0.0;
      break;
    case 4:
      value = nswap_attempts;
      nswap_attempts = 0.0;
      break;
    case 5:
      value = nswap_successes;
      nswap_successes = 0.0;
      break;
    default:
      value = 0.0;
      break;
  }
  return value;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixMCSWAP::memory_usage()
{
  double bytes = nmax * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write 
------------------------------------------------------------------------- */

void FixMCSWAP::write_restart(FILE *fp)
{
  int n = 0;
  double list[4];
  list[n++] = random_equal->state();
  list[n++] = random_unequal->state();
  list[n++] = next_reneighbor;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix  
------------------------------------------------------------------------- */

void FixMCSWAP::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  seed = static_cast<int> (list[n++]);
  random_equal->reset(seed);
  
  seed = static_cast<int> (list[n++]);
  random_unequal->reset(seed);
  
  next_reneighbor = static_cast<int> (list[n++]);
}
