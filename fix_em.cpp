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
   Contributing authors: Xiaoyin Ji
------------------------------------------------------------------------- */

#include "stdio.h"
#include "string.h"
#include "fix_em.h"
#include "atom.h"
#include "input.h"
#include "variable.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "group.h"
#include "domain.h"
#include "update.h"
#include "pointers.h"
#include "memory.h"

#include "petsc.h"
#include "petscsys.h"
#include "petscmat.h"
#include "petscpc.h"
#include "udspc.h"
#include "array.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,SPHERE,ELLIPSOID,PYRAMID};

/* ---------------------------------------------------------------------- */

FixEM::FixEM(LAMMPS *lmp, int narg, char **arg) :
Fix(lmp, narg, arg)
{
  static char help[] = "Appends to an ASCII file.\n\n";
  PetscInitialize(&narg,&arg,(char *)0,help);
  MPI_Comm_rank(world,&rank);
  MPI_Comm_size(world,&size);
  
  amu2kg = 1.66054E-27;  
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  boost = 1;
  restart_global = 1;
  igroup = 0;
  bcflag = iflag = 0;
  xstr = ystr = zstr = rstr = pstr = NULL;
  astr = bstr = cstr = NULL;
  xvalue = yvalue = zvalue = rvalue = pvalue = 0.0;
  avalue = bvalue = cvalue = 0.0;
  ctr[0] = ctr[1] = ctr[2] = 0.0;
  radiusr = 0.0;
  atc = NULL;
  grid = NULL;
  ngridcnt = ngridcnt_thermo = NULL;
  pa = NULL;
  gridKE = gridCP = gridV = NULL;  
  intTemp1 = intTemp2 = NULL;
  doubleTemp = NULL;
  TgridSAVE = TgridREF = TgridRT = NULL;
  eres = tres = NULL;
  gridI = gridE = gridVec = NULL;
  gridVSAVE = NULL;
  FF = DFF = NULL;
  restarted = false;
  periodic_yz = true; // default: FD periodic boundary on  
  gridVecflag = false; // default: translational velocity off
  fdonly = false;
  for (int i = 0; i < 16; i++) {
    para_select[i] = false;
    para_dump[i] = 1000; // default: paraview dump every 1000 timesteps
  }  
  delay = -1; // default: start FD steps right now  
  fort71 = 10; // default: fort.71 every 10 timesteps  
  fddump = 1000; // default: special dump every 1000 timesteps
  ctypel = ctypeu = 1;
  time_integrate = 1;
  atype1 = 1;
  atype2 = 2;
  istyle = NONE;
}

/* ---------------------------------------------------------------------- */

FixEM::~FixEM()
{
  destroy(ngridcnt);
  destroy(ngridcnt_thermo);
  destroy(pa);
  destroy(gridCP);
  destroy(gridKE);
  destroy(TgridSAVE);
  destroy(TgridREF);
  destroy(TgridRT);
  destroy(gridV);
  destroy(gridVSAVE);
  destroy(eres);
  destroy(tres);
  destroy(gridI);
  destroy(gridE);
  destroy(gridVec);  
  destroy(intTemp1);
  destroy(intTemp2);
  destroy(doubleTemp);
  
  if (!fdonly) {
    destroy(FF);
    destroy(DFF);
    destroy(grid);
  }

  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] rstr;
  delete [] pstr;

  PetscFinalize();
}

/* ---------------------------------------------------------------------- */

int FixEM::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write 
------------------------------------------------------------------------- */

void FixEM::write_restart(FILE *fp)
{
  int nsize = size_restart_global();
  
  double *list;
  memory->create(list,nsize,"em:list");
  
  pack_restart_data(list);
  
  if (rank == 0) {
    int size = nsize * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),nsize,fp);
  }
  
  memory->destroy(list);
}

/* ----------------------------------------------------------------------
   calculate the number of data to be packed
------------------------------------------------------------------------- */

int FixEM::size_restart_global()
{
  int nsize = ngridtot * 4 + 1;
  
  return nsize;
}

/* ----------------------------------------------------------------------
   pack restart data 
------------------------------------------------------------------------- */

int FixEM::pack_restart_data(double *list)
{
  int i,n;
  i = 0; n = 0;
  list[n++] = static_cast<double>(ngridtot);
  while (i < ngridtot) {
    list[n++] = TgridSAVE[i];
    list[n++] = pa[i];
    list[n++] = static_cast<double>(ngridcnt[i]);
    list[n++] = static_cast<double>(ngridcnt_thermo[i]);
    i++;
  }
  
  return n;
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix 
------------------------------------------------------------------------- */

void FixEM::restart(char *buf)
{
  int i,n;
  double *list = (double *) buf;
  n = 0;
  ngridtot_restart = static_cast<int>(list[n++]);
  restarted = true;  
  create(TgridSAVE, ngridtot_restart, "electromeam:TgridSAVE");
  create(TgridREF,  ngridtot_restart, "electromeam:TgridREF");
  create(TgridRT,   ngridtot_restart, "electromeam:TgridRT");
  create(ngridcnt,        ngridtot_restart, "electromeam:ngridcnt");
  create(ngridcnt_thermo, ngridtot_restart, "electromeam:ngridcnt_thermo");
  create(pa,              ngridtot_restart, "electromeam:pa");
  for (i = 0; i < ngridtot_restart; i++) {
    // set Tgrid* to t_request for the first run
    TgridSAVE[i] = TgridREF[i] = TgridRT[i] = list[n++];
    pa[i] = list[n++];
    ngridcnt[i] = static_cast<int>(list[n++]);
    ngridcnt_thermo[i] = static_cast<int>(list[n++]);
  }
}

/* ---------------------------------------------------------------------- */

int FixEM::modify_param(int narg, char** arg)
{
  int iarg = 1;
  int idata = 0;
  
  // FD only
  if (strcmp(arg[0],"fdonly") == 0) {
    if (narg != 2) 
      error->all(FLERR, "ELECTROMEAM: Illegal fdonly command");
    if (strcmp(arg[1],"yes") == 0 || strcmp(arg[1],"on") == 0) {
      fdonly = true;
    } else if (strcmp(arg[1],"no") == 0 || strcmp(arg[1],"off") == 0) {
      fdonly = false;
    } else {
      error->all(FLERR, "ELECTROMEAM: Illegal fdonly command");
    }
    return narg;

  // Run FD for every N MD steps
  } else if (strcmp(arg[0],"boost") == 0) {
    if (narg != 2) 
      error->all(FLERR, "ELECTROMEAM: Illegal boost command");
    boost = atoi(arg[1]);
    if (boost < 1)
      error->all(FLERR, "ELECTROMEAM: Boost level of FD cannot below 1");
    dtv *= static_cast<double>(boost);
    return narg;

  // alloy type
  } else if (strcmp(arg[0],"type") == 0) {
    if (narg != 3)
      error->all(FLERR, "ELECTROMEAM: Illegal atom type command");
    atype1 = atoi(arg[1]);
    atype2 = atoi(arg[2]);
    if (rank == 0) {
      if (screen) fprintf(screen,"primary metal type: %d, secondary metal type: %d\n", atype1, atype2);
      if (logfile) fprintf(logfile,"primary metal type: %d, secondary metal type: %d\n", atype1, atype2);
    }
    return narg;

  // group settings
  } else if (strcmp(arg[0],"em_all") == 0) {
    if (narg != 2) 
      error->all(FLERR, "ELECTROMEAM: Illegal em_all command");
    int findgroup = group->find(arg[iarg++]);
    if (findgroup == -1)
      error->all(FLERR, "Cound not find FD group ID");
    groupbit[0] = group->bitmask[findgroup];
    return narg;
  } else if (strcmp(arg[0],"em_free") == 0) {
    if (narg != 2) 
      error->all(FLERR, "ELECTROMEAM: Illegal em_free command");
    int findgroup = group->find(arg[iarg++]);
    if (findgroup == -1) 
      error->all(FLERR, "Cound not find FD group ID");
    groupbit[1] = group->bitmask[findgroup];
    return narg;
  } else if (strcmp(arg[0],"em_BC") == 0) {
    if (narg != 2) 
      error->all(FLERR, "ELECTROMEAM: Illegal em_BC command");
    int findgroup = group->find(arg[iarg++]);
    if (findgroup == -1) 
      error->all(FLERR, "Cound not find FD group ID");
    groupbit[2] = group->bitmask[findgroup];
	bcflag = 1;
	return narg;

  // Ghost indenter
  // automatically enable dynamic BC
  } else if (strcmp(arg[0],"indent") == 0) {
	if (narg > 8)
	  error->all(FLERR, "ELECTROMEAM: Illegal indent command");
    
      // spherical indenter
	if (strcmp(arg[iarg],"sphere") == 0) {
	  if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
	    int n = strlen(&arg[iarg+1][2]) + 1;
	    xstr = new char[n];
        strcpy(xstr,&arg[iarg+1][2]);
      } else xvalue = atof(arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
        int n = strlen(&arg[iarg+2][2]) + 1;
        ystr = new char[n];
        strcpy(ystr,&arg[iarg+2][2]);
      } else yvalue = atof(arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) {
        int n = strlen(&arg[iarg+3][2]) + 1;
        zstr = new char[n];
        strcpy(zstr,&arg[iarg+3][2]);
      } else zvalue = atof(arg[iarg+3]);
      if (strstr(arg[iarg+4],"v_") == arg[iarg+4]) {
        int n = strlen(&arg[iarg+4][2]) + 1;
        rstr = new char[n];
        strcpy(rstr,&arg[iarg+4][2]);
      } else rvalue = atof(arg[iarg+4]);
    
      istyle = SPHERE;
    
      // elliptical indenter
	} else if (strcmp(arg[iarg],"ellipsoid") == 0) {
	  if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        xstr = new char[n];
        strcpy(xstr,&arg[iarg+1][2]);
      } else xvalue = atof(arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
        int n = strlen(&arg[iarg+2][2]) + 1;
        ystr = new char[n];
        strcpy(ystr,&arg[iarg+2][2]);
      } else yvalue = atof(arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) {
        int n = strlen(&arg[iarg+3][2]) + 1;
        zstr = new char[n];
        strcpy(zstr,&arg[iarg+3][2]);
      } else zvalue = atof(arg[iarg+3]);
      if (strstr(arg[iarg+4],"v_") == arg[iarg+4]) {
        int n = strlen(&arg[iarg+4][2]) + 1;
        astr = new char[n];
        strcpy(astr,&arg[iarg+4][2]);
      } else avalue = atof(arg[iarg+4]);
      if (strstr(arg[iarg+5],"v_") == arg[iarg+5]) {
        int n = strlen(&arg[iarg+5][2]) + 1;
        bstr = new char[n];
        strcpy(bstr,&arg[iarg+5][2]);
      } else bvalue = atof(arg[iarg+5]);
      if (strstr(arg[iarg+6],"v_") == arg[iarg+6]) {
        int n = strlen(&arg[iarg+6][2]) + 1;
        cstr = new char[n];
        strcpy(cstr,&arg[iarg+6][2]);
      } else cvalue = atof(arg[iarg+6]);
    
      istyle = ELLIPSOID;
    
      // pyramid indenter
	} else if (strcmp(arg[iarg],"pyramid") == 0) {
	  if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        xstr = new char[n];
        strcpy(xstr,&arg[iarg+1][2]);
      } else xvalue = atof(arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
        int n = strlen(&arg[iarg+2][2]) + 1;
        ystr = new char[n];
        strcpy(ystr,&arg[iarg+2][2]);
      } else yvalue = atof(arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) {
        int n = strlen(&arg[iarg+3][2]) + 1;
        zstr = new char[n];
        strcpy(zstr,&arg[iarg+3][2]);
      } else zvalue = atof(arg[iarg+3]);
      if (strstr(arg[iarg+4],"v_") == arg[iarg+4]) {
        int n = strlen(&arg[iarg+4][2]) + 1;
        rstr = new char[n];
        strcpy(rstr,&arg[iarg+4][2]);
      } else rvalue = atof(arg[iarg+4]);
      if (strstr(arg[iarg+5],"v_") == arg[iarg+5]) {
        int n = strlen(&arg[iarg+5][2]) + 1;
        astr = new char[n];
        strcpy(astr,&arg[iarg+5][2]);
      } else avalue = atof(arg[iarg+5]);

      istyle = PYRAMID;

    }

	bcflag = 2;
	iflag = 1;
	return narg;

  // FD data output
  } else if (strcmp(arg[0],"dataout") == 0) {
    if (narg > 4) 
      error->all(FLERR, "ELECTROMEAM: Illegal dataout command");
    if (strcmp(arg[1],"fort71") == 0) {
      fort71 = atoi(arg[2]);
    } else if (strcmp(arg[1],"fddump") == 0) {
      fddump = atoi(arg[2]);
    } else if (strcmp(arg[1],"para") == 0) { 
      if (strcmp(arg[2], "temp") == 0) {
        para_select[0] = true;
        para_dump[0] = atoi(arg[3]);
      } else if (strcmp(arg[2],"eres") == 0) {
        para_select[1] = true;
        para_dump[1] = atoi(arg[3]);
      } else if (strcmp(arg[2],"tres") == 0) {
        para_select[2] = true;
        para_dump[2] = atoi(arg[3]);
      } else if (strcmp(arg[2],"current") == 0) {
        para_select[3] = true;
        para_dump[3] = atoi(arg[3]);
      } else if (strcmp(arg[2],"voltage") == 0) {
        para_select[4] = true;
        para_dump[4] = atoi(arg[3]);
      } else if (strcmp(arg[2],"gridke") == 0) {
        para_select[5] = true;
        para_dump[5] = atoi(arg[3]);
      } else if (strcmp(arg[2],"gridcp") == 0) {
        para_select[6] = true;
        para_dump[6] = atoi(arg[3]);
      }
    } else {
      error->all(FLERR, "ELECTROMEAM: thermo property undefined");
      return 0;
    }
    return narg;

  // FD periodic
  } else if (strcmp(arg[0],"fd_periodic") == 0) {
    if (narg != 2) 
      error->all(FLERR, "ELECTROMEAM: Illegal fd_periodic sub command");
    if (strcmp(arg[1],"yes") == 0 || strcmp(arg[1],"on") == 0) 
      periodic_yz = true;
    else if (strcmp(arg[1],"no") == 0 || strcmp(arg[1],"off") == 0) 
      periodic_yz = false;
    else 
      error->warning(FLERR, "ELECTROMEAM: incorrect perioidc_yz sub command, periodic off by default");
    return narg;

  // FD grid velocity
  } else if (strcmp(arg[0],"fd_gridVc") == 0) {
    if (narg != 2) 
      error->all(FLERR, "ELECTROMEAM: Illegal fd_gridVc sub command");
    if (strcmp(arg[1],"no") == 0 || strcmp(arg[1],"off") == 0) 
      gridVecflag = true;
    else if (strcmp(arg[1],"no") == 0 || strcmp(arg[1],"off") == 0) 
      gridVecflag = false;
    else 
      error->warning(FLERR, "ELECTROMEAM: incorrect fd_gridVc sub command, fd_gridVc on by default");
    return narg;

  // FD setup
  } else if (strcmp(arg[0],"bounds") == 0) {
    if (narg != 7) 
      error->all(FLERR, "ELECTROMEAM: Illegal bounds sub command");
    grid_lb[0] = atof(arg[iarg++]);
    grid_ub[0] = atof(arg[iarg++]);
    grid_lb[1] = atof(arg[iarg++]);
    grid_ub[1] = atof(arg[iarg++]);
    grid_lb[2] = atof(arg[iarg++]);
    grid_ub[2] = atof(arg[iarg++]);
    return narg;
  } else if (strcmp(arg[0],"continuum_xbounds") == 0) {
    if (narg != 5) 
      error->all(FLERR, "ELECTROMEAM: Illegal continuum_xbounds sub command");
    cxlower = atof(arg[iarg++]);
    cxupper = atof(arg[iarg++]);
    ctypel = atoi(arg[iarg++]);
    ctypeu = atoi(arg[iarg++]);
    if ((ctypel != atype1 && ctypel != atype2) || (ctypeu != atype1 && ctypeu != atype2))
      error->all(FLERR, "Illegal ctypel and/or ctypeu value");
    return narg;
  } else if (strcmp(arg[0],"data") == 0) {
    if (narg != 8) 
      error->all(FLERR, "ELECTROMEAM: Illegal data command");
    idata = atoi(arg[iarg++]);
    if (idata == atype1) {
      FDtc[0] = atof(arg[iarg++]);
      FDcp[0] = atof(arg[iarg++]);
      FDdens[0] = atof(arg[iarg++]);
      FDrho[0] = atof(arg[iarg++]);
      FDrhoT[0] = atof(arg[iarg++]);
      FDalpharho[0] = atof(arg[iarg++]);
    } else if (idata == atype2) {
      FDtc[1] = atof(arg[iarg++]);
      FDcp[1] = atof(arg[iarg++]);
      FDdens[1] = atof(arg[iarg++]);
      FDrho[1] = atof(arg[iarg++]);
      FDrhoT[1] = atof(arg[iarg++]);
      FDalpharho[1] = atof(arg[iarg++]);
    }
    if (rank == 0) {
      if (idata == atype1) {
        if (screen) {
           fprintf(screen,"primary metal properties: \n");
           fprintf(screen,"thermal conductivity: %10.4f heat capacity: %10.4f\n", FDtc[0],FDcp[0]);
           fprintf(screen,"density:              %10.4f rho:           %10.4f\n", FDdens[0],FDrho[0]);
           fprintf(screen,"rhoT:                 %10.4f alpharho     : %10.4f\n", FDrhoT[0],FDalpharho[0]);
        }
        if (logfile)  {
           fprintf(logfile,"primary metal properties: \n");
           fprintf(logfile,"thermal conductivity: %10.4f heat capacity: %10.4f\n", FDtc[0],FDcp[0]);
           fprintf(logfile,"density:              %10.4f rho:           %10.4f\n", FDdens[0],FDrho[0]);
           fprintf(logfile,"rhoT:                 %10.4f alpharho     : %10.4f\n", FDrhoT[0],FDalpharho[0]);
         }
      } else if (idata == atype2) {
        if (screen) {
          fprintf(screen,"secondary metal properties: \n");
          fprintf(screen,"thermal conductivity: %10.4f heat capacity: %10.4f\n", FDtc[1],FDcp[1]);
          fprintf(screen,"density:              %10.4f rho:           %10.4f\n", FDdens[1],FDrho[1]);
          fprintf(screen,"rhoT:                 %10.4f alpharho     : %10.4f\n", FDrhoT[1],FDalpharho[1]);
        }
        if (logfile)  {
          fprintf(logfile,"secondary metal properties: \n");
          fprintf(logfile,"thermal conductivity: %10.4f heat capacity: %10.4f\n", FDtc[1],FDcp[1]);
          fprintf(logfile,"density:              %10.4f rho:           %10.4f\n", FDdens[1],FDrho[1]);
          fprintf(logfile,"rhoT:                 %10.4f alpharho     : %10.4f\n", FDrhoT[1],FDalpharho[1]);
        }
      } else {
        error->all(FLERR, "ELECTROMEAM: metal type incorrect");
      }
    }
    return narg;
  } else if (strcmp(arg[0],"grid_ra") == 0) {
    if (narg != 4) 
      error->all(FLERR, "ELECTROMEAM: Illegal grid_ra command");
    ngrids[0][0] = atoi(arg[iarg++]);
    ngrids[0][1] = atoi(arg[iarg++]);
    ngrids[0][2] = atoi(arg[iarg++]);
    return narg;
  } else if (strcmp(arg[0],"grid_rb") == 0) {
    if (narg != 4) 
      error->all(FLERR, "ELECTROMEAM: Illegal grid_ra command");
    ngrids[1][0] = atoi(arg[iarg++]);
    ngrids[1][1] = atoi(arg[iarg++]);
    ngrids[1][2] = atoi(arg[iarg++]);
    return narg;
  } else if (strcmp(arg[0],"grid_rc") == 0) {
    if (narg != 4) 
      error->all(FLERR, "ELECTROMEAM: Illegal grid_ra command");
    ngrids[2][0] = atoi(arg[iarg++]);
    ngrids[2][1] = atoi(arg[iarg++]);
    ngrids[2][2] = atoi(arg[iarg++]);
    return narg;
  } else if (strcmp(arg[0],"t_request") == 0) {
    if (narg != 2) 
      error->all(FLERR, "ELECTROMEAM: Illegal t_request command");
    t_request = atof(arg[iarg++]);
    return narg;
  } else if (strcmp(arg[0],"voltage") == 0) {
    if (narg != 3) 
      error->all(FLERR, "ELECTROMEAM: Illegal voltage command");
    fd_vlv = atof(arg[iarg++]);
    fd_vuv = atof(arg[iarg++]);
    return narg;
  } else {
    error->all(FLERR, "ELECTROMEAM: ? Illegal command");
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixEM::init()
{
  if (strstr(update->integrate_style,"respa")) {
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
    step_respa = ((Respa *) update->integrate)->step;
  }
  
  ngridtot = 0;
  int i = 0;
  for (i = 0; i < 3; i++) {
    ngrideach[i] = ngrids[i][0]*ngrids[i][1]*ngrids[i][2];
    ngridtot += ngrideach[i];
  }
  
  if (restarted) {
    if (ngridtot != ngridtot_restart) {
      char line[256];
      sprintf(line, "restarted ngridtot mismatch: new%10d, old%10d", 
              ngridtot, ngridtot_restart);
      error->all(FLERR, line);
    }
  }
  
  double *mass = atom->mass;
  MASS_VAL[0] = mass[atype1];
  MASS_VAL[1] = mass[atype2];
  
  //ngridcnt = new int [ngridtot];
  create(atc,        ngridtot, "electromeam:atc");
  create(gridKE,     ngridtot, "electromeam:gridKE");
  create(gridCP,     ngridtot, "electromeam:gridCP");
  create(gridV,      ngridtot, "electromeam:gridV");
  create(gridVSAVE,  ngridtot, "electromeam:gridV");
  create(intTemp1,   ngridtot, "electromeam:intTemp1");
  create(intTemp2,   ngridtot, "electromeam:intTemp2");
  create(doubleTemp, ngridtot, "electromeam:doubleTemp");
  create(eres,       ngridtot, 3, "electromeam:eres");
  create(tres,       ngridtot, 3, "electromeam:tres");
  create(gridI,      ngridtot, 3, "electromeam:gridI");
  create(gridE,      ngridtot, 3, "electromeam:gridE");
  
  if (!restarted) {
    create(TgridSAVE, ngridtot, "electromeam:TgridSAVE");
    create(TgridREF,  ngridtot, "electromeam:TgridREF");
    create(TgridRT,   ngridtot, "electromeam:TgridRT");
    create(ngridcnt,        ngridtot, "electromeam:ngridcnt");
    create(ngridcnt_thermo, ngridtot, "electromeam:ngridcnt_thermo");
    create(pa,              ngridtot, "electromeam:pa");
  }
  
  if (!fdonly) {
    create(FF, ngridtot, "electromeam:FF");
    create(DFF,ngridtot, "electromeam:DFF");
  }
  if (gridVecflag)
    create(gridVec, ngridtot,3,"electromeam::gridVec");
  
  // warning is not enough, boundaries must match!
  if (ngrids[1][1]%ngrids[0][1]!=0)
    error->all(FLERR, "Y SIZE MISSMATCH BETWEEN REGION A AND B");
  if (ngrids[1][1]%ngrids[0][2]!=0)
    error->all(FLERR, "Z SIZE MISSMATCH BETWEEN REGION A AND B");
  if (ngrids[1][1]%ngrids[2][1]!=0)
    error->all(FLERR, "Y SIZE MISSMATCH BETWEEN REGION C AND B");
  if (ngrids[1][1]%ngrids[2][2]!=0)
    error->all(FLERR, "Z SIZE MISSMATCH BETWEEN REGION C AND B");
  
  // grid set parameters
  // cwp boundary mehtod ?
  double plen[3];
  double small = 1.0e-4;
  
  for (int j = 0; j < 3; j++) {
    plen[j] = grid_ub[j] - grid_lb[j];
    if (j == 0) plen[j] = grid_ub[j] - cxupper;
    grid_size[2][j] = (plen[j]+small*plen[j])/(double)ngrids[2][j];
  }
  for (int j = 0; j < 3; j++) {
    plen[j] = grid_ub[j] - grid_lb[j];
    if (j == 0) plen[j] = cxupper - cxlower;
    grid_size[1][j] = (plen[j]+small*plen[j])/(double)ngrids[1][j];
  }
  for (int j = 0; j < 3; j++) {
    plen[j] = grid_ub[j] - grid_lb[j];
    if (j == 0) plen[j] = cxlower - grid_lb[j];
    grid_size[0][j] = (plen[j]+small*plen[j])/(double)ngrids[0][j];
  }

  // dynamic ATC boundary BC region
  // default value: cxupper
  BC = cxupper;
  
  // preparing grid R factor and grid volume
  for (int i = 0; i < 3; i++) {
    grid_Rfactor[i][0] = grid_size[i][0]/grid_size[i][1]/grid_size[i][2]*1.0E+10;
    grid_Rfactor[i][1] = grid_size[i][1]/grid_size[i][0]/grid_size[i][2]*1.0E+10;
    grid_Rfactor[i][2] = grid_size[i][2]/grid_size[i][0]/grid_size[i][1]*1.0E+10;
    grid_volume[i] = grid_size[i][0]*grid_size[i][1]*grid_size[i][2]*1.0E-30;
  }
  
  // region matching constants
  BOAY = ngrids[1][1] / ngrids[0][1];
  BOAZ = ngrids[1][2] / ngrids[0][2];
  BOCY = ngrids[1][1] / ngrids[2][1];
  BOCZ = ngrids[1][2] / ngrids[2][2];
  ABCouple = BOAY*BOAZ;
  BCCouple = BOCY*BOCZ;
  ABCmax = MAX(ABCouple,BCCouple);
  
  // grid voltage guess for preconditioner
  double gridVSAVE_default = (fd_vlv+fd_vuv)/2;
  
  for (i = 0; i < ngridtot; i++) {    
    // set Tgrid* to t_request for the first run
    if (!restarted)
      TgridSAVE[i] = TgridREF[i] = TgridRT[i] = t_request;    
      // grid voltage guess for preconditioner
      gridVSAVE[i] = gridVSAVE_default;
      // grid velocity
      if (gridVecflag) 
        gridVec[i][0] = gridVec[i][1] = gridVec[i][2] = 0.0;    
      // grid kinetic energy
      gridKE[i] = 0.;
  }
  
  if (iflag == 1) {
    if (xstr) {
      xvar = input->variable->find(xstr);
      if (xvar < 0) 
        error->all(FLERR,"Variable name for fix indent does not exist");
      if (!input->variable->equalstyle(xvar))
        error->all(FLERR,"Variable for fix indent is invalid style");
    }
    if (ystr) {
      yvar = input->variable->find(ystr);
      if (yvar < 0) 
        error->all(FLERR,"Variable name for fix indent does not exist");
      if (!input->variable->equalstyle(yvar))
        error->all(FLERR,"Variable for fix indent is not equal style");
    }
    if (zstr) {
      zvar = input->variable->find(zstr);
      if (zvar < 0) 
        error->all(FLERR,"Variable name for fix indent does not exist");
      if (!input->variable->equalstyle(zvar))
        error->all(FLERR,"Variable for fix indent is not equal style");
    }
    if (rstr) {
      rvar = input->variable->find(rstr);
      if (rvar < 0) 
        error->all(FLERR,"Variable name for fix indent does not exist");
      if (!input->variable->equalstyle(rvar))
        error->all(FLERR,"Variable for fix indent is not equal style");
    }
    if (pstr) {
      pvar = input->variable->find(pstr);
      if (pvar < 0) 
        error->all(FLERR,"Variable name for fix indent does not exist");
      if (!input->variable->equalstyle(pvar))
        error->all(FLERR,"Variable for fix indent is not equal style");
    }
    if (astr) {
      avar = input->variable->find(astr);
      if (avar < 0) 
        error->all(FLERR,"Variable name for fix indent does not exist");
      if (!input->variable->equalstyle(avar))
        error->all(FLERR,"Variable for fix indent is not equal style");
    }
    if (bstr) {
      bvar = input->variable->find(bstr);
      if (bvar < 0) 
        error->all(FLERR,"Variable name for fix indent does not exist");
      if (!input->variable->equalstyle(bvar))
        error->all(FLERR,"Variable for fix indent is not equal style");
    }
    if (cstr) {
      cvar = input->variable->find(cstr);
      if (cvar < 0) 
        error->all(FLERR,"Variable name for fix indent does not exist");
      if (!input->variable->equalstyle(cvar))
        error->all(FLERR,"Variable for fix indent is not equal style");
    }
  }
  
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  
  gridmax = nlocal;
  create(grid, gridmax, "electromeam:grid");
  
  if (fdonly) {
    if (!restarted) 
      check_atoms();
    destroy(DFF);
    destroy(FF);
    destroy(grid);
  }
  return;
}

/* ---------------------------------------------------------------------- */

void FixEM::setup(int vflag)
{
  FILE *fpfort71;  
  fpfort71 = fopen("fort.71","w");
  fclose(fpfort71);
  
  // delete atoms and create one atom
  if (fdonly && !restarted) {
    
    int *dlist;
    // store state before delete    
    bigint natoms_previous = atom->natoms;
    
    // delete the atoms    
    int igroup = group->find("all");
    
    // allocate and initialize deletion list    
    int nlocal = atom->nlocal;
    create(dlist,nlocal,"electromeam:dlist");
    for (int i = 0; i < nlocal; i++) dlist[i] = 0;
    
    int *mask = atom->mask;
    int groupbit = group->bitmask[igroup];
    
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) dlist[i] = 1;
    
    // keep one atom on each proc
    if (nlocal)
      dlist[0] = 0;
    
    // delete local atoms flagged in dlist
    // reset nlocal
    
    AtomVec *avec = atom->avec;
    
    int i = 0;
    while (i < nlocal) {
      if (dlist[i]) {
        avec->copy(nlocal-1,i,1);
        dlist[i] = dlist[nlocal-1];
        nlocal--;
      } else i++;
    }
    
    atom->nlocal = nlocal;
    destroy(dlist);
    
    // reset atom->natoms
    // reset atom->map if it exists
    // set nghost to 0 so old ghosts of deleted atoms won't be mapped
    
    bigint nblocal = atom->nlocal;
    MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
    if (atom->map_style) {
      atom->nghost = 0;
      atom->map_init();
      atom->map_set();
    }
    
    // print before and after atom count
    
    bigint ndelete = natoms_previous - atom->natoms;
    
    if (comm->me == 0) {
      if (screen) fprintf(screen,"Deleted " BIGINT_FORMAT 
                          " atoms, new total = " BIGINT_FORMAT "\n",
                          ndelete,atom->natoms);
      if (logfile) fprintf(logfile,"Deleted " BIGINT_FORMAT 
                           " atoms, new total = " BIGINT_FORMAT "\n",
                           ndelete,atom->natoms);
    }
  }
}

/* ----------------------------------------------------------------------
   return the gridnumber
------------------------------------------------------------------------- */

int FixEM::get_grid_num(int i, int j, int k, int region)
{
  int num = -1;
  switch (region) {
    case 0:
      num = i*ngrids[0][1]*ngrids[0][2]+j*ngrids[0][2]+k;
      break;
    case 1:
      num = i*ngrids[1][1]*ngrids[1][2]+j*ngrids[1][2]+k+ngrideach[0];
      break;
    case 2:
      num = i*ngrids[2][1]*ngrids[2][2]+j*ngrids[2][2]+k+ngrideach[0]+ngrideach[1];
      break;
    default:
      break;
  }
  return num;
}

/* ---------------------------------------------------------------------- */

void FixEM::check_atoms()
{
  int i,j,BCtempN,BCtempNsum;
  int *type = atom->type;
  int *mask = atom->mask;
  int xperiodic = domain->xperiodic;
  int yperiodic = domain->yperiodic;
  int zperiodic = domain->zperiodic;
  double BCtemp,BCtemp1,BCtemp2;
  double ada[3];
  double *lo,*hi,*period;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  BCtempN = BCtempNsum = 0;
  BCtemp = BCtemp1 = BCtemp2 = 0.0;

  if (gridmax < nlocal) {
    gridmax = nlocal;
    grow(grid, gridmax, "electromeam:grid");
  }
  
  if (domain->triclinic == 0) {
    lo = domain->boxlo;
    hi = domain->boxhi;
    period = domain->prd;
  } else {
    lo = domain->boxlo_lamda;
    hi = domain->boxhi_lamda;
    period = domain->prd_lamda;
  }

  for (i = 0; i < ngridtot; i++) {
    intTemp1[i] = ngridcnt[i] = 0;
    intTemp2[i] = ngridcnt_thermo[i] = 0;
    doubleTemp[i] = pa[i] = 0.0;
  }
  
  // tricky part, lammps do not always check pbc
  // domain->pbc();
  // Atom Loop for all
  for (j = 0; j < nlocal; j++) {
    
    // if atom is outside, i value set to -1
    grid[j] = -1;
    if (!(mask[j] & groupbit[0])) continue;
    // assume 1 and 2 are for alloy	  
    if (type[j] != atype1 && type[j] != atype2) continue;
    // save atom->x to ada
    ada[0] = x[j][0]; ada[1] = x[j][1]; ada[2] = x[j][2];
    
    if (xperiodic) {
      if (ada[0] <  lo[0]) ada[0] += period[0];
      if (ada[0] >= hi[0]) {
        ada[0] -= period[0];
        ada[0] = MAX(ada[0],lo[0]);
      }
    }
    
    if (yperiodic) {
      if (ada[1] <  lo[1]) ada[1] += period[1];
      if (ada[1] >= hi[1]) {
        ada[1] -= period[1];
        ada[1] = MAX(ada[1],lo[1]);
      }
    }
    
    if (zperiodic) {
      if (ada[2] <  lo[2]) ada[2] += period[2];
      if (ada[2] >= hi[2])  {
        ada[2] -= period[2];
        ada[2] = MAX(ada[2],lo[2]);
      }
    }
    
    // cwp skips atoms out size box !modifided for multiple grid reagions
    if (ada[0] < grid_lb[0] || ada[0] > grid_ub[0]) continue;
    if (ada[1] < grid_lb[1] || ada[1] > grid_ub[1]) continue;
    if (ada[2] < grid_lb[2] || ada[2] > grid_ub[2]) continue;
    
    if (ada[0] >= grid_lb[0] && ada[0] < cxlower)    iabc = 0;
    if (ada[0] >= cxlower    && ada[0] < cxupper)    iabc = 1;
    if (ada[0] >= cxupper    && ada[0] < grid_ub[0]) iabc = 2;
    
    //  int rounds down!
    if (iabc == 0) {
      ixx = (int)((ada[0] - grid_lb[0])/grid_size[0][0]);
      iyy = (int)((ada[1] - grid_lb[1])/grid_size[0][1]);
      izz = (int)((ada[2] - grid_lb[2])/grid_size[0][2]);
    }
    
    if (iabc == 1) {
      ixx = (int)((ada[0] - cxlower   )/grid_size[1][0]);
      iyy = (int)((ada[1] - grid_lb[1])/grid_size[1][1]);
      izz = (int)((ada[2] - grid_lb[2])/grid_size[1][2]);
    }
    
    if (iabc == 2) {
      ixx = (int)((ada[0] - cxupper   )/grid_size[2][0]);
      iyy = (int)((ada[1] - grid_lb[1])/grid_size[2][1]);
      izz = (int)((ada[2] - grid_lb[2])/grid_size[2][2]);
    }
    
    i = get_grid_num(ixx,iyy,izz,iabc);
    grid[j] = i;
    
    intTemp1[i]++;
    if (mask[j] & groupbit[1]) intTemp2[i]++;
    if (type[j] == atype2) pa[i]++; // percentage of atype2 !!

    if (bcflag == 1) {
	    // dynamic ATC boundary BC region
    	// buggy and don't use for now
      if (mask[j] & groupbit[2]) {
        if (BCtempN == 0) {
          BCtemp1 = ada[0];
          BCtempN++;
        }
        BCtemp2 = ada[0];
        if (BCtemp2 != BCtemp1)
          error->all(FLERR, "inconsistant x position of dynamic ATC boundary on region BC");
      }
    }
  }

  if (bcflag == 1) {
    MPI_Allreduce(&BCtempN,&BCtempNsum,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&BCtemp1,&BCtemp,1,MPI_DOUBLE,MPI_SUM,world);
    BCtemp /= (double)BCtempNsum;

    // dynamic ATC boundary inside region B!
    if (BCtemp < cxupper) BC = BCtemp;
  }

  // finalize ngridcnt
  MPI_Allreduce(intTemp1,ngridcnt,ngridtot,MPI_INT,MPI_SUM,world);
  // finalize ngridcnt_thermo
  MPI_Allreduce(intTemp2,ngridcnt_thermo,ngridtot,MPI_INT,MPI_SUM,world);
  
  for (i = 0; i < ngridtot; i++) {
    if(ngridcnt[i] != 0) 
      doubleTemp[i] = pa[i]/ngridcnt[i];
  }
  
  MPI_Allreduce(doubleTemp,pa,ngridtot,MPI_DOUBLE,MPI_SUM,world);
  
  // if continuum regions are the same, just return here
  if (ctypel == atype1 && ctypeu == atype1) return;
  
  for (i = 0; i < ngridtot; i++) {
    
    // set iabc only
    if (i < ngrideach[0]) iabc = 0;
    else if (i >= ngrideach[0] && i < ngrideach[0] + ngrideach[1]) iabc = 1;
    else iabc = 2;
    
  	if (ngridcnt[i] == 0) {
      if (iabc == 0) pa[i] = ctypel == atype1 ? 0.0 : 1.0;
      if (iabc == 2) pa[i] = ctypeu == atype1 ? 0.0 : 1.0;
  	}
  }

  return;
}

/* ---------------------------------------------------------------------- */

void FixEM::check_atc()
{
  int set,i,ixx,iyy,izz;
  double dx,dy,dz;
  
  // initialize atc
  for (i = 0; i < ngridtot; i++) {
    if (i < ngrideach[0]) {
      iabc = 0;
      atc[i] = 1;
    } else if (i >= ngrideach[0] && i < ngrideach[0] + ngrideach[1]) {
      iabc = 1;
      atc[i] = 0;
    } else {
      iabc = 2;
      atc[i] = 1;
    }
  }
  
  // spherical indenter
  if (istyle == SPHERE) {
    // ctr = current indenter center
    // remap into periodic box

    if (xstr) ctr[0] = input->variable->compute_equal(xvar);
    else ctr[0] = xvalue;
    if (ystr) ctr[1] = input->variable->compute_equal(yvar);
    else ctr[1] = yvalue;
    if (zstr) ctr[2] = input->variable->compute_equal(zvar);
    else ctr[2] = zvalue;
    domain->remap(ctr);

    double radius;
    if (rstr) radius = input->variable->compute_equal(rvar);
    else radius = rvalue;
    radiusr = radius*radius;
  
    if (iflag == 1) { BC = ctr[0]; }

    for (i = 0; i < ngridtot; i++) {
    
	    // region A
      if (i < ngrideach[0]) {
	      iabc = 0;
	      atc[i] = 1;

      // region B
      } else if (i >= ngrideach[0] && i < ngrideach[0] + ngrideach[1]) { 
	      iabc = 1;
	      atc[i] = 0;
	      // check dynamic bc
	      if (bcflag > 0) {
          // check dynamic ATC boundary BC region
	        if (BC < cxupper) {
	          set = i - ngrideach[0];
            ixx = set/(ngrids[iabc][1]*ngrids[iabc][2]);
		        if (BC < cxlower + (ixx+1)*grid_size[1][0]) {
		          atc[i] = 1;
		        }
	        }
	      }

        // check dummy indenter
        if (iflag == 1) {
          set = i - ngrideach[0];
          ixx = set/(ngrids[iabc][1]*ngrids[iabc][2]);
          iyy = (set - ixx*ngrids[iabc][1]*ngrids[iabc][2])/ngrids[iabc][2];
          izz = set - ixx*ngrids[iabc][1]*ngrids[iabc][2] - iyy*ngrids[iabc][2];
          dx = ((double)ixx + 0.5)*grid_size[1][0]-ctr[0]+cxlower;
          dy = ((double)iyy + 0.5)*grid_size[1][1]-ctr[1]+grid_lb[1];
          dz = ((double)izz + 0.5)*grid_size[1][2]-ctr[2]+grid_lb[2];
          // inside indenter, convert into continuum
          if (dx*dx+dy*dy+dz*dz <= radiusr) {
            atc[i] = 1;
          }
        }
     
        // check frozen atoms
        if (ngridcnt[i] != 0 && ngridcnt[i] != ngridcnt_thermo[i]) {
          atc[i] = 1;
        }

      // region C
      } else {
        iabc = 2;
        atc[i] = 1;
      }
    
    }
   
  // ellipsoid indenter

  } else if (istyle == ELLIPSOID) {
    // ctr = current indenter center
    // remap into periodic box
    
    if (xstr) ctr[0] = input->variable->compute_equal(xvar);
    else ctr[0] = xvalue;
    if (ystr) ctr[1] = input->variable->compute_equal(yvar);
    else ctr[1] = yvalue;
    if (zstr) ctr[2] = input->variable->compute_equal(zvar);
    else ctr[2] = zvalue;
    domain->remap(ctr);
    
    double apara;
    if (astr) apara = input->variable->compute_equal(avar);
    else apara = avalue;
    double bpara;
    if (bstr) bpara = input->variable->compute_equal(bvar);
    else bpara = bvalue;
    double cpara;
    if (cstr) cpara = input->variable->compute_equal(cvar);
    else cpara = cvalue;
    
    double rho = 0.0;

    if (iflag == 1) { BC = ctr[0]; }
    
    for (i = 0; i < ngridtot; i++) {
      
	    // region A
      if (i < ngrideach[0]) {
	      iabc = 0;
	      atc[i] = 1;
        
      // region B
      } else if (i >= ngrideach[0] && i < ngrideach[0] + ngrideach[1]) { 
	      iabc = 1;
	      atc[i] = 0;
	      // check dynamic bc
	      if (bcflag > 0) {
          // check dynamic ATC boundary BC region
	        if (BC < cxupper) {
	          set = i - ngrideach[0];
            ixx = set/(ngrids[iabc][1]*ngrids[iabc][2]);
		        if (BC < cxlower + (ixx+1)*grid_size[1][0]) {
		          atc[i] = 1;
		        }
	        }
	      }
        
        // check dummy indenter
        if (iflag == 1) {
          set = i - ngrideach[0];
          ixx = set/(ngrids[iabc][1]*ngrids[iabc][2]);
          iyy = (set - ixx*ngrids[iabc][1]*ngrids[iabc][2])/ngrids[iabc][2];
          izz = set - ixx*ngrids[iabc][1]*ngrids[iabc][2] - iyy*ngrids[iabc][2];
          dx = ((double)ixx + 0.5)*grid_size[1][0]-ctr[0]+cxlower;
          dy = ((double)iyy + 0.5)*grid_size[1][1]-ctr[1]+grid_lb[1];
          dz = ((double)izz + 0.5)*grid_size[1][2]-ctr[2]+grid_lb[2];
          rho = dx*dx/apara/apara + dy*dy/bpara/bpara + dz*dz/cpara/cpara;
          // inside indenter, convert into continuum
          if (rho <= 1.0) {
            atc[i] = 1;
          }
        }
        
        // check frozen atoms
        if (ngridcnt[i] != 0 && ngridcnt[i] != ngridcnt_thermo[i]) {
          atc[i] = 1;
        }
        
        // region C
      } else {
        iabc = 2;
        atc[i] = 1;
      }

    }

  // pyramid indenter

  } else if (istyle == PYRAMID) {
    // ctr = current indenter center
    // remap into periodic box

    if (xstr) ctr[0] = input->variable->compute_equal(xvar);
    else ctr[0] = xvalue;
    if (ystr) ctr[1] = input->variable->compute_equal(yvar);
    else ctr[1] = yvalue;
    if (zstr) ctr[2] = input->variable->compute_equal(zvar);
    else ctr[2] = zvalue;
    domain->remap(ctr);

    double radius;
    if (rstr) radius = input->variable->compute_equal(rvar);
    else radius = rvalue;

    double apara;
    if (astr) apara = input->variable->compute_equal(avar);
    else apara = avalue;

    double rho = 0.0;

    if (iflag == 1) { BC = ctr[0]; }
  
    for (i = 0; i < ngridtot; i++) {

	  // region A
      if (i < ngrideach[0]) {
	    iabc = 0;
	    atc[i] = 1;

      // region B
      } else if (i >= ngrideach[0] && i < ngrideach[0] + ngrideach[1]) {
	    iabc = 1;
	    atc[i] = 0;
	    // check dynamic bc
	    if (bcflag > 0) {
          // check dynamic ATC boundary BC region
	      if (BC < cxupper) {
	        set = i - ngrideach[0];
            ixx = set/(ngrids[iabc][1]*ngrids[iabc][2]);
		    if (BC < cxlower + (ixx+1)*grid_size[1][0]) {
		      atc[i] = 1;
		    }
	      }
	    }

        // check dummy indenter
        if (iflag == 1) {
          set = i - ngrideach[0];
          ixx = set/(ngrids[iabc][1]*ngrids[iabc][2]);
          iyy = (set - ixx*ngrids[iabc][1]*ngrids[iabc][2])/ngrids[iabc][2];
          izz = set - ixx*ngrids[iabc][1]*ngrids[iabc][2] - iyy*ngrids[iabc][2];
          dx = ((double)ixx + 0.5)*grid_size[1][0]+cxlower;
          dy = ((double)iyy + 0.5)*grid_size[1][1]+grid_lb[1];
          dz = ((double)izz + 0.5)*grid_size[1][2]+grid_lb[2];
          atc[i] = 1;
          if (dx - ctr[0] + apara <= 0) atc[i] = 0;
          if (dx - dy + radius + ctr[1] - ctr[0] <= 0) atc[i] = 0;
          if (dx + dy + radius - ctr[1] - ctr[0] <= 0) atc[i] = 0;
          if (dx - dz + radius + ctr[2] - ctr[0] <= 0) atc[i] = 0;
          if (dx + dz + radius - ctr[2] - ctr[0] <= 0) atc[i] = 0;
        }

        // check frozen atoms
        if (ngridcnt[i] != 0 && ngridcnt[i] != ngridcnt_thermo[i]) {
          atc[i] = 1;
        }

      // region C
      } else {
        iabc = 2;
        atc[i] = 1;
      }
    }
  }
  return;
}

/* ---------------------------------------------------------------------- 
   set the center of grid mass velocity
   communication dense part
------------------------------------------------------------------------- */

void FixEM::set_grid_velocity()
{
  if (!gridVecflag) return;

  int i,j;
  double **tempvecs;
  double **v = atom->v;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  double p,tmass;
  
  create(tempvecs, ngridtot, 3, "electromeam::tempvecs");
  
  // reset zero
  for (i = 0; i < ngridtot; i++) {
    tempvecs[i][0] = 0.;
    tempvecs[i][1] = 0.;
    tempvecs[i][2] = 0.;
    gridVec[i][0] = gridVec[i][1] = gridVec[i][2] = 0.;
  }
  
  for (j = 0; j < nlocal; j++) {
    i = grid[j];
    if (atc[i] == 1) continue;
	  if (i > 0 && (mask[j] & groupbit[1])) {
      tempvecs[i][0] += mass[type[j]]*v[j][0];
      tempvecs[i][1] += mass[type[j]]*v[j][1];
      tempvecs[i][2] += mass[type[j]]*v[j][2];
    }
  }
  
  MPI_Allreduce(*tempvecs,*gridVec,ngridtot*3,MPI_DOUBLE,MPI_SUM,world);
  
  for (i = 0; i < ngridtot; i++) {
    if (ngridcnt[i] > 1) {
      p = pa[i];
      tmass = ((1-p)*MASS_VAL[0] + p*MASS_VAL[1])*ngridcnt[i];
      gridVec[i][0] /= tmass;
      gridVec[i][1] /= tmass;
      gridVec[i][2] /= tmass;
    } else if (ngridcnt[i] <= 1) {
      // if there is only one atom, grid velocity is zero
      gridVec[i][0] = gridVec[i][1] = gridVec[i][2] = 0.;
    }
  }
  
  destroy(tempvecs);
}

/* ---------------------------------------------------------------------- 
  calculate grid kinetic energy basd on atoms inside
------------------------------------------------------------------------- */

void FixEM::set_grid_ke()
{
  int i,j;
  double **v = atom->v;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double mvv2e = force->mvv2e; // conversion of mv^2 to energy
  double ke;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  
  for (i = 0; i < ngridtot; i++)
    doubleTemp[i] = 0.;
  
  if (gridVecflag) {
    for (j = 0; j < nlocal; j++) {
      i = grid[j];
      // only consider unconstrained atoms which are also in FD grids
      if (i >= 0 && (mask[j] & groupbit[1])) {
        ke  = (v[j][0] - gridVec[i][0])*(v[j][0] - gridVec[i][0]);
        ke += (v[j][1] - gridVec[i][1])*(v[j][1] - gridVec[i][1]);
        ke += (v[j][2] - gridVec[i][2])*(v[j][2] - gridVec[i][2]);
        ke *= mass[type[j]]*0.5*mvv2e;
        doubleTemp[i] += ke;
      }
    }
  } else {
    for (j = 0; j < nlocal; j++) {
      i = grid[j];
      // only consider unconstrained atoms which are also in FD grids
      if (i >= 0 && (mask[j] & groupbit[1])) {
        ke = v[j][0]*v[j][0] + v[j][1]*v[j][1] + v[j][2]*v[j][2];
        ke *= mass[type[j]]*0.5*mvv2e;
        doubleTemp[i] += ke;
      }
    }
  }
  
  MPI_Allreduce(doubleTemp,gridKE,ngridtot,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- 
   update grid temperature before solving current
   also read restart TgridSAVE
------------------------------------------------------------------------- */

void FixEM::set_grid_temp()
{
  int i,ixx,iyy,izz,set;
  double iBC = cxlower;
  double t_ext = t_request;
  double boltz = force->boltz; // Boltzmann constant (eng/degree-K)
  
  for (i = 0; i < ngridtot; i++) {

    // region A
    if (i < ngrideach[0]) {
      iabc = 0;
      // continuum region, inherit grid temperature
      TgridREF[i] = TgridSAVE[i];

	// region B
    } else if (i >= ngrideach[0] && i < ngrideach[0] + ngrideach[1]) { 
      iabc = 1;
	  if (bcflag > 0) {
	    // check dynamic ATC boundary BC region
	    set = i - ngrideach[0];
        ixx = set/(ngrids[iabc][1]*ngrids[iabc][2]);
	    iBC = cxlower + (ixx+1)*grid_size[1][0];
	  } else {
	    iBC = cxlower;
	  }
	  if (BC < iBC && bcflag > 0) {
	    // continuum grid inside region B!
	    // inherit grid temperature for continuum part in region B
		t_ext = TgridSAVE[i];
		if (t_ext ==  t_request) { // vaccum, inherit from region C
		  // set i into ixx,iyy,izz,iabc
          iyy = (set - ixx*ngrids[iabc][1]*ngrids[iabc][2])/ngrids[iabc][2];
          izz = set - ixx*ngrids[iabc][1]*ngrids[iabc][2] - iyy*ngrids[iabc][2];
          if (iBC == cxupper) {
		    iyy = iyy*ngrids[iabc+1][1]/ngrids[iabc][1];
			izz = izz*ngrids[iabc+1][2]/ngrids[iabc][2];
			set = ngrideach[0] + ngrideach[1] + iyy*ngrids[iabc+1][2] + izz;
		  } else {
			set = ngrideach[0] + (ixx+1)*ngrids[iabc][1]*ngrids[iabc][2] + iyy*ngrids[iabc][2] + izz;
		  }
		  TgridREF[i] = TgridSAVE[set];
		} else {
          TgridREF[i] = t_ext;
		}
	  } else {
        if (atc[i] == 0) {
          if (ngridcnt[i] == 0) {
			// vacuum
            TgridRT[i] = TgridREF[i] = TgridSAVE[i] = t_request;
		  } else {
            // set grid temperature from unconstrained atoms
            // skip if frozen atoms detected
			if (ngridcnt_thermo[i] > 1) {
              TgridREF[i] = 2./3./boltz*gridKE[i]/ngridcnt_thermo[i];
			} else {
			  TgridREF[i] = TgridSAVE[i];
			}
		  }
		} else {
          // continuum
          TgridREF[i] = TgridSAVE[i];
        }
	  }

    // region C
    } else {
      iabc = 2;
      // continuum region, inherit grid temperature
      TgridREF[i] = TgridSAVE[i];
    }
  }
}

/* ----------------------------------------------------------------------
   updates grid temperature after thermostat
   also writes restart TgridSAVE
------------------------------------------------------------------------- */

void FixEM::update_grid_temp()
{
  int i;
  //double boltz = force->boltz; // Boltzmann constant (eng/degree-K)

  for (i = 0; i < ngridtot; i++) {
    TgridSAVE[i] = TgridRT[i];
  }
}

/* ----------------------------------------------------------------------
   setting grid electric resistivity, SI units, ohm
   Warning: changing the order of calculation may affect precision
------------------------------------------------------------------------- */

void FixEM::set_grid_eres()
{
  int i;
  double p,TFDdens,TFDrho,TFDrhoT,TFDalpharho,phoc;
  double tmass,phoratio,temprho,fd_getR;
  
  for (i = 0; i < ngridtot; i++) {
    
    if (i < ngrideach[0]) {
	  iabc = 0;
	  p = ctypel == atype1 ? 0.0 : 1.0;
	} else if (i >= ngrideach[0] && i < ngrideach[0] + ngrideach[1]) { 
	  iabc = 1;
	  p = atc[i] == 0 ? pa[i] : ctypeu;
	} else {
	  iabc = 2;
	  p = ctypeu == atype1 ? 0.0 : 1.0;
	}
    
    tmass = (1-p)*MASS_VAL[0] + p*MASS_VAL[1];
    TFDdens = (1-p)*FDdens[0] + p*FDdens[1];
    TFDrho = (1-p)*FDrho[0] + p*FDrho[1];
    TFDrhoT = (1-p)*FDrhoT[0] + p*FDrhoT[1];
    TFDalpharho = (1-p)*FDalpharho[0] + p*FDalpharho[1];
    
    phoc = ngridcnt[i]*tmass*amu2kg;
    phoc = phoc / (grid_size[iabc][0]*grid_size[iabc][1]*grid_size[iabc][2]) / 1.0e-30;
    
    phoratio = (atc[i] == 0) ? (phoc/TFDdens) : 1.0;
    temprho = TFDrho*(1.0 + TFDalpharho*(TgridREF[i] - TFDrhoT));
    if(temprho < 1.0e-12) temprho = 1.0e-12; // minimum resistivity is 1.0E-12
    
    fd_getR = temprho/(phoratio+1.0e-15);
    
    // empty grid is assumed to be filled with air
    eres[i][0] = fd_getR*grid_size[iabc][0]/(grid_size[iabc][1]*grid_size[iabc][2]) * 1.0e10;
    eres[i][1] = fd_getR*grid_size[iabc][1]/(grid_size[iabc][0]*grid_size[iabc][2]) * 1.0e10;
    eres[i][2] = fd_getR*grid_size[iabc][2]/(grid_size[iabc][0]*grid_size[iabc][1]) * 1.0e10;

 //   if (update->ntimestep == 100 && i == 12258 && rank == 0) {
 //     cout << "ngridcnt[i] " << ngridcnt[i] << endl;
 //     cout << "atc[i] " << atc[i] << endl;
 //     cout << "phoc " << phoc << endl;
 //     cout << "phoratio " << phoratio << endl;
 //     cout << "temprho " << temprho << endl;
 //     cout << "fd_getR " << fd_getR << endl;
 //     cout << "eres[i][0] " << eres[i][0] << endl;
  //  }
  }

}

/* ----------------------------------------------------------------------
   setting grid thermal resistivity ? conductivity, SI, W / (m^2*K)
   Warning: changing the order of calculation may affect precision
------------------------------------------------------------------------- */

void FixEM::set_grid_tres()
{
  int i;
  double p,TFDdens,TFDtc;
  double tmass,phoi;
  
  for (i = 0; i < ngridtot; i++) {
    
    if (i < ngrideach[0]) {
	  iabc = 0;
	  p = ctypel == atype1 ? 0.0 : 1.0;
	} else if (i >= ngrideach[0] && i < ngrideach[0] + ngrideach[1]) { 
	  iabc = 1;
	  p = atc[i] == 0 ? pa[i] : ctypeu;
	} else {
	  iabc = 2;
	  p = ctypeu == atype1 ? 0.0 : 1.0;
	}

    tmass = (1-p)*MASS_VAL[0] + p*MASS_VAL[1];
    TFDdens = (1-p)*FDdens[0] + p*FDdens[1];
    TFDtc = (1-p)*FDtc[0] + p*FDtc[1];
    
    // actual density v.s. average density
    // think about creating a new array of i to hold the densities?
    // anyway, tmass should be accurate enough
    // 1.0e-5 empty grid is assumed to be filled with air
    
    if (atc[i] == 1) {// continuum
      tres[i][0] = grid_size[iabc][0]/TFDtc/(grid_size[iabc][1]*grid_size[iabc][2])/1.0e-10;
      tres[i][1] = grid_size[iabc][1]/TFDtc/(grid_size[iabc][0]*grid_size[iabc][2])/1.0e-10;
      tres[i][2] = grid_size[iabc][2]/TFDtc/(grid_size[iabc][0]*grid_size[iabc][1])/1.0e-10;
    } else {
      phoi = ngridcnt[i]*tmass*amu2kg;
      phoi /= grid_size[1][0]*grid_size[1][1]*grid_size[1][2]*1.0e-30;
      tres[i][0] = grid_size[iabc][0]/TFDtc/(grid_size[iabc][1]*grid_size[iabc][2])/1.0e-10*(TFDdens/(phoi+1.0e-5));
      tres[i][1] = grid_size[iabc][1]/TFDtc/(grid_size[iabc][0]*grid_size[iabc][2])/1.0e-10*(TFDdens/(phoi+1.0e-5));
      tres[i][2] = grid_size[iabc][2]/TFDtc/(grid_size[iabc][0]*grid_size[iabc][1])/1.0e-10*(TFDdens/(phoi+1.0e-5));
    }
  }
}

/* ----------------------------------------------------------------------
   setting grid heat capacity, SI, J / K
------------------------------------------------------------------------- */

void FixEM::set_grid_cp()
{
  int i;
  double p,tmass,TFDdens,TFDcp,phoi;
  
  for (i = 0; i < ngridtot; i++) {
    
    if (i < ngrideach[0]) {
	  iabc = 0;
	  p = ctypel == atype1 ? 0.0 : 1.0;
	} else if (i >= ngrideach[0] && i < ngrideach[0] + ngrideach[1]) { 
	  iabc = 1;
	  p = atc[i] == 0 ? pa[i] : ctypeu;
	} else {
	  iabc = 2;
	  p = ctypeu == atype1 ? 0.0 : 1.0;
	}
	
    tmass = (1-p)*MASS_VAL[0] + p*MASS_VAL[1];
    TFDdens = (1-p)*FDdens[0] + p*FDdens[1];
    TFDcp = (1-p)*FDcp[0] + p*FDcp[1];
    
    if (atc[i] == 1) { // ingridcnt(i) /= 0
      gridCP[i] = TFDdens*TFDcp*(grid_size[iabc][0]*grid_size[iabc][1]*grid_size[iabc][2])*1.0e-30;
    } else {
      gridCP[i] = ngridcnt[i]*tmass*amu2kg*TFDcp;
    }
  }
}

/* ----------------------------------------------------------------------
   calculate resistive heat, update grid temperature
 ------------------------------------------------------------------------- */

void FixEM::jouleheat()
{
  int i;
  double Q;
  
  for (i = 0; i < ngridtot; i++) {

    Q  = eres[i][0]*gridI[i][0]*gridI[i][0];
    Q += eres[i][1]*gridI[i][1]*gridI[i][1];
    Q += eres[i][2]*gridI[i][2]*gridI[i][2];
    Q *= dtv*1e-12;
    Q = (gridCP[i] == 0) ? 0. : Q/gridCP[i];
    
    TgridRT[i] += Q;
  }
}

/* ----------------------------------------------------------------------
   velocity and force rescaling with updated grid temperature
------------------------------------------------------------------------- */

void FixEM::thermostat()
{
  int i,j;
  double tsqrt; // velocity rescaling factor
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  
  for (i = 0; i < ngridtot; i++)
    FF[i] = DFF[i] = 0.;
  
  for (j = 0; j < nlocal; j++) {  
    if (!(mask[j] & groupbit[1])) continue;
    i = grid[j];
    if (i >= 0) {
      tsqrt = sqrt(TgridRT[i]/TgridREF[i]);
      if (gridVecflag) {
        v[j][0] -= gridVec[i][0];
        v[j][1] -= gridVec[i][1];
        v[j][2] -= gridVec[i][2];
      }
      // velocity rescaling
      v[j][0] *= tsqrt;
      v[j][1] *= tsqrt;
      v[j][2] *= tsqrt;
      // sum friction forces for Hoover for each grid box
      // force*velocity part of Hoover drag
      FF[i]  += f[j][0]*v[j][0];
      FF[i]  += f[j][1]*v[j][1];
      FF[i]  += f[j][2]*v[j][2];
      // velocity^2*mass part of Hoover drag
      DFF[i] += (v[j][0]*v[j][0] + v[j][1]*v[j][1] +
                 v[j][2]*v[j][2])*mass[type[j]]; 
    }
  }
  
  for (i = 0; i < ngridtot; i++)
    doubleTemp[i] = FF[i];
  MPI_Allreduce(doubleTemp,FF,ngridtot,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < ngridtot; i++)
    doubleTemp[i] = DFF[i];
  MPI_Allreduce(doubleTemp,DFF,ngridtot,MPI_DOUBLE,MPI_SUM,world);
  
  // force rescaling
  
  double f_factor;
  for (j = 0; j < nlocal; j++) {
    if (!(mask[j] & groupbit[1])) continue;
    i = grid[j];
    if (i >= 0 && DFF[i] > 0) {
      f_factor = FF[i]/DFF[i];
      f[j][0] -= v[j][0]*mass[type[j]]*f_factor;
      f[j][1] -= v[j][1]*mass[type[j]]*f_factor;
      f[j][2] -= v[j][2]*mass[type[j]]*f_factor;
    }
    if (i >= 0 && gridVecflag) {
      v[j][0] += gridVec[i][0];
      v[j][1] += gridVec[i][1];
      v[j][2] += gridVec[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- 
   grid ergodic loop with region matching
------------------------------------------------------------------------- */

void FixEM::compute_grid_voltage()
{
  int i,j,k,set;
  // grid region matching i number
  int XL,XR,YL,YR,ZL,ZR;
  double RKsum,matb,V;
  
  // PETSc variables
  Vec x,b,xpc;              /* solution, right hand side, pc vector*/
  Mat A;                    /* linear system matrix */
  KSP ksp;                  /* linear solver context */
  PC pc;                    /* preconditioner context */
  SampleShellPC  *shell;    /* user-defined preconditioner context */
  PetscInt Istart,Iend;
  PetscScalar *array;
  //PetscErrorCode ierr;
  
  MatCreate(PETSC_COMM_WORLD,&A);
  MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,ngridtot,ngridtot);
  MatSetType(A,MATAIJ);
  MatSetFromOptions(A);
  if (size > 1){
    MatMPIAIJSetPreallocation(A,7+ABCmax,PETSC_NULL,6+ABCmax,PETSC_NULL);
  } else {
    MatSeqAIJSetPreallocation(A,7+ABCmax,PETSC_NULL);
  }
  VecCreate(PETSC_COMM_WORLD,&b);
  VecSetSizes(b,PETSC_DECIDE,ngridtot);
  VecSetFromOptions(b);
  VecDuplicate(b,&x);
  VecDuplicate(b,&xpc);
  
  MatGetOwnershipRange(A,&Istart,&Iend);

  for (i = Istart; i < Iend; i++) {
    
    // set i into ixx,iyy,izz,iabc
    set = i;
    if (set < ngrideach[0]) {
      iabc = 0;
    } else if (set >= ngrideach[0] && set < ngrideach[0] + ngrideach[1]) {
      iabc = 1;
      set -= ngrideach[0];
    } else {
      iabc = 2;
      set -= ngrideach[0] + ngrideach[1];
    }
    ixx = set/(ngrids[iabc][1]*ngrids[iabc][2]);
    iyy = (set - ixx*ngrids[iabc][1]*ngrids[iabc][2])/ngrids[iabc][2];
    izz = set - ixx*ngrids[iabc][1]*ngrids[iabc][2] - iyy*ngrids[iabc][2];
    
    // reset RKsum, matb
    RKsum = matb = 0.;
    
    // SETUP GENERAL Y DIRECTIONS LEFT SIDE
    if (iyy > 0) {
      YL = get_grid_num(ixx,iyy-1,izz,iabc);
      V = -2./(eres[YL][1] + eres[i][1]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&YL,&V,INSERT_VALUES);
    } else if (periodic_yz && iyy == 0) {
      YL = get_grid_num(ixx,ngrids[iabc][1]-1,izz,iabc);
      V = -2./(eres[YL][1] + eres[i][1]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&YL,&V,INSERT_VALUES);
    }
    
    // SETUP GENERAL Y DIRECTIONS RIGHT SIDE
    if (iyy < ngrids[iabc][1]-1) {
      YR = get_grid_num(ixx,iyy+1,izz,iabc);
      V = -2./(eres[YR][1] + eres[i][1]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&YR,&V,INSERT_VALUES);
    } else if (periodic_yz && iyy == ngrids[iabc][1]-1) {
      YR = get_grid_num(ixx,0,izz,iabc);
      V = -2./(eres[YR][1] + eres[i][1]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&YR,&V,INSERT_VALUES);
	  }
    
    // SETUP GENERAL Z DIRECTIONS LEFT SIDE
    if (izz > 0) {
      ZL = get_grid_num(ixx,iyy,izz-1,iabc);
      V = -2./(eres[ZL][2] + eres[i][2]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&ZL,&V,INSERT_VALUES);
    } else if (periodic_yz && izz == 0) {
      ZL = get_grid_num(ixx,iyy,ngrids[iabc][2]-1,iabc);
      V = -2./(eres[ZL][2] + eres[i][2]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&ZL,&V,INSERT_VALUES);
    }
    
    // SETUP GENERAL Z DIRECTIONS RIGHT SIDE
    if (izz < ngrids[iabc][2]-1) {
      ZR = get_grid_num(ixx,iyy,izz+1,iabc);
      V = -2./(eres[ZR][2] + eres[i][2]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&ZR,&V,INSERT_VALUES);
    } else if (periodic_yz && izz == ngrids[iabc][2]-1) {
      ZR = get_grid_num(ixx,iyy,0,iabc);
      V = -2./(eres[ZR][2] + eres[i][2]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&ZR,&V,INSERT_VALUES);
    }
    
    /*  region matching
     ------------!--------!-----!-----!-----!--------!------------
                 !        !     !     !     !        !
                 !        !-----!-----!-----!        !
                 !        !     !     !     !        !
           ------!--------!-----!-----!-----!--------!------
                 !        !     !     !     !        !
                 !        !-----!-----!-----!        !
                 !        !     !     !     !        !
     ------------!--------!-----!-----!-----!--------!------------
     */
    
    if (ixx > 0) {
      // SETUP GENERAL X DIRECTIONS LEFT SIDE
      XL = get_grid_num(ixx-1,iyy,izz,iabc);
      V = -2./(eres[XL][0] + eres[i][0]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&XL,&V,INSERT_VALUES);
    } else {
      if (iabc == 0) {
        // left border, | [A]
        RKsum += 2./eres[i][0];
	  } else if (iabc == 1) {
        // region matching, A [B]
        iyr = ((iyy+1)%BOAY == 0) ? (iyy+1)/BOAY-1 : (iyy+1)/BOAY;
        izr = ((izz+1)%BOAZ == 0) ? (izz+1)/BOAZ-1 : (izz+1)/BOAZ;
        XL = get_grid_num(ngrids[0][0]-1,iyr,izr,0);
        V = -2./(ABCouple*eres[XL][0] + eres[i][0]);
        RKsum += -V;
        MatSetValues(A,1,&i,1,&XL,&V,INSERT_VALUES);
      }
    }
    
    if (ixx < ngrids[iabc][0]-1) {
      // SETUP GENERAL X DIRECTIONS RIGHT SIDE
      XR = get_grid_num(ixx+1,iyy,izz,iabc);
      V = -2./(eres[XR][0] + eres[i][0]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&XR,&V,INSERT_VALUES);
    } else {
      if (iabc == 2) {
        // right border, [C] |
        RKsum += 2./eres[i][0];
	  } else if (iabc == 1) {
        // region matching, [B] C
        iyr = ((iyy+1)%BOCY == 0) ? (iyy+1)/BOCY-1 : (iyy+1)/BOCY;
        izr = ((izz+1)%BOCZ == 0) ? (izz+1)/BOCZ-1 : (izz+1)/BOCZ;
        XR = get_grid_num(0,iyr,izr,2);
        V = -2./(BCCouple*eres[XR][0] + eres[i][0]);
        RKsum += -V;
        MatSetValues(A,1,&i,1,&XR,&V,INSERT_VALUES);
      }
    }
    
    if (iabc == 0) {
      if (ixx == 0) matb = fd_vlv/(eres[i][0]/2.0);
      if (ixx == ngrids[0][0]-1) {
        // region matching, [A] B
        for (j = 1; j <= ngrids[1][1]; j++) {
          iyr = (j%BOAY == 0) ? j/BOAY-1 : j/BOAY;
          if (iyr == iyy) {
            for (k = 1; k <= ngrids[1][2]; k++) {
              izr = (k%BOAZ == 0) ? k/BOAZ-1 : k/BOAZ;
              if (izr == izz) {
                XR= get_grid_num(0,j-1,k-1,1);
                V = -2./(eres[XR][0] + ABCouple*eres[i][0]);
                MatSetValues(A,1,&i,1,&XR,&V,INSERT_VALUES);
                RKsum += -V;
              }
            }
          }
        }
      }
    }
    
    if (iabc == 2) {
      if (ixx == ngrids[2][0]-1) matb += fd_vuv/(eres[i][0]/2.0);
      if (ixx == 0) {
        // region matching, B [C]
        for (j = 1; j <= ngrids[1][1]; j++) {
          iyr = (j%BOCY == 0) ? j/BOCY-1 : j/BOCY;
          if (iyr == iyy) {
            for (k = 1; k <= ngrids[1][2]; k++) {
              izr = (k%BOCZ == 0) ? k/BOCZ-1 : k/BOCZ;
              if (izr == izz) {
                XL = get_grid_num(ngrids[1][0]-1,j-1,k-1,1);
                V = -2./(eres[XL][0] + BCCouple*eres[i][0]);
                MatSetValues(A,1,&i,1,&XL,&V,INSERT_VALUES);
                RKsum += -V;
              }
            }
          }
        }
      }
    }
    
    VecSetValues(b,1,&i,&matb,INSERT_VALUES);
    VecSetValues(xpc,1,&i,&(gridVSAVE[i]),INSERT_VALUES);
    MatSetValues(A,1,&i,1,&i,&RKsum,INSERT_VALUES);
  }
  
  // assemble A, b
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
  VecAssemblyBegin(xpc);
  VecAssemblyEnd(xpc);
  
  // create ksp solver
  KSPCreate(PETSC_COMM_WORLD,&ksp);
  // setup ksp solver
  KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);
  // Bi-conjugate gradient method: KSPBICG
  // KSPSetType(ksp,KSPBICG);
  KSPGetPC(ksp,&pc);
  // PCSetType(pc,PCBJACOBI);
  // setup tolerance here
  KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
  
  /* (Required) Indicate to PETSc that we're using a "shell" preconditioner */
  PCSetType(pc,PCSHELL);
  
  /* (Optional) Create a context for the user-defined preconditioner; this
   context can be used to contain any application-specific data. */
  SampleShellPCCreate(&shell);
  
  /* (Required) Set the user-defined routine for applying the preconditioner */
  PCShellSetApply(pc,SampleShellPCApply);
  PCShellSetContext(pc,shell);
  
  /* (Optional) Set user-defined function to free objects used by custom preconditioner */
  PCShellSetDestroy(pc,SampleShellPCDestroy);
  
  /* (Optional) Set a name for the preconditioner, used for PCView() */
  // PCShellSetName(pc,"MyPreconditioner");
  
  /* (Optional) Do any setup required for the preconditioner */
  /* Note: This function could be set with PCShellSetSetUp and it would be called when necessary */
  SampleShellPCSetUp(pc,A,xpc);
  
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);
  // solve!
  KSPSolve(ksp,b,x);
  
  //  VecNorm(x,NORM_2,&norm);
  //  KSPGetIterationNumber(ksp,&its);
  //  if (rank == size-1) printf("voltage norm %f,its %d\n",norm,its);
  
  for (i = 0; i < ngridtot; i++)
    gridV[i] = doubleTemp[i] = 0.0;
  
  VecGetOwnershipRange(x,&Istart,&Iend);
  VecGetArray(x,&array);
  for (i = Istart; i < Iend; i++)
    doubleTemp[i] = array[i-Istart];
  VecRestoreArray(x,&array);
  MPI_Allreduce(doubleTemp,gridV,ngridtot,MPI_DOUBLE,MPI_SUM,world);
  
  for (i = 0; i < ngridtot; i++)
    gridVSAVE[i] = gridV[i];
  
  MatDestroy(&A);
  VecDestroy(&x);
  VecDestroy(&b);
  VecDestroy(&xpc);
  KSPDestroy(&ksp);
}

/* ---------------------------------------------------------------------- 
   grid ergodic loop with region matching
------------------------------------------------------------------------- */

void FixEM::compute_grid_temp()
{
  int j,k,set;
  int i;
  // grid region matching i number
  int XL,XR,YL,YR,ZL,ZR;
  double RKsum,matb,V;
  
  // PETSc variables
  Vec x,b,xpc;              /* solution, right hand side, pc vector*/
  Mat A;                    /* linear system matrix */
  KSP ksp;                  /* linear solver context */
  PC pc;                    /* preconditioner context */
  SampleShellPC  *shell;    /* user-defined preconditioner context */
  PetscInt Istart,Iend;
  PetscScalar *array;
  //PetscErrorCode ierr;
  
  MatCreate(PETSC_COMM_WORLD,&A);
  MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,ngridtot,ngridtot);
  MatSetType(A,MATAIJ);
  MatSetFromOptions(A);
  if (size > 1){
    MatMPIAIJSetPreallocation(A,7+ABCmax,PETSC_NULL,6+ABCmax,PETSC_NULL);
  } else {
    MatSeqAIJSetPreallocation(A,7+ABCmax,PETSC_NULL);
  }
  VecCreate(PETSC_COMM_WORLD,&b);
  VecSetSizes(b,PETSC_DECIDE,ngridtot);
  VecSetFromOptions(b);
  VecDuplicate(b,&x);
  VecDuplicate(b,&xpc);
  
  // get i range for this CPU
  MatGetOwnershipRange(A,&Istart,&Iend);
  
  for (i = Istart; i < Iend; i++) {
    
    // set i into ixx,iyy,izz,iabc
    set = i;
    if (set < ngrideach[0]) {
      iabc = 0;
    } else if (set >= ngrideach[0] && set < ngrideach[0] + ngrideach[1]) {
      iabc = 1;
      set -= ngrideach[0];
    } else {
      iabc = 2;
      set -= ngrideach[0] + ngrideach[1];
    }
    ixx = set/(ngrids[iabc][1]*ngrids[iabc][2]);
    iyy = (set - ixx*ngrids[iabc][1]*ngrids[iabc][2])/ngrids[iabc][2];
    izz = set - ixx*ngrids[iabc][1]*ngrids[iabc][2] - iyy*ngrids[iabc][2];
    
    // reset RKsum, matb
    RKsum = matb = 0.;
    
    // SETUP GENERAL Y DIRECTIONS LEFT SIDE
    if (iyy > 0) {
      YL = get_grid_num(ixx,iyy-1,izz,iabc);
      V = -2./(tres[YL][1] + tres[i][1]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&YL,&V,INSERT_VALUES);
    } else if (periodic_yz && iyy == 0) {
      YL = get_grid_num(ixx,ngrids[iabc][1]-1,izz,iabc);
      V = -2./(tres[YL][1] + tres[i][1]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&YL,&V,INSERT_VALUES);
    }
    
    // SETUP GENERAL Y DIRECTIONS RIGHT SIDE
    if (iyy < ngrids[iabc][1]-1) {
      YR = get_grid_num(ixx,iyy+1,izz,iabc);
      V = -2./(tres[YR][1] + tres[i][1]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&YR,&V,INSERT_VALUES);
    } else if (periodic_yz && iyy == ngrids[iabc][1]-1) {
      YR = get_grid_num(ixx,0,izz,iabc);
      V = -2./(tres[YR][1] + tres[i][1]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&YR,&V,INSERT_VALUES);
	  }
    
    // SETUP GENERAL Z DIRECTIONS LEFT SIDE
    if (izz > 0) {
      ZL = get_grid_num(ixx,iyy,izz-1,iabc);
      V = -2./(tres[ZL][2] + tres[i][2]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&ZL,&V,INSERT_VALUES);
    } else if (periodic_yz && izz == 0) {
      ZL = get_grid_num(ixx,iyy,ngrids[iabc][2]-1,iabc);
      V = -2./(tres[ZL][2] + tres[i][2]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&ZL,&V,INSERT_VALUES);
    }
    
    // SETUP GENERAL Z DIRECTIONS RIGHT SIDE
    if (izz < ngrids[iabc][2]-1) {
      ZR = get_grid_num(ixx,iyy,izz+1,iabc);
      V = -2./(tres[ZR][2] + tres[i][2]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&ZR,&V,INSERT_VALUES);
    } else if (periodic_yz && izz == ngrids[iabc][2]-1) {
      ZR = get_grid_num(ixx,iyy,0,iabc);
      V = -2./(tres[ZR][2] + tres[i][2]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&ZR,&V,INSERT_VALUES);
    }
    
    /*  region matching
     ------------!--------!-----!-----!-----!--------!------------
                 !        !     !     !     !        !
                 !        !-----!-----!-----!        !
                 !        !     !     !     !        !
           ------!--------!-----!-----!-----!--------!------
                 !        !     !     !     !        !
                 !        !-----!-----!-----!        !
                 !        !     !     !     !        !
     ------------!--------!-----!-----!-----!--------!------------
     */
    
    if (ixx > 0) {
      // SETUP GENERAL X DIRECTIONS LEFT SIDE
      XL = get_grid_num(ixx-1,iyy,izz,iabc);
      V = -2./(tres[XL][0] + tres[i][0]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&XL,&V,INSERT_VALUES);
    } else {
      if (iabc == 0) {
        // left border, | [A]
        RKsum += 2./tres[i][0];
	    } else if (iabc == 1) {
        // region matching, A [B]
        iyr = ((iyy+1)%BOAY == 0) ? (iyy+1)/BOAY-1 : (iyy+1)/BOAY;
        izr = ((izz+1)%BOAZ == 0) ? (izz+1)/BOAZ-1 : (izz+1)/BOAZ;
        XL = get_grid_num(ngrids[0][0]-1,iyr,izr,0);
        V = -2./(ABCouple*tres[XL][0] + tres[i][0]);
        RKsum += -V;
        MatSetValues(A,1,&i,1,&XL,&V,INSERT_VALUES);
      }
    }
    
    if (ixx < ngrids[iabc][0]-1) {
      // SETUP GENERAL X DIRECTIONS RIGHT SIDE
      XR = get_grid_num(ixx+1,iyy,izz,iabc);
      V = -2./(tres[XR][0] + tres[i][0]);
      RKsum += -V;
      MatSetValues(A,1,&i,1,&XR,&V,INSERT_VALUES);
    } else {
      if (iabc == 2) {
        // right border, [C] |
        RKsum += 2./tres[i][0];
      } else if (iabc == 1) {
        // region matching, [B] C
        iyr = ((iyy+1)%BOCY == 0) ? (iyy+1)/BOCY-1 : (iyy+1)/BOCY;
        izr = ((izz+1)%BOCZ == 0) ? (izz+1)/BOCZ-1 : (izz+1)/BOCZ;
        XR = get_grid_num(0,iyr,izr,2);
        V = -2./(BCCouple*tres[XR][0] + tres[i][0]);
        RKsum += -V;
        MatSetValues(A,1,&i,1,&XR,&V,INSERT_VALUES);
      }
    }
    
    if (iabc == 0) {
      if (ixx == 0) matb = t_request/(tres[i][0]/2.0);
      if (ixx == ngrids[0][0]-1) {
        // region matching, [A] B
        for (j = 1; j <= ngrids[1][1]; j++) {
          iyr = (j%BOAY == 0) ? j/BOAY-1 : j/BOAY;
          if (iyr == iyy) {
            for (k = 1; k <= ngrids[1][2]; k++) {
              izr = (k%BOAZ == 0) ? k/BOAZ-1 : k/BOAZ;
              if (izr == izz) {
                XR = get_grid_num(0,j-1,k-1,1);
                V = -2./(tres[XR][0] + ABCouple*tres[i][0]);
                MatSetValues(A,1,&i,1,&XR,&V,INSERT_VALUES);
                RKsum += -V;
              }
            }
          }
        }
      }
    }
    
    if (iabc == 2) {
      if (ixx == ngrids[2][0]-1) matb += t_request/(tres[i][0]/2.0);
      if (ixx == 0) {
        // region matching, B [C]
        for (j = 1; j <= ngrids[1][1]; j++) {
          iyr = (j%BOCY == 0) ? j/BOCY-1 : j/BOCY;
          if (iyr == iyy) {
            for (k = 1; k <= ngrids[1][2]; k++) {
              izr = (k%BOCZ == 0) ? k/BOCZ-1 : k/BOCZ;
              if (izr == izz) {
                XL = get_grid_num(ngrids[1][0]-1,j-1,k-1,1);
                V = -2./(tres[XL][0] + BCCouple*tres[i][0]);
                MatSetValues(A,1,&i,1,&XL,&V,INSERT_VALUES);
                RKsum += -V;
              }
            }
          }
        }
      }
    }
    
    matb += TgridRT[i]*gridCP[i]/(dtv*1e-12);
    RKsum += gridCP[i]/(dtv*1e-12);
    
    VecSetValues(b,1,&i,&matb,INSERT_VALUES);
    VecSetValues(xpc,1,&i,&(TgridREF[i]),INSERT_VALUES);
    MatSetValues(A,1,&i,1,&i,&RKsum,INSERT_VALUES);
  }
  
  // assemble A, b
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
  VecAssemblyBegin(xpc);
  VecAssemblyEnd(xpc);
  
  // create ksp solver
  KSPCreate(PETSC_COMM_WORLD,&ksp);
  // setup ksp solver
  KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);
  // Bi-conjugate gradient method: KSPBICG
  // KSPSetType(ksp,KSPBICG);
  KSPGetPC(ksp,&pc);
  // PCSetType(pc,PCBJACOBI);
  // setup tolerance here
  KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
  
  /* (Required) Indicate to PETSc that we're using a "shell" preconditioner */
  PCSetType(pc,PCSHELL);
  
  /* (Optional) Create a context for the user-defined preconditioner; this
     context can be used to contain any application-specific data. */
  SampleShellPCCreate(&shell);
  
  /* (Required) Set the user-defined routine for applying the preconditioner */
  PCShellSetApply(pc,SampleShellPCApply);
  PCShellSetContext(pc,shell);
  
  /* (Optional) Set user-defined function to free objects used by custom preconditioner */
  PCShellSetDestroy(pc,SampleShellPCDestroy);
  
  /* (Optional) Set a name for the preconditioner, used for PCView() */
  // PCShellSetName(pc,"MyPreconditioner");
  
  /* (Optional) Do any setup required for the preconditioner */
  /* Note: This function could be set with PCShellSetSetUp and it would be called when necessary */
  SampleShellPCSetUp(pc,A,xpc);
  
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);
  // solve!
  KSPSolve(ksp,b,x);
  
  //  VecNorm(x,NORM_2,&norm);
  //  KSPGetIterationNumber(ksp,&its);
  //  if (rank == size-1) printf("temperature norm %f,its %d\n",norm,its);
  
  for (i = 0; i < ngridtot; i++)
    TgridRT[i] = doubleTemp[i] = 0.0;
  
  VecGetOwnershipRange(x,&Istart,&Iend);
  VecGetArray(x,&array);
  for (i = Istart; i < Iend; i++)
    doubleTemp[i] = array[i-Istart];
  VecRestoreArray(x,&array);
  MPI_Allreduce(doubleTemp,TgridRT,ngridtot,MPI_DOUBLE,MPI_SUM,world);
  
  MatDestroy(&A);
  VecDestroy(&x);
  VecDestroy(&b);
  VecDestroy(&xpc);
  KSPDestroy(&ksp);
}

/* ---------------------------------------------------------------------- 
   grid ergodic loop with region matching
------------------------------------------------------------------------- */

void FixEM::compute_grid_current()
{
  int i,j,k,set;
  // grid region matching i number
  int XL,XR,YL,YR,ZL,ZR;
  double delta_V,delta_R;
  
  for (i = 0; i < ngridtot; i++) {
    
    // set i into ixx,iyy,izz,iabc
    set = i;
    if (set < ngrideach[0]) {
      iabc = 0;
    } else if (set >= ngrideach[0] && set < ngrideach[0] + ngrideach[1]) {
      iabc = 1;
      set -= ngrideach[0];
    } else {
      iabc = 2;
      set -= ngrideach[0] + ngrideach[1];
    }
    ixx = set/(ngrids[iabc][1]*ngrids[iabc][2]);
    iyy = (set - ixx*ngrids[iabc][1]*ngrids[iabc][2])/ngrids[iabc][2];
    izz = set - ixx*ngrids[iabc][1]*ngrids[iabc][2] - iyy*ngrids[iabc][2];
    
    // initialize gridI gridE	
    gridI[i][0] = gridE[i][0] = 0;
    gridI[i][1] = gridE[i][1] = 0;
    gridI[i][2] = gridE[i][2] = 0;
    
    // SETUP GENERAL Y DIRECTIONS LEFT SIDE
    if (iyy > 0) {
      YL = get_grid_num(ixx,iyy-1,izz,iabc);
      delta_R = (eres[i][1] + eres[YL][1])/2.;
      delta_V = gridV[i] - gridV[YL];
      gridI[i][1] += delta_V / delta_R;
      gridE[i][1] += delta_V / grid_size[iabc][1];
    } else if (periodic_yz && iyy == 0) {
      YL = get_grid_num(ixx,ngrids[iabc][1]-1,izz,iabc);
      delta_R = (eres[i][1] + eres[YL][1])/2.;
      delta_V = gridV[i] - gridV[YL];
      gridI[i][1] += delta_V / delta_R;
      gridE[i][1] += delta_V / grid_size[iabc][1];
    }
    
    // SETUP GENERAL Y DIRECTIONS RIGHT SIDE
    if (iyy < ngrids[iabc][1]-1) {
      YR = get_grid_num(ixx,iyy+1,izz,iabc);
      delta_R = (eres[i][1] + eres[YR][1])/2.;
      delta_V = gridV[YR] - gridV[i];
      gridI[i][1] += delta_V / delta_R;
      gridE[i][1] += delta_V / grid_size[iabc][1];
    } else if (periodic_yz && iyy == ngrids[iabc][1]-1) {
      YR = get_grid_num(ixx,0,izz,iabc);
      delta_R = (eres[i][1] + eres[YR][1])/2.;
      delta_V = gridV[YR] - gridV[i];
      gridI[i][1] += delta_V / delta_R;
      gridE[i][1] += delta_V / grid_size[iabc][1];
    }
    
    // SETUP GENERAL Z DIRECTIONS LEFT SIDE
    if (izz > 0) {
      ZL = get_grid_num(ixx,iyy,izz-1,iabc);
      delta_R = (eres[i][2] + eres[ZL][2])/2.;
      delta_V = gridV[i] - gridV[ZL];
      gridI[i][2] += delta_V / delta_R;
      gridE[i][2] += delta_V / grid_size[iabc][2];
    } else if (periodic_yz && izz == 0) {
      ZL = get_grid_num(ixx,iyy,ngrids[iabc][2]-1,iabc);
      delta_R = (eres[i][2] + eres[ZL][2])/2.;
      delta_V = gridV[i] - gridV[ZL];
      gridI[i][2] += delta_V / delta_R;
      gridE[i][2] += delta_V / grid_size[iabc][2];
    }
    
    // SETUP GENERAL Z DIRECTIONS RIGHT SIDE
    if (izz < ngrids[iabc][2]-1) {
      ZR = get_grid_num(ixx,iyy,izz+1,iabc);
      delta_R = (eres[i][2] + eres[ZR][2])/2.;
      delta_V = gridV[ZR] - gridV[i];
      gridI[i][2] += delta_V / delta_R;
      gridE[i][2] += delta_V / grid_size[iabc][2];
    } else if (periodic_yz && izz == ngrids[iabc][2]-1) {
      ZR = get_grid_num(ixx,iyy,0,iabc);
      delta_R = (eres[i][2] + eres[ZR][2])/2.;
      delta_V = gridV[ZR] - gridV[i];
      gridI[i][2] += delta_V / delta_R;
      gridE[i][2] += delta_V / grid_size[iabc][2];
    }
    
    /*  region matching in X direction
     ------------!--------!-----!-----!-----!--------!------------
                 !        !     !     !     !        !
                 !        !-----!-----!-----!        !
                 !        !     !     !     !        !
           ------!--------!-----!-----!-----!--------!------
                 !        !     !     !     !        !
                 !        !-----!-----!-----!        !
                 !        !     !     !     !        !
     ------------!--------!-----!-----!-----!--------!------------
     */
    
    if (ixx > 0) {
      // SETUP GENERAL X DIRECTIONS LEFT SIDE
      XL = get_grid_num(ixx-1,iyy,izz,iabc);
      delta_R = (eres[i][0] + eres[XL][0])/2.;
      delta_V = gridV[i] - gridV[XL];
      gridI[i][0] += delta_V / delta_R;
      gridE[i][0] += delta_V / grid_size[iabc][0];
    } else {
      if (iabc == 1) {
        // region matching, A [B]
        iyr = ((iyy+1)%BOAY == 0) ? (iyy+1)/BOAY-1 : (iyy+1)/BOAY;
        izr = ((izz+1)%BOAZ == 0) ? (izz+1)/BOAZ-1 : (izz+1)/BOAZ;
        XL = get_grid_num(ngrids[0][0]-1,iyr,izr,0);
        delta_R = (ABCouple*eres[XL][0] + eres[i][0])/2.;
        delta_V = gridV[i] - gridV[XL];
        gridI[i][0] += delta_V / delta_R;
        gridE[i][0] += delta_V / grid_size[iabc][0];
      }
    }
    
    if (ixx < ngrids[iabc][0]-1) {
      // SETUP GENERAL X DIRECTIONS RIGHT SIDE
      XR = get_grid_num(ixx+1,iyy,izz,iabc);
      delta_R = (eres[i][0] + eres[XR][0])/2.;
      delta_V = gridV[XR] - gridV[i];
      gridI[i][0] += delta_V / delta_R;
      gridE[i][0] += delta_V / grid_size[iabc][0];
    } else {
      if (iabc == 1) {
        // region matching, [B] C
        iyr = ((iyy+1)%BOCY == 0) ? (iyy+1)/BOCY-1 : (iyy+1)/BOCY;
        izr = ((izz+1)%BOCZ == 0) ? (izz+1)/BOCZ-1 : (izz+1)/BOCZ;
        XR = get_grid_num(0,iyr,izr,2);
        delta_R = (BCCouple*eres[XR][0] + eres[i][0])/2.0;
        delta_V = gridV[XR] - gridV[i];
        gridI[i][0] += delta_V / delta_R;
        gridE[i][0] += delta_V / grid_size[iabc][0];
      }
    }
    
    if (iabc == 0) {
      if (ixx == 0) gridI[i][0] *= 2;
      if (ixx == ngrids[0][0]-1) {
        // region matching, [A] B
        XR = -1;
        for (j = 1; j <= ngrids[1][1]; j++) {
          iyr = (j%BOAY == 0) ? j/BOAY-1 : j/BOAY;
          if (iyr == iyy) {
            for (k = 1; k <= ngrids[1][2]; k++) {
              izr = (k%BOAZ == 0) ? k/BOAZ-1 : k/BOAZ;
              if (izr == izz) {
                XR = get_grid_num(0,j-1,k-1,1);
                delta_R = (ABCouple*eres[i][0] + eres[XR][0])/2.0;
                delta_V = gridV[XR] - gridV[i];
                gridI[i][0] += delta_V / delta_R;
                gridE[i][0] += delta_V / grid_size[iabc][0];
              }
            }
          }
        }
      }
    }
    
    if (iabc == 2) {
      if (ixx == ngrids[2][0]-1) gridI [i][0] *= 2;
      if (ixx == 0) {
        // region matching, B [C]
        for (j = 1; j <= ngrids[1][1]; j++) {
          iyr = (j%BOCY == 0) ? j/BOCY-1 : j/BOCY;
          if (iyr == iyy) {
            for (k = 1; k <= ngrids[1][2]; k++) {
              izr = (k%BOCZ == 0) ? k/BOCZ-1 : k/BOCZ;
              if (izr == izz) {
                XL = get_grid_num(ngrids[1][0]-1,j-1,k-1,1);
                delta_R = (BCCouple*eres[i][0] + eres[XL][0])/2.;
                delta_V = gridV[i] - gridV[XL];
                gridI[i][0] += delta_V / delta_R;
                gridE[i][0] += delta_V / grid_size[iabc][0];
              }
            }
          }
        }
      }
    }
    
    // correct overlap
    gridI[i][0] /= 2.;
    gridE[i][0] /= 2.;
    gridI[i][1] /= 2.;
    gridE[i][1] /= 2.;
    gridI[i][2] /= 2.;
    gridE[i][2] /= 2.;
  }
}

/* ----------------------------------------------------------------------
   ELECTROMEAM Plugin called inside LAMMPS MD loops
------------------------------------------------------------------------- */

void FixEM::initial_integrate(int vflag)
{
  int mdstep = update->ntimestep;
  if (mdstep < delay) return;
  if (mdstep % boost != 0) return; 
  if (!fdonly) {
    check_atoms();
    check_atc();
    set_grid_velocity();
    set_grid_ke();
  }
  set_grid_temp();
  set_grid_eres();
  set_grid_tres();
  set_grid_cp();
  compute_grid_voltage();
  compute_grid_current();
  jouleheat();
  compute_grid_temp();
  if (!fdonly)
    thermostat();
  if (rank == 0) 
    dataout(mdstep);
  update_grid_temp();
}

/* ---------------------------------------------------------------------- */

void FixEM::final_integrate()
{

}

/* ---------------------------------------------------------------------- */

void FixEM::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  
  // fix em always take outermost level
  if (ilevel == nlevels_respa-1) initial_integrate(vflag);
}

/* ---------------------------------------------------------------------- */

void FixEM::final_integrate_respa(int ilevel, int iloop)
{
  
  // fix em always take outermost level
  if (ilevel == nlevels_respa-1) final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixEM::reset_dt()
{
  dtv = update->dt * boost;
  dtf = update->dt * force->ftm2v * boost;
}

/* ---------------------------------------------------------------------- */

void FixEM::dataout(int mdstep)
{  
  // fort.71
  if (mdstep%fort71 == 0) {
    FILE *fpfort71;
    int sum_gridNATOMS;
    int i,ave;
    double T_gridRT,T_gridREF,T_gridSAVE,maxT;
    double centerX;
    
    maxT = 0.;
    
    fpfort71 = fopen("fort.71","a");
    fprintf(fpfort71,"Step: %010d\n",mdstep);
    
    for (iabc = 0; iabc < 3; iabc++) {
      
      if (iabc == 0) 
        centerX = grid_lb[0]-grid_size[0][0]/2;
      else if (iabc == 1) 
        centerX = cxlower-grid_size[1][0]/2;
      else 
        centerX = cxupper-grid_size[2][0]/2;
      
      for (ixx = 0; ixx < ngrids[iabc][0]; ixx++) {
        T_gridRT = T_gridREF = T_gridSAVE = 0;
        sum_gridNATOMS = 0;
        ave = 0;
        centerX += grid_size[iabc][0];
        for (iyy = 0; iyy < ngrids[iabc][1]; iyy++) {
          for (izz = 0; izz < ngrids[iabc][2]; izz++) {
            i = get_grid_num(ixx,iyy,izz,iabc);
            if (ngridcnt[i] == 0 && iabc == 1) continue;
            T_gridRT += TgridRT[i];
            T_gridREF += TgridREF[i];
            T_gridSAVE += TgridSAVE[i];
            sum_gridNATOMS += ngridcnt[i];
            ave++;
          }
        }
        if (ave == 0) {
          for (iyy = 0; iyy < ngrids[iabc][1]; iyy++) {
            for (izz = 0; izz < ngrids[iabc][2]; izz++) {
              i = get_grid_num(ixx,iyy,izz,iabc);
              T_gridRT += TgridRT[i];
              T_gridREF += TgridREF[i];
              T_gridSAVE += TgridSAVE[i];
              ave++;
            }
          }
        }
        T_gridRT /= ave;
        T_gridREF /= ave;
        T_gridSAVE /= ave;      
        i = get_grid_num(ixx,ngrids[iabc][1]/2,ngrids[iabc][2]/2,iabc);
        
        fprintf(fpfort71,"%+15.6e %15.6e %15.6e %15.6e %08d\n", \
                centerX,T_gridRT,T_gridREF,T_gridSAVE,sum_gridNATOMS);      
      }
    }
    fclose(fpfort71);
  }
  
  // fddump
  if (mdstep%fddump == 0) {
    FILE *fpfddump;
    int i;
    char fddumpname[64];
    // properties to be dumped
    sprintf(fddumpname,"FDDUMP-LAMMPS-EM_%010d.txt", mdstep);
    fpfddump = fopen(fddumpname,"w");
    fprintf(fpfddump, \
            "//----------------- TgridRT,TgridREF,TgirdSAVE ---------------\n");
    for (i = 0; i < ngridtot; i++) {
      fprintf(fpfddump, "%20.11e %20.11e %20.11e %06d\n", \
               TgridRT[i],TgridREF[i],TgridSAVE[i],ngridcnt[i]);
    }
    fprintf(fpfddump, \
            "//---------------------- eres[3] ------------------\n");
    for (i = 0; i < ngridtot; i++) {
      fprintf(fpfddump, "%20.11e %20.11e %20.11e %02d\n", \
              eres[i][0], eres[i][1],eres[i][2],atc[i]);
    }
    fprintf(fpfddump, \
            "//----------------------- tres[3] ------------------\n");
    for (i = 0; i < ngridtot; i++) {
      fprintf(fpfddump, "%20.11e %20.11e %20.11e\n", \
              tres[i][0],tres[i][1],tres[i][2]);
    }
    fprintf(fpfddump, \
            "//-------------------------- gridI[3] ------------------------\n");
    for (i = 0; i < ngridtot; i++) {
      fprintf(fpfddump, "%20.11e %20.11e %20.11e\n", \
              gridI[i][0],gridI[i][1],gridI[i][2]);
    }
    fprintf(fpfddump, \
            "//-------------------- gridV gridKE gridCP -------------------\n");
    for (i = 0; i < ngridtot; i++) {
      fprintf(fpfddump, "%20.11e %20.11e %20.11e\n", \
              gridV[i],gridKE[i],gridCP[i]);
    }
    fclose(fpfddump);
  }
  
  for (int i = 0; i < 16; i++) {
    if (para_select[i]) paraout(i,mdstep);
  }
}

/* ---------------------------------------------------------------------- */

void FixEM::paraout(int vflag, int mdstep)
{
  if (mdstep%para_dump[vflag] != 0) return;
  FILE *fppara;
  int i,j,k;
  int Ntotal = (ngrids[0][0]+1)*(ngrids[0][1]+1)*(ngrids[0][2]+1) + 
               (ngrids[1][0]+1)*(ngrids[1][1]+1)*(ngrids[1][2]+1) + 
               (ngrids[2][0]+1)*(ngrids[2][1]+1)*(ngrids[2][2]+1);
  double paradata;
  char paraname[64];
  switch (vflag) {
    // TgridRT
    case 0:
      sprintf(paraname, "para_TgridRT_%010d.vtk", mdstep);
      break;
    // eres
    case 1:
      sprintf(paraname, "para_eres_%010d.vtk", mdstep);
      break;
    // tres
    case 2:
      sprintf(paraname, "para_tres_%010d.vtk", mdstep);
      break;
    // gridI
    case 3:
      sprintf(paraname, "para_i_%010d.vtk", mdstep);
      break;
    // gridV
    case 4:
      sprintf(paraname, "para_v_%010d.vtk", mdstep);
      break;
    // gridKE
    case 5:
      sprintf(paraname, "para_ke_%010d.vtk", mdstep);
      break;
    // gridcp
    case 6:
      sprintf(paraname, "para_cp_%010d.vtk", mdstep);
      break;
    default:
      break;
  }
  fppara = fopen(paraname, "w");
  fprintf(fppara, "# vtk DataFile Version 2.0\n");
  fprintf(fppara, "Solution: electromeam mesh, timestep %d\n",mdstep);
	fprintf(fppara, "ASCII\n\n");	
	fprintf(fppara, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fppara, "POINTS %d float\n", Ntotal);
  for (i = 0; i <= ngrids[0][0]; i++) {
		for (j = 0; j <= ngrids[0][1]; j++) {
			for (k = 0; k <= ngrids[0][2]; k++) {
        fprintf(fppara, "%17.10f%17.10f%17.10f\n", 
                grid_lb[0]+i*grid_size[0][0], 
                grid_lb[1]+j*grid_size[0][1], 
                grid_lb[2]+k*grid_size[0][2]);
			}
		}
	}
	for (i = 0; i <= ngrids[1][0]; i++) {
		for (j = 0; j <= ngrids[1][1]; j++) {
			for (k = 0; k <= ngrids[1][2]; k++) {
        fprintf(fppara, "%17.10f%17.10f%17.10f\n", 
                cxlower+i*grid_size[1][0], 
                grid_lb[1]+j*grid_size[1][1], 
                grid_lb[2]+k*grid_size[1][2]);
			}
		}
	}
	for (i = 0; i <= ngrids[2][0]; i++) {
		for (j = 0; j <= ngrids[2][1]; j++) {
			for (k = 0; k <= ngrids[2][2]; k++) {
        fprintf(fppara, "%17.10f%17.10f%17.10f\n", 
                cxupper+i*grid_size[2][0], 
                grid_lb[1]+j*grid_size[2][1], 
                grid_lb[2]+k*grid_size[2][2]);
			}
		}
	}
	fprintf(fppara, "CELLS %10d%10d\n",ngridtot,9*ngridtot);
  //       3 --- 2
  //      /|    /|
  //     / 0 --/ 1     y
  //    7 --- 6 /      |
  //    |     |/       |
  //    4 --- 5         ---> x
  //                  /
  //                 /
  //                z
	int icell = 0;
  
	for (i = 0; i <= ngrids[0][0]; i++) {
		for (j = 0; j <= ngrids[0][1]; j++) {
			for (k = 0; k <= ngrids[0][2]; k++) {
				if (!(i == ngrids[0][0] || j == ngrids[0][1] || k == ngrids[0][2])) {
          fprintf(fppara, "8%8d%8d%8d%8d%8d%8d%8d%8d\n", 
                  icell,  // 0
                  icell + (ngrids[0][1]+1)*(ngrids[0][2]+1),  // 1
					        icell + (ngrids[0][1]+2)*(ngrids[0][2]+1),  // 2
                  icell + ngrids[0][2] + 1,  // 3
					        icell + 1,  // 4
					        icell + (ngrids[0][1]+1)*(ngrids[0][2]+1) + 1,  // 5
					        icell + (ngrids[0][1]+2)*(ngrids[0][2]+1) + 1,  // 6
                  icell + ngrids[0][2] + 2); // 7
				}
				icell++;
			}
		}
	}
	for (i = 0; i <= ngrids[1][0]; i++) {
		for (j = 0; j <= ngrids[1][1]; j++) {
			for (k = 0; k <= ngrids[1][2]; k++) {
				if (!(i == ngrids[1][0] || j == ngrids[1][1] || k == ngrids[1][2])) {
					fprintf(fppara, "8%8d%8d%8d%8d%8d%8d%8d%8d\n", 
					        icell,  // 0
					        icell + (ngrids[1][1]+1)*(ngrids[1][2]+1),  // 1
					        icell + (ngrids[1][1]+2)*(ngrids[1][2]+1),  // 2
					        icell + ngrids[1][2] + 1,  // 3
                  icell + 1,  // 4
					        icell + (ngrids[1][1]+1)*(ngrids[1][2]+1) + 1,  // 5
					        icell + (ngrids[1][1]+2)*(ngrids[1][2]+1) + 1,  // 6
					        icell + ngrids[1][2] + 2); // 7
				}
				icell++;
			}
		}
	}
	for (i = 0; i <= ngrids[2][0]; i++) {
		for (j = 0; j <= ngrids[2][1]; j++) {
			for (k = 0; k <= ngrids[2][2]; k++) {
				if (!(i == ngrids[2][0] || j == ngrids[2][1] || k == ngrids[2][2])) {
					fprintf(fppara, "8%8d%8d%8d%8d%8d%8d%8d%8d\n", 
					        icell,  // 0
					        icell + (ngrids[2][1]+1)*(ngrids[2][2]+1),  // 1
					        icell + (ngrids[2][1]+2)*(ngrids[2][2]+1),  // 2
					        icell + ngrids[2][2] + 1,  // 3
					        icell + 1,  // 4
					        icell + (ngrids[2][1]+1)*(ngrids[2][2]+1) + 1,  // 5
					        icell + (ngrids[2][1]+2)*(ngrids[2][2]+1) + 1,  // 6
					        icell + ngrids[2][2] + 2); // 7
				}
				icell++;
			}
		}
	}
	fprintf(fppara, "CELL_TYPES%10d\n", ngridtot);
	for (i = 0; i < ngridtot; i++)
		fprintf(fppara, "12\n");
	fprintf(fppara, "POINT_DATA%10d\n", Ntotal);
	fprintf(fppara, "SCALARS pressure float 1\n");
	fprintf(fppara, "LOOKUP_TABLE default\n");
  for (i = 0; i <= ngrids[0][0]; i++) {
		for (j = 0; j <= ngrids[0][1]; j++) {
			for (k = 0; k <= ngrids[0][2]; k++) {
				int ii = i; int jj = j; int kk = k;
				if (ii == ngrids[0][0]) ii--;
				if (jj == ngrids[0][1]) jj--;
				if (kk == ngrids[0][2]) kk--;
        int gridnum = ii*ngrids[0][1]*ngrids[0][2] + jj*ngrids[0][2] + kk;
        switch (vflag) {
          case 0:
            paradata = TgridRT[gridnum];
            break;
          case 1:
            paradata = sqrt(eres[gridnum][0]*eres[gridnum][0] + 
                            eres[gridnum][1]*eres[gridnum][1] + 
                            eres[gridnum][2]*eres[gridnum][2]);
            paradata = paradata > 100000 ? 1 : paradata;
            break;
          case 2:
            paradata = sqrt(tres[gridnum][0]*tres[gridnum][0] + 
                            tres[gridnum][1]*tres[gridnum][1] + 
                            tres[gridnum][2]*tres[gridnum][2]);
            break;
          case 3:
            paradata = sqrt(gridI[gridnum][0]*gridI[gridnum][0] + 
                            gridI[gridnum][1]*gridI[gridnum][1] + 
                            gridI[gridnum][2]*gridI[gridnum][2]);
            break;
          case 4:
            paradata = gridV[gridnum];
            break;
          case 5:
            paradata = gridKE[gridnum];
            break;
          case 6:
            paradata = gridCP[gridnum];
            break;
          default:
            break;
        }
        fprintf(fppara, "%20.11e\n", paradata);
			}
		}
	}
	for (i = 0; i <= ngrids[1][0]; i++) {
		for (j = 0; j <= ngrids[1][1]; j++) {
			for (k = 0; k <= ngrids[1][2]; k++) {
				int ii = i; int jj = j; int kk = k;
				if (ii == ngrids[1][0]) ii--;
				if (jj == ngrids[1][1]) jj--;
				if (kk == ngrids[1][2]) kk--;
        int gridnum = ii*ngrids[1][1]*ngrids[1][2] + jj*ngrids[1][2] + kk + ngrideach[0];
        switch (vflag) {
          case 0:
            paradata = TgridRT[gridnum];
            break;
          case 1:
            paradata = sqrt(eres[gridnum][0]*eres[gridnum][0] + 
                            eres[gridnum][1]*eres[gridnum][1] + 
                            eres[gridnum][2]*eres[gridnum][2]);
            paradata = paradata > 100000 ? 1 : paradata;
            break;
          case 2:
            paradata = sqrt(tres[gridnum][0]*tres[gridnum][0] + 
                            tres[gridnum][1]*tres[gridnum][1] + 
                            tres[gridnum][2]*tres[gridnum][2]);
            break;
          case 3:
            paradata = sqrt(gridI[gridnum][0]*gridI[gridnum][0] + 
                            gridI[gridnum][1]*gridI[gridnum][1] + 
                            gridI[gridnum][2]*gridI[gridnum][2]);
            break;
          case 4:
            paradata = gridV[gridnum];
            break;
          case 5:
            paradata = gridKE[gridnum];
            break;
          case 6:
            paradata = gridCP[gridnum];
            break;
          default:
            break;
        }
        fprintf(fppara, "%20.11e\n", paradata);
			}
		}
	}
	for (i = 0; i <= ngrids[2][0]; i++) {
		for (j = 0; j <= ngrids[2][1]; j++) {
			for (k = 0; k <= ngrids[2][2]; k++) {
				int ii = i; int jj = j; int kk = k;
				if (ii == ngrids[2][0]) ii--;
				if (jj == ngrids[2][1]) jj--;
				if (kk == ngrids[2][2]) kk--;
        int gridnum = ii*ngrids[2][1]*ngrids[2][2] + jj*ngrids[2][2] + kk + ngrideach[0] + ngrideach[1];
        switch (vflag) {
          case 0:
            paradata = TgridRT[gridnum];
            break;
          case 1:
            paradata = sqrt(eres[gridnum][0]*eres[gridnum][0] + 
                            eres[gridnum][1]*eres[gridnum][1] + 
                            eres[gridnum][2]*eres[gridnum][2]);
            paradata = paradata > 100000 ? 1 : paradata;
            break;
          case 2:
            paradata = sqrt(tres[gridnum][0]*tres[gridnum][0] + 
                            tres[gridnum][1]*tres[gridnum][1] + 
                            tres[gridnum][2]*tres[gridnum][2]);
            break;
          case 3:
            paradata = sqrt(gridI[gridnum][0]*gridI[gridnum][0] + 
                            gridI[gridnum][1]*gridI[gridnum][1] + 
                            gridI[gridnum][2]*gridI[gridnum][2]);
            break;
          case 4:
            paradata = gridV[gridnum];
            break;
          case 5:
            paradata = gridKE[gridnum];
            break;
          case 6:
            paradata = gridCP[gridnum];
            break;
          default:
            break;
        }
        fprintf(fppara, "%20.11e\n", paradata);
			}
		}
	}

  
  
  
  fclose(fppara);  
}

/* ---------------------------------------------------------------------- */
