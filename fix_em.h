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

#ifdef FIX_CLASS

FixStyle(em,FixEM)

#else

#ifndef LMP_FIX_EM_H
#define LMP_FIX_EM_H

#include "fix.h"
#include "petsc.h"

namespace LAMMPS_NS {
  
  class FixEM : public Fix {
  public:
    int igroup,groupbit[8];
    int delay; // run MD only for this amount of steps
    int fort71,fddump;
    bool para_select[16];
    bool fdonly;
    bool restarted;
    int para_dump[16];
    int *intTemp1, *intTemp2;
    double *doubleTemp;

    int *grid,gridmax;
    int *atc; // define grid property of atomistic or continuum
    int *ngridcnt;  // number of atoms IN EACH Grid BOX
    int *ngridcnt_thermo; // number of atoms contribute to grid temperature
    bool periodic_yz; // set periodic boundary conditions for Alice Loop
    bool gridVecflag;
    double **eres;  // grid electric resistivity, global
    double **tres;  // grid thermal resistivity , global    
    double **gridI; // grid electric current
    double **gridE; // grid E field
    double **gridVec;
    double *TgridSAVE; // grid temperature from last time step
    double *TgridREF;  // grid temperature from MD, current step
    double *TgridRT;   // grid temperature, most up-to-date
    double *pa;        // percentage of alloy in each grid, based on atomic info.
    double *gridV,*gridKE,*gridCP;
    double *gridVSAVE; // gridV from last time step
    double *FF,*DFF;
    
    // atomic mass unit to kg, the same as ELECTROMEAM
    double amu2kg; 
        
    FixEM(class LAMMPS *, int, char **);
    virtual ~FixEM();
    int setmask();
    virtual void init();
    virtual void setup(int);
    virtual void initial_integrate(int);
    virtual void final_integrate();
    /** Fake initial_integrate_fd */
    virtual void initial_integrate_respa(int, int, int);
    /** Fake final_integrate_fd */
    virtual void final_integrate_respa(int, int);
    virtual void reset_dt();
    
    void write_restart(FILE *);
    virtual int pack_restart_data(double *); // pack restart data
    virtual void restart(char *);
    /** modify em parameters (parser) */
    int modify_param(int narg, char** arg);
    
    // electromeam subroutines
    /** get grid number */
    int get_grid_num(int, int, int, int);
    /** count num. atoms for each grid, calculate percent alloy */
    void check_atoms();
    /** check atc */
    void check_atc();
    /** calculate translational velocity */
    void set_grid_velocity();
    /** kinetic energy of atoms */
    void set_grid_ke();
    /** temp before FD, based on ke, tgridsave, etc. */
    void set_grid_temp();
    /** update temp after thermostat */
    void update_grid_temp();
    /** electric resistivity */
    void set_grid_eres();
    /** thermal conductivity */
    void set_grid_tres();
    /** heat capacity */
    void set_grid_cp();
    /** increase grid temp via joule heating */
    void jouleheat();
    /** Noose thermostat */
    void thermostat();
    /** PETSc subroutine: calculate grid voltage */
    void compute_grid_voltage();
    /** calculate current based on voltage */
    void compute_grid_current();
    /** PETSc subroutine: calculate grid temp */
    void compute_grid_temp();
    /** dump FD data files */
    void dataout(int);
    /** dump paraview data files */
    void paraout(int,int);
    
  protected:
    int rank,size;
    int nlevels_respa;
    double dtv,dtf;
    double *step_respa;
    /** total number of grids; */
    int ngridtot,ngridtot_restart;
    /** INPUT.D :: grid number, ngrids[region number][x,y,z] */
    int ngrids[3][3];
    /** total number of grids for each region, ngrideach[region number] */
    int ngrideach[3];
    /** grid size, grid_size[region number][x,y,z] */
    double grid_size[3][3];
    /** grid resistor factor, 
      * for calculating thermal and electrical resistivity, 
      * grid_Rfactor[region number][x,y,z] */
    double grid_Rfactor[3][3];
    /** grid volume, grid_volume[region number] */
    double grid_volume[3];
    /** INPUT.D :: the upper boundarys for over all regions */
    double grid_ub[3];
    /** INPUT.D :: the lower boundarys for over all regions */
    double grid_lb[3];
    /** INPUT.D :: continuum xbounds 1 */
    double cxlower;
    /** INPUT.D :: continuum xbounds 2 */
    double cxupper;
    int ctypel,ctypeu;
    int atype1,atype2;
    
    int boost;
    
    /** use get_gridcoord to update these coordinate numbers */
    int ixx,iyy,izz;
    int iyr,izr;
    /** identify the region   */
    int iabc;
    int BOAY,BOAZ,BOCY,BOCZ;
    int ABCouple,BCCouple,ABCmax;
    
    double t_request;    
    double FDdens[2];
    double FDrho[2],FDrhoT[2],FDalpharho[2];
    double MASS_VAL[2];
    double FDcp[2],FDtc[2];
    double fd_vuv,fd_vlv;
    
    virtual int size_restart_global();

  private:		
    // custom indenter
    int iflag,istyle;
    char *xstr,*ystr,*zstr,*rstr,*pstr;
    char *astr,*bstr,*cstr;
    int xvar,yvar,zvar,rvar,pvar;
    int avar,bvar,cvar;
    double xvalue,yvalue,zvalue,rvalue,pvalue;
    double avalue,bvalue,cvalue;
    double ctr[3],radiusr;

    // dynamic ATC boundary of BC region
    int bcflag;
    double BC;
  };
  
}

#endif
#endif
