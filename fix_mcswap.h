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

FixStyle(mcswap,FixMCSWAP)

#else

#ifndef LMP_FIX_MCSWAP_H
#define LMP_FIX_MCSWAP_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixMCSWAP : public Fix {
 public:
   
  class NeighList *list;
  
  FixMCSWAP(class LAMMPS *, int, char **);
  ~FixMCSWAP();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void update_atom_list();
  void pre_exchange();
  void attempt_move();
  void attempt_swap();
  double compute_vector(int);
  double memory_usage();
  void write_restart(FILE *);
  void restart(char *);

 private:
   
  void attempt_postswap();
  double denergy(int, int, double *);

  int rank,size;
  int ntype,nevery,seed;
  int ncycles,nexchanges,nmcmoves;
  int ngas;           // # of gas molecules (or atoms) on all procs 
  int ngas_local;     // # of gas molecules (or atoms) on this proc 
  int ngas_before;    // # of gas molecules (or atoms) on procs < this proc
  int itype1,itype2;
  int ntype1,ntype2;
  int ntype1_local,ntype2_local;
  int ntype1_before,ntype2_before;
  int nswap,iswap;
  int npostswap;

  double nmove_attempts;
  double nmove_successes;
  double npostswap_attempts;
  double npostswap_successes;
  double nswap_attempts;
  double nswap_successes;
  
  double *balance_all,*balance_type1,*balance_type2;
  int *ngas_local_all,*ntype1_local_all,*ntype2_local_all;

  int nmax;
  double cutoff;
  double reservoir_temperature;
  double chemical_potential;
  double displace;
  double beta,zz,sigma,volume;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double mcxlo,mcxhi;
  double *sublo,*subhi;
  int *local_gas_list;
  int *local_type1_list,*local_type2_list;
  
  int scn;
  double *sc1x,*sc1y,*sc1z;
  double *sc2x,*sc2y,*sc2z;
  double *svc1x,*svc1y,*svc1z;
  double *svc2x,*svc2y,*svc2z;
  double **cutsq;
  
  int *sc1tag,*sc2tag;
  int *tag1,*tag2;
  int *dtag1,*dtag2,*dtag1_all,*dtag2_all;
  int *tag1list,*tag2list,*taglist;
  double *coordx1,*coordy1,*coordz1;
  double *coordx2,*coordy2,*coordz2;
  double *velx1,*vely1,*velz1;
  double *velx2,*vely2,*velz2;  
  double *denergy1_all,*denergy2_all;

  class Pair *pair;
  class RanPark *random_equal;
  class RanPark *random_unequal;
  
  void options(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid atom type in fix GCMC command

The atom type specified in the GCMC command does not exist.

E: Cannot do GCMC on atoms in atom_modify first group

This is a restriction due to the way atoms are organized in a list to
enable the atom_modify first command.

W: Fix GCMC may delete atom with non-zero molecule ID

This is probably an error, since you should not delete only one atom
of a molecule. The GCMC molecule exchange feature does not yet work.

E: Fix GCMC molecule command requires atom attribute molecule

Should not choose the GCMC molecule feature if no molecules are being
simulated. The general molecule flag is off, but GCMC's molecule flag
is on.

E: Fix GCMC molecule feature does not yet work

Fix GCMC cannot (yet) be used to exchange molecules, only atoms.

E: Fix GCMC incompatible with given pair_style

Some pair_styles do not provide single-atom energies, which are needed
by fix GCMC.

*/
