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

#ifdef PAIR_CLASS

PairStyle(ctip,PairCtip)

#else

#ifndef LMP_PAIR_CTIP_H
#define LMP_PAIR_CTIP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCtip : public Pair {
 public:
  PairCtip(class LAMMPS *);
  virtual ~PairCtip();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double single(int, int, int, int, double, double, double, double& fforce) { fforce = 0.0; return 0.0; }
  double singleembedding(int, int, double *) {return 0.0;}
  double memory_usage();
  
  virtual double yasu_char(double *, int &);

 protected:
  struct Param {
    int type;
    double q_min,q_max;
    double chi,J,xi,Z;
  };
  
  double alf;
  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  int ***elem2param;            // mapping from element triplets to parameters
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  double precision;
  Param *params;                // parameter set

  int nmax;
  double *qf;

  double *esm, **fafb, **dfafb, **ddfafb, **phin, **dphin, **erpaw;
  double **afb, **dafb, **ddafb;
  double *charge;
  int **intype, *typeno;
  int *NCo, cor_flag, cuo_flag, cuo_flag1, cuo_flag2;
  double **bbij;

  void allocate();
  void setup();

  double self(int, double, double);
  void sm_table();
  void potal_calc(double &, double &, double &);
  
  void tri_point(double, int &, int &, int &, double &, double &, double &);
  void direct(int,int,int,int,int,double,double,double,double,double,double,
        double,double,double,double,double,double &,double &);  
  
  double qfo_self(int, double, double);
  void qfo_direct (int, int, int, int, int, double, double, double, double,
        double, double, double, double, double, double &, double &);

  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  int pack_comm(int , int *, double *, int, int *);
  void unpack_comm(int , int , double *);

  // short range neighbor list

  void add_pages(int howmany = 1);
  int maxpage, pgsize, oneatom, **pages;
  double cutmin;

  // vector functions, inline for efficiency

  inline double vec3_dot(const double x[3], const double y[3]) const {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
  }

  inline void vec3_add(const double x[3], const double y[3],
                       double * const z) const {
    z[0] = x[0]+y[0];  z[1] = x[1]+y[1];  z[2] = x[2]+y[2];
  }

  inline void vec3_scale(const double k, const double x[3],
                         double y[3]) const {
    y[0] = k*x[0];  y[1] = k*x[1];  y[2] = k*x[2];
  }

  inline void vec3_scaleadd(const double k, const double x[3],
                            const double y[3], double * const z) const {
    z[0] = k*x[0]+y[0];
    z[1] = k*x[1]+y[1];
    z[2] = k*x[2]+y[2];
  }
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style COMB requires atom IDs

This is a requirement to use the AIREBO potential.

E: Pair style COMB requires newton pair on

See the newton command.  This is a restriction to use the COMB
potential.

E: Pair style COMB requires atom attribute q

Self-explanatory.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open COMB potential file %s

The specified COMB potential file cannot be opened.  Check that the
path and name are correct.

E: Incorrect format in COMB potential file

Incorrect number of words per line in the potential file.

E: Illegal COMB parameter

One or more of the coefficients defined in the potential file is
invalid.

E: Potential file has duplicate entry

The potential file for a SW or Tersoff potential has more than
one entry for the same 3 ordered elements.

E: Potential file is missing an entry

The potential file for a SW or Tersoff potential does not have a
needed entry.

W: Pair COMB charge %.10f with force %.10f hit min barrier

Something is possibly wrong with your model.

W: Pair COMB charge %.10f with force %.10f hit max barrier

Something is possibly wrong with your model.

*/
