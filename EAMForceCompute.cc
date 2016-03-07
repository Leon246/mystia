/*
Highly Optimized Object-Oriented Molecular Dynamics (HOOMD) Open
Source Software License
Copyright (c) 2008 Ames Laboratory Iowa State University
All rights reserved.

Redistribution and use of HOOMD, in source and binary forms, with or
without modification, are permitted, provided that the following
conditions are met:

* Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names HOOMD's
contributors may be used to endorse or promote products derived from this
software without specific prior written permission.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND
CONTRIBUTORS ``AS IS''  AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 

IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS  BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef WIN32
#pragma warning( push )
#pragma warning( disable : 4103 4244 )
#endif

#include <boost/python.hpp>
using namespace boost::python;

#include "EAMForceCompute.h"
#include <stdexcept>

/*! \file EAMForceCompute.cc
	\brief Defines the EAMForceCompute class
*/

using namespace std;

/*! \param pdata Particle Data to compute forces on
 	\param nlist Neighborlist to use for computing the forces
	\param r_cut Cuttoff radius beyond which the force is 0
	\post memory is allocated and all parameters lj1 and lj2 are set to 0.0
*/
EAMForceCompute::EAMForceCompute(boost::shared_ptr<ParticleData> pdata, boost::shared_ptr<NeighborList> nlist, Scalar r_cut) 
	: ForceCompute(pdata), m_nlist(nlist), m_r_cut(r_cut), m_shift_mode(no_shift), m_xplor_fraction(Scalar(2.0/3.0))
	{
	assert(m_pdata);
	assert(m_nlist);
	
	if (r_cut < 0.0)
		{
		cerr << endl << "***Error! Negative r_cut in EAMForceCompute makes no sense" << endl << endl;
		throw runtime_error("Error initializing EAMForceCompute");
		}
	
	// initialize the number of types value
	m_ntypes = m_pdata->getNTypes();
	assert(m_ntypes > 0);
	
	// allocate data for lj1 and lj2
	m_lj1 = new Scalar[m_ntypes*m_ntypes];
	m_lj2 = new Scalar[m_ntypes*m_ntypes];
	
	// sanity check
	assert(m_lj1 != NULL && m_lj2 != NULL);
	
	// initialize the parameters to 0;
	memset((void*)m_lj1, 0, sizeof(Scalar)*m_ntypes*m_ntypes);
	memset((void*)m_lj2, 0, sizeof(Scalar)*m_ntypes*m_ntypes);
  
	// allocate data on the host
	unsigned int num_particles = m_pdata->getN();
	m_rho = new Scalar[num_particles];
	m_fp  = new Scalar[num_particles];
	}
	

EAMForceCompute::~EAMForceCompute()
	{
	// deallocate our memory
	delete[] m_lj1;
	delete[] m_lj2;
	m_lj1 = NULL;
	m_lj2 = NULL;
  // deallocate rho, fp
  delete[] m_rho;
  delete[] m_fp;
	}
		
int EAMForceCompute::add_eam(const std::string &filename)
{
  char *ptr;
  char line[MAXLINE];
  int nrho,nr;
  double drho,dr,cut,mass;
  vector<Scalar> tmpfrho;
  vector<Scalar> tmprhor;
  vector<Scalar> tmpzr;

  ifstream ifs;
  ifs.open(filename.c_str(), ifstream::in);
  ifs.getline(line,MAXLINE);
  ifs.getline(line,MAXLINE);
  int tmp;
  sscanf(line,"%d %lg",&tmp,&mass);
  ifs.getline(line,MAXLINE);
  sscanf(line,"%d %lg %d %lg %lg",&nrho,&drho,&nr,&dr,&cut);

  for (int i = 0; i < nrho; i++) {
    ifs.getline(line,MAXLINE);
    ptr = strtok(line," \t\n\r\f");
    tmpfrho.push_back(Scalar(atof(ptr)));
    while (ptr = strtok(NULL," \t\n\r\f")) {
      tmpfrho.push_back(Scalar(atof(ptr)));
      i++;
    }
  }
  for (int i = 0; i < nr; i++) {
    ifs.getline(line,MAXLINE);
    ptr = strtok(line," \t\n\r\f");
    tmpzr.push_back(Scalar(atof(ptr)));
    while (ptr = strtok(NULL," \t\n\r\f")) {
      tmpzr.push_back(Scalar(atof(ptr)));
      i++;
    }
  }
  for (int i = 0; i < nr; i++) {
    ifs.getline(line,MAXLINE);
    ptr = strtok(line," \t\n\r\f");
    tmprhor.push_back(Scalar(atof(ptr)));
    while (ptr = strtok(NULL," \t\n\r\f")) {
      tmprhor.push_back(Scalar(atof(ptr)));
      i++;
    }
  }
  funcfl.push_back(Funcfl(filename, nrho, nr, Scalar(drho), Scalar(dr), Scalar(cut), Scalar(mass), tmpfrho, tmprhor, tmpzr));
  nfuncfl = funcfl.size();
  ifs.close();
  return funcfl.size()-1;
}

void EAMForceCompute::file2array()
{
  int i,j,k,m,n;
  Scalar sixth = 1.0/6.0;
  vector<Scalar> tmp;

  // determine max function params from all active funcfl files

  Scalar rmax,rhomax;
  dr = drho = rmax = rhomax = 0.0;

  for (int i = 0; i < nfuncfl; i++) {
    Funcfl *file = &funcfl[i];
    dr = MAX(dr,file->dr);
    drho = MAX(drho,file->drho);
    rmax = MAX(rmax,(file->nr-1) * file->dr);
    rhomax = MAX(rhomax,(file->nrho-1) * file->drho);
  }

  // set nr,nrho from cutoff and spacings
  // 0.5 is for round-off in divide

  nr = static_cast<int> (rmax/dr + 0.5);
  nrho = static_cast<int> (rhomax/drho + 0.5);

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------
  
  nfrho = nfuncfl;

  // interpolate each file's frho to a single grid and cutoff

  Scalar r,p,cof1,cof2,cof3,cof4;
  
  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *file = &funcfl[i];
    for (m = 0; m < nrho; m++) {
      r = Scalar(m)*drho;
      p = r/file->drho + 1.0;
      k = static_cast<int> (p);
      k = MIN(k,file->nrho-2);
      k = MAX(k,2);
      p -= k;
      p = MIN(p,2.0);
      cof1 = -sixth*p*(p-1.0)*(p-2.0);
      cof2 = 0.5*(p*p-1.0)*(p-2.0);
      cof3 = -0.5*p*(p+1.0)*(p-2.0);
      cof4 = sixth*p*(p*p-1.0);
      tmp.push_back(cof1*file->frho[k-2] + cof2*file->frho[k-1] + cof3*file->frho[k] + cof4*file->frho[k+1]);
    }
    frho.push_back(tmp);
    tmp.clear();
    n++;
  }

  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------

  // allocate rhor arrays
  // nrhor = # of funcfl files

  nrhor = nfuncfl;

  // interpolate each file's rhor to a single grid and cutoff

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *file = &funcfl[i];
    for (m = 0; m < nr; m++) {
      r = Scalar(m)*dr;
      p = r/file->dr + 1.0;
      k = static_cast<int> (p);
      k = MIN(k,file->nr-2);
      k = MAX(k,2);
      p -= k;
      p = MIN(p,2.0);
      cof1 = -sixth*p*(p-1.0)*(p-2.0);
      cof2 = 0.5*(p*p-1.0)*(p-2.0);
      cof3 = -0.5*p*(p+1.0)*(p-2.0);
      cof4 = sixth*p*(p*p-1.0);
      tmp.push_back(cof1*file->rhor[k-2] + cof2*file->rhor[k-1] + cof3*file->rhor[k] + cof4*file->rhor[k+1]);
    }
    rhor.push_back(tmp);
    tmp.clear();
    n++;
  }

  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------

  // allocate z2r arrays
  // nz2r = N*(N+1)/2 where N = # of funcfl files

  nz2r = nfuncfl*(nfuncfl+1)/2;

  // create a z2r array for each file against other files, only for I >= J
  // interpolate zri and zrj to a single grid and cutoff

  Scalar zri,zrj;

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *ifile = &funcfl[i];
    for (j = 0; j <= i; j++) {
      Funcfl *jfile = &funcfl[j];
      for (m = 0; m < nr; m++) {
	      r = Scalar(m)*dr;

	      p = r/ifile->dr + 1.0;
	      k = static_cast<int> (p);
	      k = MIN(k,ifile->nr-2);
	      k = MAX(k,2);
	      p -= k;
	      p = MIN(p,2.0);
	      cof1 = -sixth*p*(p-1.0)*(p-2.0);
	      cof2 = 0.5*(p*p-1.0)*(p-2.0);
	      cof3 = -0.5*p*(p+1.0)*(p-2.0);
	      cof4 = sixth*p*(p*p-1.0);
	      zri = cof1*ifile->zr[k-2] + cof2*ifile->zr[k-1] + cof3*ifile->zr[k] + cof4*ifile->zr[k+1];

	      p = r/jfile->dr + 1.0;
	      k = static_cast<int> (p);
	      k = MIN(k,jfile->nr-2);
	      k = MAX(k,2);
	      p -= k;
	      p = MIN(p,2.0);
	      cof1 = -sixth*p*(p-1.0)*(p-2.0);
	      cof2 = 0.5*(p*p-1.0)*(p-2.0);
	      cof3 = -0.5*p*(p+1.0)*(p-2.0);
	      cof4 = sixth*p*(p*p-1.0);
	      zrj = cof1*jfile->zr[k-2] + cof2*jfile->zr[k-1] + cof3*jfile->zr[k] + cof4*jfile->zr[k+1];

        tmp.push_back(27.2*0.529 * zri*zrj);
      }
      z2r.push_back(tmp);
      tmp.clear();
      n++;
    }
  }
}

void EAMForceCompute::array2spline()
{
  rdr = 1.0/dr;
  rdrho = 1.0/drho;
  
  std::vector<Scalar> tmp1d;
  std::vector< std::vector<Scalar> > tmp2d;

  frho_spline.clear();
  rhor_spline.clear();
  z2r_spline.clear();

  tmp1d.clear();
  tmp2d.clear();
  frho_spline.clear();
  for (int i = 0; i < nfrho; i++) {
    // construct 2d vector
    for (int j = 0; j < nrho+1; j++) {
      for (int k = 0; k < 7; k++)
        tmp1d.push_back(Scalar(0.0));
      tmp2d.push_back(tmp1d);
      tmp1d.clear();
    }

    int n = nrho;
    Scalar delta = drho;

    // interpolate()
    for (int m = 1; m <= n; m++) tmp2d[m][6] = frho[i][m-1];

    tmp2d[1][5] = tmp2d[2][6] - tmp2d[1][6];
    tmp2d[2][5] = 0.5 * (tmp2d[3][6]-tmp2d[1][6]);
    tmp2d[n-1][5] = 0.5 * (tmp2d[n][6]-tmp2d[n-2][6]);
    tmp2d[n][5] = tmp2d[n][6] - tmp2d[n-1][6];
    
    for (int m = 3; m <= n-2; m++)
      tmp2d[m][5] = ((tmp2d[m-2][6]-tmp2d[m+2][6]) + 8.0*(tmp2d[m+1][6]-tmp2d[m-1][6])) / 12.0;
    
    for (int m = 1; m <= n-1; m++) {
      tmp2d[m][4] = 3.0*(tmp2d[m+1][6]-tmp2d[m][6]) - 2.0*tmp2d[m][5] - tmp2d[m+1][5];
      tmp2d[m][3] = tmp2d[m][5] + tmp2d[m+1][5] - 2.0*(tmp2d[m+1][6]-tmp2d[m][6]);
    }
    
    tmp2d[n][4] = 0.0;
    tmp2d[n][3] = 0.0;
    
    for (int m = 1; m <= n; m++) {
      tmp2d[m][2] = tmp2d[m][5]/delta;
      tmp2d[m][1] = 2.0*tmp2d[m][4]/delta;
      tmp2d[m][0] = 3.0*tmp2d[m][3]/delta;
    }

    // 3d vector
    frho_spline.push_back(tmp2d);
  }

  tmp1d.clear();
  tmp2d.clear();
  rhor_spline.clear();
  for (int i = 0; i < nrhor; i++) {
    // construct 2d vector
    for (int j = 0; j < nr+1; j++) {
      for (int k = 0; k < 7; k++)
        tmp1d.push_back(Scalar(0.0));
      tmp2d.push_back(tmp1d);
      tmp1d.clear();
    }

    int n = nr;
    Scalar delta = dr;

    // interpolate()
    for (int m = 1; m <= n; m++) tmp2d[m][6] = rhor[i][m-1];

    tmp2d[1][5] = tmp2d[2][6] - tmp2d[1][6];
    tmp2d[2][5] = 0.5 * (tmp2d[3][6]-tmp2d[1][6]);
    tmp2d[n-1][5] = 0.5 * (tmp2d[n][6]-tmp2d[n-2][6]);
    tmp2d[n][5] = tmp2d[n][6] - tmp2d[n-1][6];
    
    for (int m = 3; m <= n-2; m++)
      tmp2d[m][5] = ((tmp2d[m-2][6]-tmp2d[m+2][6]) + 8.0*(tmp2d[m+1][6]-tmp2d[m-1][6])) / 12.0;
    
    for (int m = 1; m <= n-1; m++) {
      tmp2d[m][4] = 3.0*(tmp2d[m+1][6]-tmp2d[m][6]) - 2.0*tmp2d[m][5] - tmp2d[m+1][5];
      tmp2d[m][3] = tmp2d[m][5] + tmp2d[m+1][5] - 2.0*(tmp2d[m+1][6]-tmp2d[m][6]);
    }
    
    tmp2d[n][4] = 0.0;
    tmp2d[n][3] = 0.0;
    
    for (int m = 1; m <= n; m++) {
      tmp2d[m][2] = tmp2d[m][5]/delta;
      tmp2d[m][1] = 2.0*tmp2d[m][4]/delta;
      tmp2d[m][0] = 3.0*tmp2d[m][3]/delta;
    }

    // 3d vector
    rhor_spline.push_back(tmp2d);
  }
  
  tmp1d.clear();
  tmp2d.clear();
  z2r_spline.clear();
  for (int i = 0; i < nz2r; i++) {
    // construct 2d vector
    for (int j = 0; j < nr+1; j++) {
      for (int k = 0; k < 7; k++)
        tmp1d.push_back(Scalar(0.0));
      tmp2d.push_back(tmp1d);
      tmp1d.clear();
    }

    int n = nr;
    Scalar delta = dr;

    // interpolate()
    for (int m = 1; m <= n; m++) tmp2d[m][6] = z2r[i][m-1];

    tmp2d[1][5] = tmp2d[2][6] - tmp2d[1][6];
    tmp2d[2][5] = 0.5 * (tmp2d[3][6]-tmp2d[1][6]);
    tmp2d[n-1][5] = 0.5 * (tmp2d[n][6]-tmp2d[n-2][6]);
    tmp2d[n][5] = tmp2d[n][6] - tmp2d[n-1][6];
    
    for (int m = 3; m <= n-2; m++)
      tmp2d[m][5] = ((tmp2d[m-2][6]-tmp2d[m+2][6]) + 8.0*(tmp2d[m+1][6]-tmp2d[m-1][6])) / 12.0;
    
    for (int m = 1; m <= n-1; m++) {
      tmp2d[m][4] = 3.0*(tmp2d[m+1][6]-tmp2d[m][6]) - 2.0*tmp2d[m][5] - tmp2d[m+1][5];
      tmp2d[m][3] = tmp2d[m][5] + tmp2d[m+1][5] - 2.0*(tmp2d[m+1][6]-tmp2d[m][6]);
    }
    
    tmp2d[n][4] = 0.0;
    tmp2d[n][3] = 0.0;
    
    for (int m = 1; m <= n; m++) {
      tmp2d[m][2] = tmp2d[m][5]/delta;
      tmp2d[m][1] = 2.0*tmp2d[m][4]/delta;
      tmp2d[m][0] = 3.0*tmp2d[m][3]/delta;
    }

    // 3d vector
    z2r_spline.push_back(tmp2d);
  }
}

/*! \post The parameters \a lj1 and \a lj2 are set for the pairs \a typ1, \a typ2 and \a typ2, \a typ1.
	\note \a lj? are low level parameters used in the calculation. In order to specify
	these for a normal lennard jones formula (with alpha), they should be set to the following.
	- \a lj1 = 4.0 * epsilon * pow(sigma,12.0)
	- \a lj2 = alpha * 4.0 * epsilon * pow(sigma,6.0);
	
	Setting the parameters for typ1,typ2 automatically sets the same parameters for typ2,typ1: there
	is no need to call this funciton for symmetric pairs. Any pairs that this function is not called
	for will have lj1 and lj2 set to 0.0.
	
	\param typ1 Specifies one type of the pair
	\param typ2 Specifies the second type of the pair
	\param lj1 First parameter used to calcluate forces
	\param lj2 Second parameter used to calculate forces
*/
void EAMForceCompute::setParams(unsigned int typ1, unsigned int typ2, Scalar lj1, Scalar lj2)
	{
	if (typ1 >= m_ntypes || typ2 >= m_ntypes)
		{
		cerr << endl << "***Error! Trying to set EAM params for a non existant type! " << typ1 << "," << typ2 << endl << endl;
		throw runtime_error("Error setting parameters in EAMForceCompute");
		}
	
	// set lj1 in both symmetric positions in the matrix	
	m_lj1[typ1*m_ntypes + typ2] = lj1;
	m_lj1[typ2*m_ntypes + typ1] = lj1;
	
	// set lj2 in both symmetric positions in the matrix
	m_lj2[typ1*m_ntypes + typ2] = lj2;
	m_lj2[typ2*m_ntypes + typ1] = lj2;

  add_eam("Au_u3.eam");
  file2array();
  array2spline();

  Funcfl *file_ = &funcfl[0];
  m_r_cut = file_->cut;
/*
  //! verifying spline data
  ofstream ofs;
  ofs.open("spline.txt", ofstream::out);
  Funcfl *file = &funcfl[0];
  for (int j = 0; j <= nrho; j++) {
    ofs << "frho[" << j << "]: " << file->frho[j] << endl;
  }
  for (int j = 0; j <= nr; j++) {
    ofs << "rhor[" << j << "]: " << file->rhor[j] << endl;
  }
  for (int j = 0; j <= nr; j++) {
    ofs << "zr[" << j << "]:   " << file->zr[j] << endl;
  }
  ofs <<"------------------------------------------"<< endl;
  for (int j = 0; j <= nrho; j++) {
    ofs << "frho[" << j << "]: " << frho[0][j] << endl;
  }
  for (int j = 0; j <= nrho; j++) {
    ofs << "rhor[" << j << "]: " << rhor[0][j] << endl;
  }
  for (int j = 0; j <= nrho; j++) {
    ofs << "z2r[" << j << "]:  " << z2r[0][j] << endl;
  }
  ofs <<"------------------------------------------"<< endl;
  for (int i = 0; i <= nrho; i++) {
    ofs << "frho_spline[" << i << "]: " << frho_spline[0][i][0] << " ";
    ofs << frho_spline[0][i][1] << " ";
    ofs << frho_spline[0][i][2] << " ";
    ofs << frho_spline[0][i][3] << " ";
    ofs << frho_spline[0][i][4] << " ";
    ofs << frho_spline[0][i][5] << " ";
    ofs << frho_spline[0][i][6] << endl;
  }  
  for (int i = 0; i <= nr; i++) {
    ofs << "rhor_spline[" << i << "]: " << rhor_spline[0][i][0] << " ";
    ofs << rhor_spline[0][i][1] << " ";
    ofs << rhor_spline[0][i][2] << " ";
    ofs << rhor_spline[0][i][3] << " ";
    ofs << rhor_spline[0][i][4] << " ";
    ofs << rhor_spline[0][i][5] << " ";
    ofs << rhor_spline[0][i][6] << endl;
  }
  for (int i = 0; i <= nr; i++) {
    ofs << "z2r_spline[" << i << "]:  " << z2r_spline[0][i][0] << " ";
    ofs << z2r_spline[0][i][1] << " ";
    ofs << z2r_spline[0][i][2] << " ";
    ofs << z2r_spline[0][i][3] << " ";
    ofs << z2r_spline[0][i][4] << " ";
    ofs << z2r_spline[0][i][5] << " ";
    ofs << z2r_spline[0][i][6] << endl;
  }
  ofs.close();
  */
	}
	
/*! EAMForceCompute provides
	- \c pair_eam_energy
*/
std::vector< std::string > EAMForceCompute::getProvidedLogQuantities()
	{
	vector<string> list;
	list.push_back("pair_eam_energy");
	return list;
	}
	
Scalar EAMForceCompute::getLogValue(const std::string& quantity, unsigned int timestep)
	{
	if (quantity == string("pair_eam_energy"))
		{
		compute(timestep);
		return calcEnergySum();
		}
	else
		{
		cerr << endl << "***Error! " << quantity << " is not a valid log quantity for EAMForceCompute" << endl << endl;
		throw runtime_error("Error getting log value");
		}
	}

/*! \post The lennard jones forces are computed for the given timestep. The neighborlist's
 	compute method is called to ensure that it is up to date.
	
	\param timestep specifies the current time step of the simulation
*/
void EAMForceCompute::computeForces(unsigned int timestep)
	{
  int m;
  Scalar r,p,rhoip,rhojp,z2,z2p,recip,phip,psip,phi,fpair;
  vector<Scalar> coeff;
	// start by updating the neighborlist
	m_nlist->compute(timestep);
	
	// start the profile for this compute
	if (m_prof) m_prof->push("EAM pair");
	
	// depending on the neighborlist settings, we can take advantage of newton's third law
	// to reduce computations at the cost of memory access complexity: set that flag now
	bool third_law = m_nlist->getStorageMode() == NeighborList::half;
	
	// access the neighbor list
	const vector< vector< unsigned int > >& full_list = m_nlist->getList();

	// access the particle data
	const ParticleDataArraysConst& arrays = m_pdata->acquireReadOnly(); 
	// sanity check
	assert(arrays.x != NULL && arrays.y != NULL && arrays.z != NULL);
	
	// get a local copy of the simulation box too
	const BoxDim& box = m_pdata->getBox();
	// sanity check
	assert(box.xhi > box.xlo && box.yhi > box.ylo && box.zhi > box.zlo);	
	
	// create a temporary copy of r_cut sqaured
	Scalar r_cut_sq = m_r_cut * m_r_cut;
	
	// precalculate box lenghts for use in the periodic imaging
	Scalar Lx = box.xhi - box.xlo;
	Scalar Ly = box.yhi - box.ylo;
	Scalar Lz = box.zhi - box.zlo;
	
	// tally up the number of forces calculated
	int64_t n_calc = 0;
	
	// need to start from a zero force, energy and virial
	// (MEM TRANSFER: 5*N scalars)
	memset(m_fx, 0, sizeof(Scalar)*arrays.nparticles);
	memset(m_fy, 0, sizeof(Scalar)*arrays.nparticles);
	memset(m_fz, 0, sizeof(Scalar)*arrays.nparticles);
	memset(m_pe, 0, sizeof(Scalar)*arrays.nparticles);
	memset(m_virial, 0, sizeof(Scalar)*arrays.nparticles);

  // start from zero rho, fp
	memset(m_rho, 0, sizeof(Scalar)*arrays.nparticles);
	memset(m_fp, 0, sizeof(Scalar)*arrays.nparticles);

	// for each particle
	for (unsigned int i = 0; i < arrays.nparticles; i++) {
		// access the particle's position and type (MEM TRANSFER: 4 scalars)
		Scalar xi = arrays.x[i];
		Scalar yi = arrays.y[i];
		Scalar zi = arrays.z[i];
		unsigned int typei = arrays.type[i];
		// sanity check
		assert(typei < m_pdata->getNTypes());
		
		Scalar rho = 0.0;
		// loop over all of the neighbors of this particle
		const vector< unsigned int >& list = full_list[i];
		const unsigned int size = (unsigned int)list.size();
		for (unsigned int j = 0; j < size; j++)	{
			// increment our calculation counter
			n_calc++;
			
			// access the index of this neighbor (MEM TRANSFER: 1 scalar)
			unsigned int k = list[j];
			// sanity check
			assert(k < m_pdata->getN());
				
			// calculate dr (MEM TRANSFER: 3 scalars / FLOPS: 3)
			Scalar dx = xi - arrays.x[k];
			Scalar dy = yi - arrays.y[k];
			Scalar dz = zi - arrays.z[k];
			
			// access the type of the neighbor particle (MEM TRANSFER: 1 scalar
			unsigned int typej = arrays.type[k];
			// sanity check
			assert(typej < m_pdata->getNTypes());
			
			// apply periodic boundary conditions (FLOPS: 9 (worst case: first branch is missed, the 2nd is taken and the add is done)
			if (dx >= box.xhi)
				dx -= Lx;
			else
			if (dx < box.xlo)
				dx += Lx;

			if (dy >= box.yhi)
				dy -= Ly;
			else
			if (dy < box.ylo)
				dy += Ly;
	
			if (dz >= box.zhi)
				dz -= Lz;
			else
			if (dz < box.zlo)
				dz += Lz;
			
			// start computing the force
			// calculate r squared (FLOPS: 5)
			Scalar rsq = dx*dx + dy*dy + dz*dz;
		
			// only compute the force if the particles are closer than the cuttoff (FLOPS: 1)
			if (rsq < r_cut_sq)	{
        p = sqrt(rsq)*rdr + 1.0;
	      m = static_cast<int> (p);
	      m = MIN(m,nr-1);
	      p -= m;
	      p = MIN(p,1.0);
	      coeff = rhor_spline[0][m];
	      m_rho[i] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
	      if (third_law) {
	        coeff = rhor_spline[0][m];
	        m_rho[k] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
	      }
      }
		}
	}

  // embedding energy to each atom
  for (unsigned int i = 0; i < arrays.nparticles; i++) {
    p = m_rho[i]*rdrho + 1.0;
    m = static_cast<int> (p);
    m = MAX(1,MIN(m,nrho-1));
    p -= m;
    p = MIN(p,1.0);
    coeff = frho_spline[0][m];
    m_fp[i] = (coeff[0]*p + coeff[1])*p + coeff[2];
    m_pe[i] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
  }

  // for each particle
	for (unsigned int i = 0; i < arrays.nparticles; i++)
		{
		// access the particle's position and type (MEM TRANSFER: 4 scalars)
		Scalar xi = arrays.x[i];
		Scalar yi = arrays.y[i];
		Scalar zi = arrays.z[i];
		unsigned int typei = arrays.type[i];
		// sanity check
		assert(typei < m_pdata->getNTypes());

		// initialize current particle force, potential energy, and virial to 0
		Scalar fxi = 0.0;
		Scalar fyi = 0.0;
		Scalar fzi = 0.0;
		Scalar pei = 0.0;
		Scalar viriali = 0.0;
		
		// loop over all of the neighbors of this particle
		const vector< unsigned int >& list = full_list[i];
		const unsigned int size = (unsigned int)list.size();
		for (unsigned int j = 0; j < size; j++)
			{
			// increment our calculation counter
			n_calc++;
			
			// access the index of this neighbor (MEM TRANSFER: 1 scalar)
			unsigned int k = list[j];
			// sanity check
			assert(k < m_pdata->getN());
				
			// calculate dr (MEM TRANSFER: 3 scalars / FLOPS: 3)
			Scalar dx = xi - arrays.x[k];
			Scalar dy = yi - arrays.y[k];
			Scalar dz = zi - arrays.z[k];
			
			// access the type of the neighbor particle (MEM TRANSFER: 1 scalar
			unsigned int typej = arrays.type[k];
			// sanity check
			assert(typej < m_pdata->getNTypes());
			
			// apply periodic boundary conditions (FLOPS: 9 (worst case: first branch is missed, the 2nd is taken and the add is done)
			if (dx >= box.xhi)
				dx -= Lx;
			else
			if (dx < box.xlo)
				dx += Lx;

			if (dy >= box.yhi)
				dy -= Ly;
			else
			if (dy < box.ylo)
				dy += Ly;
	
			if (dz >= box.zhi)
				dz -= Lz;
			else
			if (dz < box.zlo)
				dz += Lz;
			
			// start computing the force
			// calculate r squared (FLOPS: 5)
			Scalar rsq = dx*dx + dy*dy + dz*dz;
		
			// only compute the force if the particles are closer than the cuttoff (FLOPS: 1)
			if (rsq < r_cut_sq) {
				r = sqrt(rsq);
	      p = r*rdr + 1.0;
	      m = static_cast<int> (p);
	      m = MIN(m,nr-1);
	      p -= m;
	      p = MIN(p,1.0);

	      // rhoip = derivative of (density at atom j due to atom i)
	      // rhojp = derivative of (density at atom i due to atom j)
	      // phi = pair potential energy
	      // phip = phi'
	      // z2 = phi * r
	      // z2p = (phi * r)' = (phi' r) + phi
	      // psip needs both fp[i] and fp[j] terms since r_ij appears in two
	      //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
	      //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip

	      coeff = rhor_spline[0][m];
	      rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
	      coeff = rhor_spline[0][m];
	      rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];

	      coeff = z2r_spline[0][m];
	      z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
	      z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

	      recip = 1.0/r;
	      phi = z2*recip;
	      phip = z2p*recip - phi*recip;
	      psip = m_fp[i]*rhojp + m_fp[k]*rhoip + phip;
	      fpair = -psip*recip;

	      // compute the virial (FLOPS: 2)
				// note the sign in the virial calculation, this is because dx,dy,dz are \vec{r}_{ji} thus
				// there is no - in the 1/6 to compensate	
				Scalar pair_virial = Scalar(1.0/6.0) * rsq * fpair;
				
				// add the force, potential energy and virial to the particle i
				// (FLOPS: 8)
				fxi += dx*fpair;
				fyi += dy*fpair;
				fzi += dz*fpair;
				pei += phi*0.5;
				viriali += pair_virial;
				
				// add the force to particle j if we are using the third law (MEM TRANSFER: 10 scalars / FLOPS: 8)
				if (third_law)
					{
					m_fx[k] -= dx*fpair;
					m_fy[k] -= dy*fpair;
					m_fz[k] -= dz*fpair;
					m_pe[k] += phi*0.5;
					m_virial[k] += pair_virial;
					}
				}			
			}
			
		// finally, increment the force, potential energy and virial for particle i
		// (MEM TRANSFER: 10 scalars / FLOPS: 5)
		m_fx[i] += fxi;
		m_fy[i] += fyi;
		m_fz[i] += fzi;
		m_pe[i] += pei;
		m_virial[i] += viriali;
		}

	m_pdata->release();
	
	#ifdef ENABLE_CUDA
	// the force data is now only up to date on the cpu
	m_data_location = cpu;
	#endif

	int64_t flops = m_pdata->getN() * 5 + n_calc * (3+5+9+1+9+6+8);
	if (m_shift_mode == shift)
		flops += n_calc * 5;
	else
	if (m_shift_mode == xplor)
		flops += n_calc * 16;

	if (third_law) flops += n_calc * 8;
	int64_t mem_transfer = m_pdata->getN() * (5+4+10)*sizeof(Scalar) + n_calc * (1+3+1)*sizeof(Scalar);
	if (third_law) mem_transfer += n_calc*10*sizeof(Scalar);
	if (m_prof) m_prof->pop(flops, mem_transfer);
	}

void export_EAMForceCompute()
	{
	scope in_lj = class_<EAMForceCompute, boost::shared_ptr<EAMForceCompute>, bases<ForceCompute>, boost::noncopyable >
		("EAMForceCompute", init< boost::shared_ptr<ParticleData>, boost::shared_ptr<NeighborList>, Scalar >())
		.def("setParams", &EAMForceCompute::setParams)
		.def("setXplorFraction", &EAMForceCompute::setXplorFraction)
		.def("setShiftMode", &EAMForceCompute::setShiftMode)
		;
		
	enum_<EAMForceCompute::energyShiftMode>("energyShiftMode")
		.value("no_shift", EAMForceCompute::no_shift)
		.value("shift", EAMForceCompute::shift)
		.value("xplor", EAMForceCompute::xplor)
		;
	}

#ifdef WIN32
#pragma warning( pop )
#endif
