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

#include <boost/shared_ptr.hpp>

#include "ForceCompute.h"
#include "NeighborList.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define MAXLINE 1024

/*! \file EAMForceCompute.h
	\brief Declares the EAMForceCompute class
*/

#ifndef __EAMForceCompute_H__
#define __EAMForceCompute_H__

//! Computes Lennard-Jones forces on each particle
/*! The total pair force is summed for each particle when compute() is called. Forces are only summed between
	neighboring particles with a separation distance less than \c r_cut. A NeighborList must be provided
	to identify these neighbors. Calling compute() in this class will in turn result in a call to the 
	NeighborList's compute() to make sure that the neighbor list is up to date.
	
	Usage: Construct a EAMForceCompute, providing it an already constructed ParticleData and NeighborList.
	Then set parameters for all possible pairs of types by calling setParams. 
	
	Forces can be computed directly by calling compute() and then retrieved with a call to acquire(), but 
	a more typical usage will be to add the force compute to NVEUpdater or NVTUpdater. 
	
	\ingroup computes
*/
class EAMForceCompute : public ForceCompute
	{
	public:
		//! Constructs the compute
		EAMForceCompute(boost::shared_ptr<ParticleData> pdata, boost::shared_ptr<NeighborList> nlist, Scalar r_cut);
		
		//! Destructor
		virtual ~EAMForceCompute();
		
		//! Set the parameters for a single type pair
		virtual void setParams(unsigned int typ1, unsigned int typ2, Scalar lj1, Scalar lj2);
		
		//! Returns a list of log quantities this compute calculates
		virtual std::vector< std::string > getProvidedLogQuantities(); 
		
		//! Calculates the requested log value and returns it
		virtual Scalar getLogValue(const std::string& quantity, unsigned int timestep);
		
		//! Shifting modes that can be applied to the energy
		enum energyShiftMode
			{
			no_shift = 0,
			shift,
			xplor
			};
		
		//! Set the fraction of r_cut at which the xplor shifting is enabled (if that mode is set)
		void setXplorFraction(Scalar f) { m_xplor_fraction = f; }
		
		//! Set the mode to use for shifting the energy
		void setShiftMode(energyShiftMode mode) { m_shift_mode = mode; }

    //! read eam file and save in Funcfl
    int add_eam(const std::string& filename);

    //! convert Funcfl
    void file2array();

    //! generate final spline
    void array2spline();

		Scalar * __restrict__ m_rho;	//!< rho
		Scalar * __restrict__ m_fp;   //!< fp

    //! ....
    class Funcfl;
    std::vector<Funcfl> funcfl;
    int nfuncfl;
    
    //! ....
		int nrho,nr;
		int nfrho,nrhor,nz2r;
		Scalar dr,rdr,drho,rdrho;

    //! ....
    std::vector< std::vector<Scalar> > frho;
    std::vector< std::vector<Scalar> > rhor;
    std::vector< std::vector<Scalar> > z2r;

    //! eam spline
    std::vector< std::vector< std::vector<Scalar> > > frho_spline;
    std::vector< std::vector< std::vector<Scalar> > > rhor_spline;
    std::vector< std::vector< std::vector<Scalar> > > z2r_spline;

	protected:
		boost::shared_ptr<NeighborList> m_nlist;	//!< The neighborlist to use for the computation
		Scalar m_r_cut;								//!< Cuttoff radius beyond which the force is set to 0
		unsigned int m_ntypes;						//!< Store the width and height of lj1 and lj2 here
		energyShiftMode m_shift_mode;				//!< Store the mode with which to handle the energy shift at r_cut
		Scalar m_xplor_fraction;					//!< Fraction of r_cut at which to turn on the xplor shifting
		
		// This is a low level force summing class, it ONLY sums forces, and doesn't do high
		// level concepts like mixing. That is for the caller to handle. So, I only store 
		// lj1 and lj2 here
		Scalar * __restrict__ m_lj1;	//!< Parameter for computing forces (m_ntypes by m_ntypes array)
		Scalar * __restrict__ m_lj2;	//!< Parameter for computing forces	(m_ntypes by m_ntypes array)
		
		//! Actually compute the forces
		virtual void computeForces(unsigned int timestep);
  private:
	};

class EAMForceCompute::Funcfl {
public:
  Funcfl(){}
  Funcfl(const std::string& filename, int nrho, int nr, Scalar drho, Scalar dr, Scalar cut, Scalar mass,
         const std::vector<Scalar>& frho, const std::vector<Scalar>& rhor, const std::vector<Scalar>& zr) :
        filename(filename), nrho(nrho), nr(nr), drho(drho), dr(dr), cut(cut), mass(mass), frho(frho), rhor(rhor), zr(zr) {
  }
  std::string filename;
  int nrho,nr;
  Scalar drho,dr,cut,mass;
  std::vector<Scalar> frho;
  std::vector<Scalar> rhor;
  std::vector<Scalar> zr;
};

//! Exports the EAMForceCompute class to python
void export_EAMForceCompute();

#endif
