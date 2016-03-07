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

#include "EAMForceCompute.h"
#include "NeighborList.h"
#include "EAMForceGPU.cuh"

#include <boost/shared_ptr.hpp>

/*! \file EAMForceComputeGPU.h
	\brief Declares the class EAMForceComputeGPU
*/

#ifndef __EAMForceComputeGPU_H__
#define __EAMForceComputeGPU_H__

//! Computes Lennard-Jones forces on each particle using the GPU
/*! Calculates the same forces as EAMForceCompute, but on the GPU.
	
	The GPU kernel for calculating the forces is in ljforcesum_kernel.cu.
	\ingroup computes
*/
class EAMForceComputeGPU : public EAMForceCompute
	{
	public:
		//! Constructs the compute
		EAMForceComputeGPU(boost::shared_ptr<ParticleData> pdata, boost::shared_ptr<NeighborList> nlist, Scalar r_cut);
		
		//! Destructor
		virtual ~EAMForceComputeGPU();
		
		//! Set the parameters for a single type pair
		virtual void setParams(unsigned int typ1, unsigned int typ2, Scalar lj1, Scalar lj2);
		
		//! Sets the block size to run at
		void setBlockSize(int block_size);
    
		eam_options opt;
    vector<gpu_eam_spline> spline;
	protected:
		vector<float2 *> d_coeffs;		//!< Pointer to the coefficients on the GPU
		float2 * h_coeffs;				//!< Pointer to the coefficients on the host
		int m_block_size;				//!< The block size to run on the GPU
		bool m_ulf_workaround;			//!< Stores decision made by the constructor whether to enable the ULF workaround

		//! Actually compute the forces
		virtual void computeForces(unsigned int timestep);
	};
	
//! Exports the EAMForceComputeGPU class to python
void export_EAMForceComputeGPU();

#endif
