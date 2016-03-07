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

#include "ForceCompute.cuh"
#include "NeighborList.cuh"
#include "ParticleData.cuh"

/*! \file EAMForceGPU.cuh
	\brief Declares GPU kernel code for calculating the eam pair forces. Used by EAMForceComputeGPU.
*/

#ifndef __EAMFORCEGPU_CUH__
#define __EAMFORCEGPU_CUH__

//! options struct for passing additional options to gpu_compute_eam_forces
struct eam_options
	{
	float r_cutsq;			//!< cutoff distance squared
	int block_size;			//!< block size to execute on
	bool ulf_workaround;	//!< Set to true to enable the ULF workaroudn
  
    //! ....
		int nrho,nr;
		float rdr,rdrho;
	};

struct gpu_eam_spline {
	float4 *frho_spline;
	float4 *rhor_spline;
	float4 *z2r_spline;	
	float4 *frho_spline_;
	float4 *rhor_spline_;
	float4 *z2r_spline_;
};

//! Kernel driver that computes eam forces on the GPU for EAMForceComputeGPU
cudaError_t gpu_compute_eam_forces(const gpu_force_data_arrays& force_data, const gpu_pdata_arrays &pdata, const gpu_boxsize &box,
                                   const gpu_nlist_array &nlist, const eam_options& opt, const gpu_eam_spline &spline);

#endif
