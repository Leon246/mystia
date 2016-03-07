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

#include "gpu_settings.h"
#include "EAMForceGPU.cuh"

#ifdef WIN32
#include <cassert>
#else
#include <assert.h>
#endif

//! tabulated function coeffs is set
int IsSet = 0;

//! Texture for reading particle positions
texture<float4, 1, cudaReadModeElementType> pdata_pos_tex;
// texture<float4, 1, cudaReadModeElementType> force_data_tex;

//! Texture for tabulated function coeffs
// coeff.x = spline[3]
// coeff.y = spline[4]
// coeff.z = spline[5]
// coeff.w = spline[6]
texture<float4, 1, cudaReadModeElementType> frho_coeff;
texture<float4, 1, cudaReadModeElementType> rhor_coeff;
texture<float4, 1, cudaReadModeElementType> z2r_coeff;
//! Texture for tabulated function coeffs
// coeff.x = spline[0]
// coeff.y = spline[1]
// coeff.z = spline[2]
// texture<float4, 1, cudaReadModeElementType> frho_coeff_;
// texture<float4, 1, cudaReadModeElementType> rhor_coeff_;
// texture<float4, 1, cudaReadModeElementType> z2r_coeff_;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

//! Kernel for calculating lj forces
/*! This kernel is called to calculate the lennard-jones forces on all N particles

	\param force_data Device memory array to write calculated forces to
	\param pdata Particle data on the GPU to calculate forces on
	\param nlist Neigbhor list data on the GPU to use to calculate the forces
	\param d_coeffs Coefficients to the lennard jones force.
	\param coeff_width Width of the coefficient matrix
	\param r_cutsq Precalculated r_cut*r_cut, where r_cut is the radius beyond which forces are
		set to 0
	\param rcut6inv Precalculated 1/r_cut**6
	\param xplor_denom_inv Precalculated 1/xplor denominator
	\param r_on_sq Precalculated r_on*r_on (for xplor)
	\param box Box dimensions used to implement periodic boundary conditions
	
	\a coeffs is a pointer to a matrix in memory. \c coeffs[i*coeff_width+j].x is \a lj1 for the type pair \a i, \a j.
	Similarly, .y is the \a lj2 parameter. The values in d_coeffs are read into shared memory, so 
	\c coeff_width*coeff_width*sizeof(float2) bytes of extern shared memory must be allocated for the kernel call.
	
	Developer information:
	Each block will calculate the forces on a block of particles.
	Each thread will calculate the total force on one particle.
	The neighborlist is arranged in columns so that reads are fully coalesced when doing this.
*/
template<bool ulf_workaround> __global__ void gpu_compute_eam_embed_kernel(gpu_force_data_arrays force_data, gpu_pdata_arrays pdata, gpu_boxsize box, gpu_nlist_array nlist, 
                                                                           float r_cutsq, int nr, int nrho, float rdr, float rdrho)
{	
	// start by identifying which particle we are to handle
	unsigned int idx_local = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (idx_local >= pdata.local_num)
		return;
	
	unsigned int idx_global = idx_local + pdata.local_beg;
	
	// load in the length of the list (MEM_TRANSFER: 4 bytes)
	unsigned int n_neigh = nlist.n_neigh[idx_global];

	// read in the position of our particle. Texture reads of float4's are faster than global reads on compute 1.0 hardware
	// (MEM TRANSFER: 16 bytes)
	float4 pos = tex1Dfetch(pdata_pos_tex, idx_global);
	
	// initialize the force to 0
	float4 force = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
  float rho = 0.0f;
	float p;
	unsigned int m;
	float4 coeff = make_float4(0.0f, 0.0f, 0.0f, 0.0f);

	// prefetch neighbor index
	unsigned int cur_neigh = 0;
	unsigned int next_neigh = nlist.list[idx_global];

	// loop over neighbors
	// on pre C1060 hardware, there is a bug that causes rare and random ULFs when simply looping over n_neigh
	// the workaround (activated via the template paramter) is to loop over nlist.height and put an if (i < n_neigh)
	// inside the loop
	int n_loop;
	if (ulf_workaround)
		n_loop = nlist.height;
	else
		n_loop = n_neigh;
	
	for (int neigh_idx = 0; neigh_idx < n_loop; neigh_idx++) {
		if (!ulf_workaround || neigh_idx < n_neigh) {

		  // read the current neighbor index (MEM TRANSFER: 4 bytes)
		  // prefetch the next value and set the current one
		  cur_neigh = next_neigh;
		  next_neigh = nlist.list[nlist.pitch*(neigh_idx+1) + idx_global];
  		
		  // get the neighbor's position (MEM TRANSFER: 16 bytes)
		  float4 neigh_pos = tex1Dfetch(pdata_pos_tex, cur_neigh);
  		
		  // calculate dr (with periodic boundary conditions) (FLOPS: 3)
		  float dx = pos.x - neigh_pos.x;
		  float dy = pos.y - neigh_pos.y;
		  float dz = pos.z - neigh_pos.z;
  			
		  // apply periodic boundary conditions: (FLOPS 12)
		  dx -= box.Lx * rintf(dx * box.Lxinv);
		  dy -= box.Ly * rintf(dy * box.Lyinv);
		  dz -= box.Lz * rintf(dz * box.Lzinv);
  			
		  // calculate r squard (FLOPS: 5)
		  float rsq = dx*dx + dy*dy + dz*dz;
  		
      if (rsq >= r_cutsq)
			  rsq = 0.0f;
	
		  p = sqrt(rsq)*rdr + 1.0f;
		  m = floor(p);
		  m = MIN(m,nr-1);
		  p -= m;
		  p = MIN(p,1.0f);
		  coeff = tex1Dfetch(rhor_coeff, m);
		  rho += ((coeff.x*p+coeff.y)*p+coeff.z)*p+coeff.w;
		}
	}

  // embedding energy to each atom
  p = rho*rdrho + 1.0f;
	m = floor(p);
	m = MAX(1,MIN(m,nrho-1));
	p -= m;
	p = MIN(p, 1.0f);
	coeff = tex1Dfetch(frho_coeff, m);
	float fp = ((3.0f*coeff.x*p+2.0f*coeff.y)*p+coeff.z)*rdrho;
	float phi = ((coeff.x*p+coeff.y)*p+coeff.z)*p+coeff.w;
  force.x = 0.0f;
  force.y = 0.0f;
  force.z = phi;
  force.w = fp;
  
	force_data.force[idx_local] = force;
}

template<bool ulf_workaround> __global__ void gpu_compute_eam_force_kernel(gpu_force_data_arrays force_data, gpu_pdata_arrays pdata, gpu_boxsize box, gpu_nlist_array nlist,  float r_cutsq, unsigned int nr, unsigned int nrho, float rdr, float rdrho)
{
	// start by identifying which particle we are to handle
	unsigned int idx_local = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (idx_local >= pdata.local_num)
		return;
	
	unsigned int idx_global = idx_local + pdata.local_beg;
	
	// load in the length of the list (MEM_TRANSFER: 4 bytes)
	unsigned int n_neigh = nlist.n_neigh[idx_global];

	// read in the position of our particle. Texture reads of float4's are faster than global reads on compute 1.0 hardware
	// (MEM TRANSFER: 16 bytes)
	float4 pos = tex1Dfetch(pdata_pos_tex, idx_global);
  // get fake force
	// float4 fforce = tex1Dfetch(force_data_tex, idx_global);
  float4 fforce = force_data.force[idx_global];
	
	// initialize the force to 0
	float4 coeff;
	float4 force = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
	float virial = 0.0f;

	// prefetch neighbor index
	unsigned int cur_neigh = 0;
	unsigned int next_neigh = nlist.list[idx_global];

	// loop over neighbors
	// on pre C1060 hardware, there is a bug that causes rare and random ULFs when simply looping over n_neigh
	// the workaround (activated via the template paramter) is to loop over nlist.height and put an if (i < n_neigh)
	// inside the loop
	int n_loop;
	if (ulf_workaround)
		n_loop = nlist.height;
	else
		n_loop = n_neigh;
		
	for (int neigh_idx = 0; neigh_idx < n_loop; neigh_idx++)
		{
		if (!ulf_workaround || neigh_idx < n_neigh)
		{
		// read the current neighbor index (MEM TRANSFER: 4 bytes)
		// prefetch the next value and set the current one
		cur_neigh = next_neigh;
		next_neigh = nlist.list[nlist.pitch*(neigh_idx+1) + idx_global];
		
		// get the neighbor's position (MEM TRANSFER: 16 bytes)
		float4 neigh_pos = tex1Dfetch(pdata_pos_tex, cur_neigh);
    // get the neighbor's fake fp
		// float4 neigh_fforce = tex1Dfetch(force_data_tex, cur_neigh);
    float4 neigh_fforce = force_data.force[cur_neigh];
				
		// calculate dr (with periodic boundary conditions) (FLOPS: 3)
		float dx = pos.x - neigh_pos.x;
		float dy = pos.y - neigh_pos.y;
		float dz = pos.z - neigh_pos.z;
			
		// apply periodic boundary conditions: (FLOPS 12)
		dx -= box.Lx * rintf(dx * box.Lxinv);
		dy -= box.Ly * rintf(dy * box.Lyinv);
		dz -= box.Lz * rintf(dz * box.Lzinv);
			
		// calculate r squard (FLOPS: 5)
		float rsq = dx*dx + dy*dy + dz*dz;
		
		if (rsq > r_cutsq)
			rsq = 0.0f;
	
		float p = sqrt(rsq)*rdr + 1.0f;
		unsigned int m = floor(p);
		m = MIN(m,nr-1);
		p -= m;
		p = MIN(p,1.0f);
		coeff = tex1Dfetch(rhor_coeff, m);
		float rhoip = ((3.0f*coeff.x*p+2.0f*coeff.y)*p+coeff.z)*rdr;
//		coeff = tex1Dfetch(rhor_coeff, m);
		float rhojp = ((3.0f*coeff.x*p+2.0f*coeff.y)*p+coeff.z)*rdr;
		
		coeff = tex1Dfetch(z2r_coeff, m);
		float z2p = ((3.0f*coeff.x*p+2.0f*coeff.y)*p+coeff.z)*rdr;
		float z2 = ((coeff.x*p+coeff.y)*p+coeff.z)*p+coeff.w;
		
		float recip = 1.0f/sqrt(rsq);
		float phi = z2*recip;
		float phip = z2p*recip - phi*recip;
		// pos.w -> fp[i], neigh_pos.w -> fp[j]
		float psip = fforce.w*rhojp + neigh_fforce.w*rhoip + phip;
		
		float forcemag_divr = -psip*recip;
		
		// calculate the virial (FLOPS: 3)
		virial += float(1.0/6.0) * rsq * forcemag_divr;

		// add up the force vector components (FLOPS: 7)
		force.x += dx * forcemag_divr;
		force.y += dy * forcemag_divr;
		force.z += dz * forcemag_divr;
		force.w += phi;
		}
		}
	force.w += fforce.z;
	// potential energy per particle must be halved
	force.w *= 0.5f;
	// now that the force calculation is complete, write out the result (MEM TRANSFER: 20 bytes)
	force_data.force[idx_local] = force;
	force_data.virial[idx_local] = virial;
	}

/*! \param force_data Force data on GPU to write forces to
	\param pdata Particle data on the GPU to perform the calculation on
	\param box Box dimensions (in GPU format) to use for periodic boundary conditions
	\param nlist Neighbor list stored on the gpu
	\param d_coeffs A \a coeff_width by \a coeff_width matrix of coefficients indexed by type
		pair i,j. The x-component is lj1 and the y-component is lj2.
	\param coeff_width Width of the \a d_coeffs matrix.
	\param opt More execution options bundled up in a strct
	
	\returns Any error code resulting from the kernel launch
	
	This is just a driver for calcLJForces_kernel, see the documentation for it for more information.
*/
cudaError_t gpu_compute_eam_forces(const gpu_force_data_arrays& force_data, const gpu_pdata_arrays &pdata, const gpu_boxsize &box,
                                   const gpu_nlist_array &nlist, const eam_options& opt, const gpu_eam_spline &spline)
{
  // setup the grid to run the kernel
  dim3 grid( (int)ceil((double)pdata.local_num / (double)opt.block_size), 1, 1);
  dim3 threads(opt.block_size, 1, 1);

	// bind the texture
	pdata_pos_tex.normalized = false;
	pdata_pos_tex.filterMode = cudaFilterModePoint;
	cudaError_t error = cudaBindTexture(0, pdata_pos_tex, pdata.pos, sizeof(float4) * pdata.N);
	if (error != cudaSuccess)
		return error;

  if (IsSet == 0) {
    // bind spline arrays
	  error = cudaBindTexture(0, frho_coeff,  spline.frho_spline,  (opt.nrho+1)*sizeof(float4));
	  if (error != cudaSuccess) return error;
	  error = cudaBindTexture(0, rhor_coeff,  spline.rhor_spline,  (opt.nr+1)*sizeof(float4));
	  if (error != cudaSuccess) return error;
	  error = cudaBindTexture(0, z2r_coeff,   spline.z2r_spline,   (opt.nr+1)*sizeof(float4));
	  if (error != cudaSuccess) return error;
    IsSet = 1;
  }

	// run the kernel
	if (opt.ulf_workaround)	{
		gpu_compute_eam_embed_kernel<true ><<<grid, threads>>>(force_data, pdata, box, nlist, opt.r_cutsq, opt.nr, opt.nrho, opt.rdr, opt.rdrho);
		gpu_compute_eam_force_kernel<true ><<<grid, threads>>>(force_data, pdata, box, nlist, opt.r_cutsq, opt.nr, opt.nrho, opt.rdr, opt.rdrho);
	}	else {
		gpu_compute_eam_embed_kernel<false><<<grid, threads>>>(force_data, pdata, box, nlist, opt.r_cutsq, opt.nr, opt.nrho, opt.rdr, opt.rdrho);
		gpu_compute_eam_force_kernel<false><<<grid, threads>>>(force_data, pdata, box, nlist, opt.r_cutsq, opt.nr, opt.nrho, opt.rdr, opt.rdrho);
	}

	if (!g_gpu_error_checking) {
		return cudaSuccess;
	}	else {
		cudaThreadSynchronize();
		return cudaGetLastError();
	}
}

// vim:syntax=cpp
