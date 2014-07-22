// System includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// Code-specific includes
#include <mex.h>

// Declare texture reference for 2D float texture
texture<float, cudaTextureType2D, cudaReadModeElementType> rf_tex;

// Function prototypes
void cleanup();
void xcorr(int nlhs, mxArray *plhs[], float *rfdata, int nsamps, int nsteps, int nbeams, int srchsz, int kernsz);

// CUDA kernels
#include "xcorr1d_kernel.cu"

////////////////////////////////////////////////////////////////////////////////
// Gateway function to MATLAB (main function)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	if (nrhs != 3)
		mexErrMsgTxt("Wrong number of inputs.\n");
	if (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS && mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
		mexErrMsgTxt("The input rfdata must be of class single or double.\n");
	float *rfdata;

//	if (mxGetNumberOfDimensions(prhs[0]) > 2)
//		mexErrMsgTxt("Only have functionality for 1D xcorr for now. Input must be in (fast time x slow time).\n");

	if (mxGetClassID(prhs[0]) == mxDOUBLE_CLASS) {
		int numel = mxGetNumberOfElements(prhs[0]);
		double *dat = mxGetPr(prhs[0]);
		rfdata = (float *)mxMalloc(sizeof(float)*numel);
		for (int i = 0; i < numel; i++)
			rfdata[i] = (float)dat[i];
	}
	else
		rfdata = (float *)mxGetData(prhs[0]);
	int srchsz = mxGetScalar(prhs[1]);
	int kernsz = mxGetScalar(prhs[2]);

	// Get dimensions of data
	mwSize ndims = mxGetNumberOfDimensions(prhs[0]);
	const mwSize *dims;
	dims = mxGetDimensions(prhs[0]);
	int nsamps = dims[0];
	int nsteps = dims[1];
	int nbeams = 1; for (int i = 2; i < ndims; i++) nbeams *= dims[i];

	// Run cross-correlation peak detector script
	xcorr(nlhs, plhs, rfdata, nsamps, nsteps, nbeams, srchsz, kernsz);
	
	if (mxGetClassID(prhs[0]) == mxDOUBLE_CLASS)
		mxFree(rfdata);
	
	cudaDeviceReset();
}

void xcorr(int nlhs, mxArray *plhs[], float *rfdata, int nsamps, int nsteps, int nbeams, int srchsz, int kernsz) {

	cudaChannelFormatDesc channelDescFLOAT = cudaCreateChannelDesc<float>();
	rf_tex.addressMode[0] = cudaAddressModeClamp;
	rf_tex.addressMode[1] = cudaAddressModeClamp;
	rf_tex.filterMode     = cudaFilterModeLinear;
	rf_tex.normalized     = false;

	float *rf_d, *disp_d, *ccs_d;
	size_t pitch;
	cudaMallocPitch(&rf_d, &pitch, sizeof(float)*nsamps, nsteps*nbeams);
	cudaMalloc((void **)&disp_d,sizeof(float)*nsamps*(nsteps-1)*nbeams);
	cudaMalloc((void **)&ccs_d, sizeof(float)*nsamps*(nsteps-1)*nbeams);
	cudaMemcpy2D(rf_d, pitch, rfdata, sizeof(float)*nsamps,
			sizeof(float)*nsamps, nsteps*nbeams, cudaMemcpyHostToDevice);
	cudaBindTexture2D(NULL, rf_tex, rf_d, channelDescFLOAT, nsamps, nsteps*nbeams, pitch);
	cudaMemset(disp_d, 0, sizeof(float)*nsamps*(nsteps-1)*nbeams);

	dim3 dimB(16, 16, 1);
	dim3 dimG(1, ceil(nsteps/dimB.y)+1, ceil(nbeams/dimB.z)+1);	

	for (int blk = 0; blk < ceil(nsamps/dimB.x)+1; blk++) {
		xcorr1d_kernel<<<dimG, dimB, 0>>>(disp_d, ccs_d, nsamps, nsteps, nbeams, srchsz, kernsz, blk);
		cudaDeviceSynchronize();
	}

	mwSize dims[3] = {nsamps, nsteps-1, nbeams};
	plhs[0] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);

	float *disp_h, *ccs_h;
	disp_h = (float *)mxGetData(plhs[0]);
	ccs_h  = (float *)mxGetData(plhs[1]);
	
	cudaUnbindTexture(rf_tex);
	cudaFree(rf_d);
	cudaMemcpy(disp_h, disp_d,sizeof(float)*nsamps*(nsteps-1)*nbeams, cudaMemcpyDeviceToHost);
	cudaMemcpy(ccs_h,  ccs_d, sizeof(float)*nsamps*(nsteps-1)*nbeams, cudaMemcpyDeviceToHost);

}



















































