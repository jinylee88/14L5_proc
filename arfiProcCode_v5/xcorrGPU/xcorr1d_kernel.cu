__global__ void xcorr1d_kernel(float *disp, float *cc, int nsamps, int nsteps, int nbeams, int srchsz, int kernsz, int blk) {

	int samp = threadIdx.x +    blk    *blockDim.x; 
	int step = threadIdx.y + blockIdx.y*blockDim.y;
	int beam = threadIdx.z + blockIdx.z*blockDim.z;

	int halfk = (kernsz-1)/2;
	int range = (srchsz - kernsz)/2;

	if (samp < nsamps && step < nsteps-1 && beam < nbeams) {

		disp += nsamps*(step + (nsteps-1)*beam);
		
		float x, y, sum_X, sum_Y, sumXX, sumYY, sumXY, tmp, maxrho, d;
		d = 0.0f; maxrho = -1.0f;

		// Build off of previous displacement, if available
		float prev_d = 0.0f;
		if (blk > 0) prev_d = disp[blk*blockDim.x-1];

		// Load post-track data
		sum_Y = sumYY = 0.0f;
		for (int k = -halfk; k <= halfk; k++) {
		
			// Fetch data from texture
			y = tex2D(rf_tex, samp+k + 0.5f, step+1 + beam*nsteps + 0.5f);
			
			// Accumulate sums
			sum_Y += y;
			sumYY += y*y;
		}

		// Loop through search region of reference track and compute CC's
		for (float s = samp-range+prev_d; s <= samp+range+prev_d; s+=0.05f) {
			sum_X = sumXX = sumXY = 0.0f;

			// Loop over kernel
			for (int k = -halfk; k <= halfk; k++) {

				// Fetch data from texture
				x = tex2D(rf_tex, s+k    + 0.5f, step+0 + beam*nsteps + 0.5f);
				y = tex2D(rf_tex, samp+k + 0.5f, step+1 + beam*nsteps + 0.5f);

				// Accumulate sums
				sum_X += x;
				sumXX += x*x;
				sumXY += x*y;

			}
			// Compute normalized cross correlation coefficient
  	        tmp = (kernsz*sumXY - sum_X*sum_Y) *
 			 rsqrt((kernsz*sumXX - sum_X*sum_X) *
 			       (kernsz*sumYY - sum_Y*sum_Y));
			
			if (!isinf(tmp) && tmp > maxrho) {
				maxrho = tmp;
				d      = s-samp;
			}
		}

		// Store outputs
		disp[samp] = d;
		cc  [samp] = maxrho;
	}
} 
/*
__global__ void xcorr1d_kernel(float *disp, float *cc, int nsamps, int nsteps, int srchsz, int kernsz) {

	int samp = threadIdx.x + blockIdx.x*blockDim.x; 
	int step = threadIdx.y + blockIdx.y*blockDim.y;

	int halfk = (kernsz-1)/2;
	int range = (srchsz - kernsz)/2;

	if (samp < nsamps && step < nsteps) {
		
		float x, y, sum_X, sum_Y, sumXX, sumYY, sumXY, tmp, maxrho, d;
		d = 0.0f; maxrho = -1.0f;
		
		for (float s = samp-range; s <= samp+range; s+=0.05f) {

			sum_X = sum_Y = sumXX = sumYY = sumXY = 0.0f;
			for (int k = -halfk; k <= halfk; k++) {
				
				// Fetch data from texture
				x = tex2D(rf_tex, s+k    + 0.5f, 0    + 0.5f);
				y = tex2D(rf_tex, samp+k + 0.5f, step + 0.5f);

				// Accumulate sums
				sum_X += x;
				sum_Y += y;
				sumXX += x*x;
				sumYY += y*y;
				sumXY += x*y;

			}
			// Compute normalized cross correlation coefficient
  	        tmp = (kernsz*sumXY - sum_X*sum_Y) *
 			 rsqrt((kernsz*sumXX - sum_X*sum_X) *
 			       (kernsz*sumYY - sum_Y*sum_Y));
			
			if (!isinf(tmp) && tmp > maxrho) {
				maxrho = tmp;
				d      = s-samp;
			}
		}

		disp[samp + step*nsamps] = d;
		cc  [samp + step*nsamps] = maxrho;

	}
} 
*/
