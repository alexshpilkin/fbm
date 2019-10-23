#include <assert.h>
#include <fftw3.h>
#include <getopt.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define BARRIER 0.1

static double inner(double const *vec1, double const *vec2, size_t size) {
	double sum = 0.0;
	for (size_t i = 0; i < size; i++)
		sum += vec1[i] * vec2[i];
	return sum;
}

static void scavec(double *out, double coe, double const *vec, size_t size) {
	for (size_t i = 0; i < size; i++)
		out[i] = coe * vec[i];
}

static void symrk1(double *restrict out, double coe, double const *vec,
                   size_t size) {
	size_t k = 0;
	for (size_t i = 0; i < size; i++) {
		double tmp = coe * vec[i];
		for (size_t j = 0; j <= i; j++)
			out[k++] += tmp * vec[j];
	}
}

static void matvec(double *restrict out, double const *mat, double const *vec,
                   size_t size) {
	for (size_t i = 0; i < size; i++) {
		double sum = 0.0; size_t j, k;
		for (j = 0, k = i*(i+1)/2; j < i; j++, k++)
			sum += mat[k] * vec[j];
		for (/* j = i, k = (i+1)*(i+2)/2 */; j < size; j++, k += j)
			sum += mat[k] * vec[j];
		out[i] = sum;
	}
}

static void printvec(double *vec, size_t size) { /* FIXME debugging only */
	printf("[ ");
	for (size_t row = 0; row < size; row++)
		printf("% .4f ", *vec++);
	printf("]^T\n");
}

static void printmat(double *mat, size_t size) { /* FIXME debugging only */
	for (size_t row = 0; row < size; row++) {
		printf("[ ");
		for (size_t col = 0; col <= row; col++)
			printf("% .4f ", *mat++);
		printf("]\n");
	}
}

static gsl_rng *rng;
static double hurst = 0.5, lindrift = 0.0, fracdrift = 0.0;

static void addpoint(double *restrict cinv, double *restrict times, double time,
                     size_t size) {
	/* FIXME Inner-loop allocation */
	double *gamma = malloc(size * sizeof(*gamma)),
	       *g = malloc(size * sizeof(*g));

	times[size] = time;

	assert(time >= 0.0);
	for (size_t i = 0; i < size; i++) {
		assert(times[i] >= 0.0);
		gamma[i] = pow(times[i], 2*hurst) + pow(time, 2*hurst) -
		           pow(fabs(times[i] - time), 2*hurst);
	}
	matvec(g, cinv, gamma, size);

	double sigsq = 2.0 * pow(time, 2*hurst) - inner(gamma, g, size);
	symrk1(cinv, 1.0/sigsq, g, size);
	scavec(cinv + size*(size+1)/2, -1.0/sigsq, g, size);
	cinv[(size+1)*(size+2)/2 - 1] = 1.0/sigsq;

	free(gamma); free(g);
}

int main(int argc, char **argv) {
	unsigned logn = 8;
	unsigned long seed = 0;
	unsigned iters = 0;

	int c;
	while ((c = getopt(argc, argv, "h:g:m:n:I:S:")) != -1) {
		switch (c) {
		case 'h': hurst = atof(optarg); break;
		case 'g': logn = atoi(optarg); break;
		case 'm': lindrift = -atof(optarg); break; /* NB sign */
		case 'n': fracdrift = -atof(optarg); break; /* NB sign */
		case 'I': iters = atoi(optarg); break;
		case 'S': seed = atol(optarg); break;
		}
	}

	unsigned n = 1u << logn;
	double dt = 1.0 / n;

	printf("# First passage times of fractional Brownian Motion with drift "
	       "using Davies-Harte method\n"
	       "# H: %g\n"
	       "# System size: %i\n"
	       "# Linear drift: %g\n"
	       "# Fractional drift: %g\n"
	       "# Iterations: %i\n",
	       hurst, n, lindrift, fracdrift, iters);

	gsl_rng_env_setup();
	rng = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(rng, (unsigned)seed);
	printf("# RNG seed: %lu\n", seed);

	/* Compute circulant eigenvalues */

	double *eigen = fftw_alloc_real(n + 1);
	double prevexp, currexp, nextexp;
	currexp = pow(1.0 / n, 2 * hurst);
	nextexp = 0.0;
	for (unsigned i = 0; i < n; i++) {
		prevexp = currexp;
		currexp = nextexp;
		nextexp = pow((double)(i+1) / (double)n, 2 * hurst);
		eigen[i] = prevexp + nextexp - 2.0*currexp;
	}
	eigen[n] = 0.0;

	fftw_plan eigenplan = fftw_plan_r2r_1d(n + 1, eigen, eigen,
	                                       FFTW_REDFT00, FFTW_ESTIMATE);
	fftw_execute(eigenplan);
	fftw_destroy_plan(eigenplan);

	/* Compute inverse correlation matrices */

	double **cinvs = malloc(n * sizeof(*cinvs)),
	        *times = malloc(n * sizeof(*times));

	for (unsigned i = 0; i < n; i++) {
		cinvs[i] = malloc((i+1)*(i+2)/2 * sizeof(*cinvs[i]));
		if (i > 0) { /* Avoid undefined behaviour */
			memcpy(cinvs[i], cinvs[i-1],
			       i*(i+1)/2 * sizeof(*cinvs[i]));
		}
		addpoint(cinvs[i], times, (i+1)*dt, i);
	}

	/* Iterate */

	double *noise = fftw_alloc_real(2*n), *trace = noise;
	fftw_plan noiseplan = fftw_plan_r2r_1d(2*n, noise, noise,
	                                       FFTW_HC2R, FFTW_ESTIMATE);

	for (unsigned iter = 0; iter < iters; iter++) {
		/* Generate noise */

		gsl_ran_gaussian_ziggurat(rng, 1.0); /* FIXME Sync with  */
		gsl_ran_gaussian_ziggurat(rng, 1.0); /* Walter’s version */
		for (unsigned i = 1; i < n; i++) {
			noise[i] = sqrt(0.25 * eigen[i] / n) *
			           gsl_ran_gaussian_ziggurat(rng, 1.0);
			noise[2*n-i] = sqrt(0.25 * eigen[i] / n) *
			               gsl_ran_gaussian_ziggurat(rng, 1.0);
		}
		noise[0] = sqrt(0.5 * eigen[0] / n) *
		           gsl_ran_gaussian_ziggurat(rng, 1.0);
		noise[n] = sqrt(0.5 * eigen[n] / n) *
		           gsl_ran_gaussian_ziggurat(rng, 1.0);
		fftw_execute(noiseplan);

		/* Integrate */

		/* FIXME Kahan summation ? */
		double sum = 0.0;
		for (unsigned i = 0; i < n; i++) {
			double a = noise[i] + dt * (lindrift + fracdrift * pow((i + 0.5)*dt, 2.0*hurst - 1.0));
			trace[i] = sum += a;
		}

		/* Find first passage */

		unsigned i;
		double tprev = 0.0;
		for (i = 0; i < n; i++) {
			if (trace[i] > BARRIER)
				break;
			tprev = trace[i];
		}
		printf("%g\n", i >= n ? 1.0 : dt * (i + (BARRIER - tprev) / (trace[i] - tprev)));
	}

	fftw_free(noise);
	fftw_destroy_plan(noiseplan);

	for (unsigned i = 0; i < n; i++)
		free(cinvs[i]);
	free(cinvs); free(times);

	fftw_free(eigen);

	fftw_cleanup();
	gsl_rng_free(rng);
	return EXIT_SUCCESS;
}
