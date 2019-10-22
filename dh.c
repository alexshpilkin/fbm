#include <assert.h>
#include <fftw3.h>
#include <getopt.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>

#define BARRIER 0.1

gsl_rng *rng;

int main(int argc, char **argv) {
	double hurst = 0.5, lindrift = 0.0, fracdrift = 0.0;
	unsigned logn = 8;
	int seed = -1;
	int iters = 0;

	int c;
	while ((c = getopt(argc, argv, "h:g:m:n:I:S:")) != -1) {
		switch (c) {
		case 'h':
			hurst = atof(optarg);
			break;
		case 'g':
			logn = atoi(optarg);
			break;
		case 'm':
			lindrift = -atof(optarg); /* NB sign */
			break;
		case 'n':
			fracdrift = -atof(optarg); /* NB sign */
			break;
		case 'I':
			iters = atoi(optarg);
			break;
		case 'S':
			seed = atoi(optarg);
			break;
		default:
			__builtin_unreachable();
		}
	}

	assert(hurst > 0.0 && hurst < 1.0);

	unsigned n;
	assert(logn < sizeof(n) * CHAR_BIT);
	n = 1u << logn;
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
	assert(seed >= 0);
	gsl_rng_set(rng, (unsigned)seed);
	printf("# RNG seed: %i\n", seed);

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

	/* Iterate */

	fftw_complex *temp = fftw_alloc_complex(n + 1);
	double *noise = (double *)temp, *trace = noise;
	fftw_plan noiseplan = fftw_plan_dft_c2r_1d(2 * n, temp, noise,
	                                           FFTW_ESTIMATE);

	for (int iter = 0; iter < iters; iter++) {
		/* Generate noise */

		for (unsigned i = 0 /* FIXME 1 */; i < n; i++) {
			temp[i][0] = sqrt(0.25 * eigen[i] / n) *
				     gsl_ran_gaussian_ziggurat(rng, 1.0);
			temp[i][1] = sqrt(0.25 * eigen[i] / n) *
				     gsl_ran_gaussian_ziggurat(rng, 1.0);
		}
		temp[0][0] = sqrt(0.5 * eigen[0] / n) * gsl_ran_gaussian_ziggurat(rng, 1.0);
		temp[n][0] = sqrt(0.5 * eigen[n] / n) * gsl_ran_gaussian_ziggurat(rng, 1.0);
		temp[0][1] = temp[n][1] = 0.0;
		fftw_execute(noiseplan);

		/* Integrate */

		/* FIXME Kahan summation ? */
		double sum = 0.0;
		for (unsigned i = 0; i < n; i++) {
			double a = noise[i] + dt * (lindrift + fracdrift * pow((i + 0.5)*dt, 2.0*hurst - 1.0));
			trace[i] = sum;
			sum += a;
		}
		trace[n] = sum;

		/* Find first passage */

		unsigned i;
		for (i = 1; i <= n; i++) {
			if (trace[i] > BARRIER)
				break;
		}
		printf("%g\n", i > n ? 1.0 : dt * (i - 1 + (BARRIER - trace[i-1]) / (trace[i] - trace[i-1])));
	}

	fftw_free(temp);
	fftw_destroy_plan(noiseplan);
	fftw_free(eigen);

	fftw_cleanup();
	gsl_rng_free(rng);
	return EXIT_SUCCESS;
}
