#include <assert.h>
#include <fftw3.h>
#include <getopt.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>

gsl_rng *rng;

int main(int argc, char **argv) {
	double hurst = 0.5;
	unsigned g = 8;
	int seed = -1;

	int c;
	while ((c = getopt(argc, argv, "h:g:S:")) != -1) {
		switch (c) {
		case 'h':
			hurst = atof(optarg);
			break;
		case 'g':
			g = atoi(optarg);
		case 'S':
			seed = atoi(optarg);
			break;
		default:
			__builtin_unreachable();
		}
	}

	assert(hurst > 0.0 && hurst < 1.0);

	unsigned n;
	assert(g < sizeof(n) * CHAR_BIT);
	n = 1u << g;

	printf("# First passage times of fractional Brownian Motion with drift "
	       "using Davies-Harte method\n"
	       "# H: %g\n"
	       "# System size: %i\n"
	       "# Linear drift: %g\n"
	       "# Fractional drift: %g\n"
	       "# Iterations: %i\n",
	       hurst, n, 0.0 /* FIXME */, 0.0 /* FIXME */, 0 /* FIXME */);

	gsl_rng_env_setup();
	rng = gsl_rng_alloc(gsl_rng_default);
	assert(seed >= 0);
	gsl_rng_set(rng, (unsigned)seed);
	printf("# RNG seed: %i\n", seed);

	/* Compute circulant eigenvalues */

	double *corr = fftw_alloc_real(2*n + 2);
	double prevexp, currexp, nextexp;
	currexp = pow(1.0 / n, 2 * hurst);
	nextexp = 0.0;
	for (unsigned i = 0; i < n; i++) {
		prevexp = currexp;
		currexp = nextexp;
		nextexp = pow((double)(i+1) / (double)n, 2 * hurst);
		corr[i] = corr[2*n-i] = prevexp + nextexp - 2.0*currexp;
	}
	corr[n] = 0.0;

	fftw_complex *eigen = (fftw_complex *)corr;
	fftw_plan eigenplan = fftw_plan_dft_r2c_1d(2*n, corr, eigen, FFTW_ESTIMATE);
	fftw_execute(eigenplan);
	/* corr[] invalid from this point on */
	fftw_destroy_plan(eigenplan);

	fftw_free(eigen);

	fftw_cleanup();
	gsl_rng_free(rng);
	return EXIT_SUCCESS;
}
