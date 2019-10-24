#include <assert.h>
#include <fftw3.h>
#include <getopt.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#define BARRIER 0.1

#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define SQRTPI 1.77245385090551602729816748334114518

static double erfcinv(double x) {
	/* Blair, Edwards, Johnson. "Rational Chebyshev approximations for the
	 * inverse of the error function". Math.Comp. 30 (1976), 827--830.
	 */
	double eta = -log(SQRTPI * x), logeta = log(eta);
	return sqrt(eta - 0.5*logeta + (0.25*logeta - 0.5)/eta);
}

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

static double hurst = 0.5, lindrift = 0.0, fracdrift = 0.0, epsilon = 1e-9;
static gsl_rng *rng;
static size_t size, alloc;
static double *cinv, *times, *values, *gamma_, *g;

static void extend(double *restrict cinv, double *restrict times, double time,
                   size_t size) {
	times[size] = time;

	assert(time >= 0.0);
	for (size_t i = 0; i < size; i++) {
		assert(times[i] >= 0.0);
		gamma_[i] = pow(times[i], 2*hurst) + pow(time, 2*hurst) -
		            pow(fabs(times[i] - time), 2*hurst);
	}
	matvec(g, cinv, gamma_, size);

	double sigsq = 2.0 * pow(time, 2*hurst) - inner(gamma_, g, size);
	assert(sigsq >= 0.0);
	symrk1(cinv, 1.0/sigsq, g, size);
	scavec(cinv + size*(size+1)/2, -1.0/sigsq, g, size);
	cinv[(size+1)*(size+2)/2 - 1] = 1.0/sigsq;
}

bool visitfpt(double *fpt, double ltime, double lpos, double rtime, double rpos,
              unsigned level, double strip) {
	if (level == 0) {
		if (rpos < BARRIER)
			return false;
		*fpt = ltime + (rtime - ltime)*(BARRIER - lpos)/(rpos - lpos);
		return true;
	}
	if (MAX(lpos, rpos) < BARRIER - strip)
		return false;

	if (size + 1 > alloc) {
		alloc *= 2;
		cinv = realloc(cinv, alloc*(alloc+1)/2 * sizeof(*cinv));
		times = realloc(times, alloc * sizeof(*times));
		values = realloc(values, alloc * sizeof(*values));
		gamma_ = realloc(gamma_, alloc * sizeof(*gamma_));
		g = realloc(g, alloc * sizeof(*g));
	}

	double mtime = (ltime + rtime) / 2;
	extend(cinv, times, mtime, size);
	double var = 1.0 / cinv[(size+1)*(size+2)/2 - 1],
	       mean = -var * inner(values, cinv + size*(size+1)/2, size);
	       /* FIXME this doesn’t work with drift */
	double mval = values[size++] = mean + gsl_ran_gaussian_ziggurat(rng, sqrt(var)),
	       mpos = mval + lindrift * mtime + fracdrift * pow(mtime, 2*hurst);

	strip /= pow(2, hurst);
	return visitfpt(fpt, ltime, lpos, mtime, mpos, level - 1, strip) ||
	       visitfpt(fpt, mtime, mpos, rtime, rpos, level - 1, strip);
}

int main(int argc, char **argv) {
	unsigned logn = 8;
	unsigned long seed = 0;
	unsigned iters = 0, levels = 0;

	int c;
	while ((c = getopt(argc, argv, "g:h:m:n:G:I:S:")) != -1) {
		switch (c) {
		case 'g': logn = atoi(optarg); break;
		case 'h': hurst = atof(optarg); break;
		case 'm': lindrift = -atof(optarg); break; /* NB sign */
		case 'n': fracdrift = -atof(optarg); break; /* NB sign */
		case 'E': epsilon = atof(optarg); break;
		case 'G': levels = atoi(optarg); break;
		case 'I': iters = atoi(optarg); break;
		case 'S': seed = atol(optarg); break;
		}
	}

	unsigned n = 1u << logn;
	double dt = 1.0 / n;
	double strip = erfcinv(2*epsilon) * sqrt(4 / pow(2, 2*hurst) - 1) *
	               pow(dt, hurst);

	printf("# Hurst parameter: %g\n"
	       "# Linear drift: %g\n"
	       "# Fractional drift: %g\n"
	       "# Grid size: 2^{%u}\n"
	       "# Effective grid size: 2^{%u}\n"
	       "# Error tolerance: %g\n"
	       "# Iterations: %u\n",
	       hurst, lindrift, fracdrift, logn, logn + levels, epsilon, iters);

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

	alloc = n;
	cinv = malloc(n*(n+1)/2 * sizeof(*cinv));
	times = malloc(n * sizeof(*times));
	values = malloc(n * sizeof(*values));
	gamma_ = malloc(n * sizeof(*gamma_));
	g = malloc(n * sizeof(*g));

	double **cinvs = malloc(n * sizeof(*cinvs));
	for (unsigned i = 0; i < n; i++) {
		cinvs[i] = malloc((i+1)*(i+2)/2 * sizeof(*cinvs[i]));
		if (i > 0) { /* Avoid undefined behaviour */
			memcpy(cinvs[i], cinvs[i-1],
			       i*(i+1)/2 * sizeof(*cinvs[i]));
		}
		extend(cinvs[i], times, (i+1)*dt, i);
	}

	/* Iterate */

	double *noise = fftw_alloc_real(2*n);
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
		for (size = 0; size < n; size++) {
			if (sum + lindrift * size*dt + fracdrift * pow(size*dt, 2*hurst) >= BARRIER)
				break;
			times[size] = (size + 1) * dt;
			values[size] = sum += noise[size];
		}
		memcpy(cinv, cinvs[size-1], size*(size+1)/2 * sizeof(*cinv));

		/* Find first passage */

		double fpt = 1.0, prevtime = 0.0, prevval = 0.0;
		for (unsigned i = 0; i < n; i++) {
			double val = values[i] + lindrift * (i+1)*dt + fracdrift * pow((i+1)*dt, 2*hurst);
			if (visitfpt(&fpt, prevtime, prevval, (i+1)*dt, val,
			             levels, strip))
				break;
			prevtime = (i + 1)*dt;
			prevval  = val;
		}
		printf("%g\n", fpt);
	}

	fftw_free(noise); fftw_destroy_plan(noiseplan);

	free(cinv); free(times); free(values); free(gamma_); free(g);

	for (unsigned i = 0; i < n; i++)
		free(cinvs[i]);
	free(cinvs);

	fftw_free(eigen);

	fftw_cleanup();
	gsl_rng_free(rng);
	return EXIT_SUCCESS;
}
