#define __STDC_WANT_IEC_60559_TYPES_EXT__ 1

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

#if defined(USE_LDBL)
typedef long double real_t;
#define REAL(TOK)   TOK ## l
#define FFTWR(TOK)  fftwl_ ## TOK
#elif defined(USE_FLT128)
typedef _Float128 real_t;
#define REAL(TOK)   TOK ## f128
#define FFTWR(TOK)  fftwq_ ## TOK
#else
typedef double real_t;
#define REAL(TOK)   TOK
#define FFTWR(TOK)  fftw_ ## TOK
#endif

#define fabsr  REAL(fabs)
#define logr   REAL(log)
#define powr   REAL(pow)
#define sqrtr  REAL(sqrt)

#define fftwr_alloc_real    FFTWR(alloc_real)
#define fftwr_cleanup       FFTWR(cleanup)
#define fftwr_destroy_plan  FFTWR(destroy_plan)
#define fftwr_execute       FFTWR(execute)
#define fftwr_free          FFTWR(free)
#define fftwr_plan          FFTWR(plan)
#define fftwr_plan_r2r_1d   FFTWR(plan_r2r_1d)

#define countof(A) (sizeof((A)) / sizeof((A)[0]))
#define MAX(A, B) ((A) > (B) ? (A) : (B))
/* #define PI 3.14159265358979323846264338327950288 */
#define SQRTPI REAL(1.77245385090551602729816748334114518)

#ifdef __GNUC__
#define slower(X) (__builtin_expect((X), 0))
#else
#define slower(X) (X)
#endif

static real_t erfcinv(real_t x) {
	/* Blair, Edwards, Johnson. "Rational Chebyshev approximations for the
	 * inverse of the error function". Math.Comp. 30 (1976), 827--830.
	 */
	real_t eta = -logr(SQRTPI * x), logeta = logr(eta);
	return sqrtr(eta - 0.5*logeta + (0.25*logeta - 0.5)/eta);
}

static real_t inner(real_t const *vec1, real_t const *vec2, size_t size) {
	real_t sum = 0.0;
	for (size_t i = 0; i < size; i++)
		sum += vec1[i] * vec2[i];
	return sum;
}

static void scavec(real_t *out, real_t coe, real_t const *vec, size_t size) {
	for (size_t i = 0; i < size; i++)
		out[i] = coe * vec[i];
}

static void symrk1(real_t *restrict out, real_t coe, real_t const *vec,
                   size_t size) {
	size_t k = 0;
	for (size_t i = 0; i < size; i++) {
		real_t tmp = coe * vec[i];
		for (size_t j = 0; j <= i; j++)
			out[k++] += tmp * vec[j];
	}
}

static void matvec(real_t *restrict out, real_t const *mat, real_t const *vec,
                   size_t size) {
	for (size_t i = 0; i < size; i++) {
		real_t sum = 0.0; size_t j, k;
		for (j = 0, k = i*(i+1)/2; j < i; j++, k++)
			sum += mat[k] * vec[j];
		for (/* j = i, k = (i+1)*(i+2)/2 */; j < size; j++, k += j)
			sum += mat[k] * vec[j];
		out[i] = sum;
	}
}

static void printvec(real_t *vec, size_t size) { /* FIXME debugging only */
	printf("[ ");
	for (size_t row = 0; row < size; row++)
		printf("% .4f ", (double)*vec++);
	printf("]^T\n");
}

static void printmat(real_t *mat, size_t size) { /* FIXME debugging only */
	for (size_t row = 0; row < size; row++) {
		printf("[ ");
		for (size_t col = 0; col <= row; col++)
			printf("% .4f ", (double)*mat++);
		printf("]\n");
	}
}

typedef enum {
	TVARIANCE = 1,
} tflag_t;

static struct {tflag_t flag; char const *name;} tflags[] = {
	{TVARIANCE, "variance"},
};

static tflag_t trace = 0;
static real_t hurst = 0.5, lindrift = 0.0, fracdrift = 0.0, epsilon = 1e-9;
static gsl_rng *rng;
static size_t size, alloc;
static real_t *cinv, *times, *values, *gamma_, *g;

static void extend(real_t *restrict cinv, real_t *restrict times, real_t time,
                   size_t size) {
	times[size] = time;

	assert(time >= 0.0);
	for (size_t i = 0; i < size; i++) {
		assert(times[i] >= 0.0);
		gamma_[i] = powr(times[i], 2*hurst) + powr(time, 2*hurst) -
		            powr(fabsr(times[i] - time), 2*hurst);
	}
	matvec(g, cinv, gamma_, size);

	real_t sigsq = 2.0 * powr(time, 2*hurst) - inner(gamma_, g, size);
	assert(sigsq >= 0.0);
	symrk1(cinv, 1.0/sigsq, g, size);
	scavec(cinv + size*(size+1)/2, -1.0/sigsq, g, size);
	cinv[(size+1)*(size+2)/2 - 1] = 1.0/sigsq;
}

bool visitfpt(real_t *fpt, real_t ltime, real_t lpos, real_t rtime, real_t rpos,
              unsigned level, real_t strip) {
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

	real_t mtime = (ltime + rtime) / 2;
	extend(cinv, times, mtime, size);
	real_t var = 1.0 / cinv[(size+1)*(size+2)/2 - 1],
	       mean = -var * inner(values, cinv + size*(size+1)/2, size);
	real_t mval = values[size++] = mean + gsl_ran_gaussian_ziggurat(rng, sqrt(var)),
	       mpos = mval + lindrift * mtime + fracdrift * pow(mtime, 2*hurst);
	if slower(trace & TVARIANCE)
		printf("# variance %u %g\n", level, (double)var);

	strip /= pow(2, hurst);
	return visitfpt(fpt, ltime, lpos, mtime, mpos, level - 1, strip) ||
	       visitfpt(fpt, mtime, mpos, rtime, rpos, level - 1, strip);
}

int main(int argc, char **argv) {
	unsigned logn = 8;
	unsigned long seed = 0;
	unsigned iters = 0, levels = 0;

	int c;
	while ((c = getopt(argc, argv, "g:h:m:n:t:G:I:S:")) != -1) {
		switch (c) {
		case 'g': logn = atoi(optarg); break;
		case 'h': hurst = atof(optarg); break;
		case 'm': lindrift = -atof(optarg); break; /* NB sign */
		case 'n': fracdrift = -atof(optarg); break; /* NB sign */
		case 'E': epsilon = atof(optarg); break;
		case 'G': levels = atoi(optarg); break;
		case 'I': iters = atoi(optarg); break;
		case 'S': seed = atol(optarg); break;

		case 't':
			for (size_t i = 0; i < countof(tflags); i++) {
				if (!strcmp(optarg, tflags[i].name))
					trace |= tflags[i].flag;
			}
			break;
		}
	}

	unsigned n = 1u << logn;
	real_t dt = 1.0 / (real_t)n;
	real_t strip = erfcinv(2*epsilon) * sqrtr(4 / powr(2, 2*hurst) - 1) *
	               powr(dt, hurst);

	gsl_rng_env_setup();
	rng = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(rng, (unsigned)seed);

	printf("# Hurst parameter: %.17e\n"
	       "# Linear drift: %.17e\n"
	       "# Fractional drift: %.17e\n"
	       "# Log of grid size: %u\n"
	       "# Levels to descend: %u\n"
	       "# Error tolerance: %.17e\n"
	       "# Iterations: %u\n"
	       "# RNG seed: %lu\n"
	       "#\n",
	       (double)hurst, (double)lindrift, (double)fracdrift, logn,
	       levels, (double)epsilon, iters, seed);

	/* Compute circulant eigenvalues */

	real_t *eigen = fftwr_alloc_real(n + 1);
	real_t prevexp, currexp, nextexp;
	currexp = pow(1.0 / n, 2 * hurst);
	nextexp = 0.0;
	for (unsigned i = 0; i < n; i++) {
		prevexp = currexp;
		currexp = nextexp;
		nextexp = powr((real_t)(i+1) / (real_t)n, 2*hurst);
		eigen[i] = prevexp + nextexp - 2.0*currexp;
	}
	eigen[n] = 0.0;

	fftwr_plan eigenplan = fftwr_plan_r2r_1d(n + 1, eigen, eigen,
	                                         FFTW_REDFT00, FFTW_ESTIMATE);
	fftwr_execute(eigenplan);
	fftwr_destroy_plan(eigenplan);

	/* Compute inverse correlation matrices */

	alloc = n;
	cinv = malloc(n*(n+1)/2 * sizeof(*cinv));
	times = malloc(n * sizeof(*times));
	values = malloc(n * sizeof(*values));
	gamma_ = malloc(n * sizeof(*gamma_));
	g = malloc(n * sizeof(*g));

	real_t **cinvs = malloc(n * sizeof(*cinvs));
	for (unsigned i = 0; i < n; i++) {
		cinvs[i] = malloc((i+1)*(i+2)/2 * sizeof(*cinvs[i]));
		if (i > 0) { /* Avoid undefined behaviour */
			memcpy(cinvs[i], cinvs[i-1],
			       i*(i+1)/2 * sizeof(*cinvs[i]));
		}
		extend(cinvs[i], times, (i+1)*dt, i);
	}

	/* Iterate */

	real_t *noise = fftwr_alloc_real(2*n);
	fftwr_plan noiseplan = fftwr_plan_r2r_1d(2*n, noise, noise,
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
		fftwr_execute(noiseplan);

		/* Integrate */

		/* FIXME Kahan summation ? */
		real_t sum = 0.0;
		for (size = 0; size < n; size++) {
			if (sum + lindrift * size*dt + fracdrift * pow(size*dt, 2*hurst) >= BARRIER)
				break;
			times[size] = (size + 1) * dt;
			values[size] = sum += noise[size];
		}
		memcpy(cinv, cinvs[size-1], size*(size+1)/2 * sizeof(*cinv));

		/* Find first passage */

		real_t fpt = 1.0, prevtime = 0.0, prevpos = 0.0;
		for (unsigned i = 0; i < n; i++) {
			real_t time = (i + 1) * dt,
			       pos = values[i] + lindrift * (i+1)*dt +
			             fracdrift * powr((i+1)*dt, 2*hurst);
			if (visitfpt(&fpt, prevtime, prevpos, time, pos,
			             levels, strip))
				break;
			prevtime = time;
			prevpos  = pos;
		}
		printf("%g\n", (double)fpt);
	}

	fftwr_free(noise); fftwr_destroy_plan(noiseplan);

	free(cinv); free(times); free(values); free(gamma_); free(g);

	for (unsigned i = 0; i < n; i++)
		free(cinvs[i]);
	free(cinvs);

	fftwr_free(eigen);

	fftwr_cleanup();
	gsl_rng_free(rng);
	return EXIT_SUCCESS;
}
