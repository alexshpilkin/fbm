#include <assert.h>
#include <fftw3.h>
#include <getopt.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <math.h>
#ifdef USE_QUAD
#include <quadmath.h>
#endif
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#if defined(USE_LDBL)
typedef long double real_t;
#define REAL(TOK)   TOK ## l
#define FFTWR(TOK)  fftwl_ ## TOK
#elif defined(USE_QUAD)
typedef __float128 real_t;
#define REAL(TOK)   TOK ## q
#define FFTWR(TOK)  fftwq_ ## TOK
#else
typedef double real_t;
#define REAL(TOK)   TOK
#define FFTWR(TOK)  fftw_ ## TOK
#endif

#define fabsr   REAL(fabs)
#define log2r   REAL(log2)
#define logr    REAL(log)
#define powr    REAL(pow)
#define lrintr  REAL(lrint)
#define sqrtr   REAL(sqrt)

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
#define faster(X) (__builtin_expect(!!(X), 1))
#define slower(X) (__builtin_expect((X), 0))
#else
#define faster(X) (!!(X))
#define slower(X) (X)
#endif

static real_t erfcinv(real_t x) {
	/* Blair, Edwards, Johnson. "Rational Chebyshev approximations for the
	 * inverse of the error function". Math.Comp. 30 (1976), 827--830.
	 */
	real_t eta = -logr(SQRTPI * x), logeta = logr(eta);
	return sqrtr(eta - 0.5*logeta + (0.25*logeta - 0.5)/eta);
}

#ifndef DO_PHONEBOOK
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
#endif /* !DO_PHONEBOOK */

typedef enum {
	TBISECTS = 1,
	TBRIDGES = 2,
	TSTRIP = 4,
#ifndef DO_PHONEBOOK
	TVARIANCE = 8,
#endif /* !DO_PHONEBOOK */
} tflag_t;

static struct {tflag_t flag; char const *name;} tflags[] = {
	{TBISECTS, "bisects"},
	{TBRIDGES, "bridges"},
	{TSTRIP, "strip"},
#ifndef DO_PHONEBOOK
	{TVARIANCE, "variance"},
#endif /* !DO_PHONEBOOK */
};

typedef struct {
	real_t ltime, lpos, rtime, rpos;
} bridge_t;

static tflag_t trace = 0;
static real_t hurst = 0.5, lindrift = 0.0, fracdrift = 0.0, epsilon = 1e-9,
              stripfac;
#ifdef DO_FPT
static real_t barrier = 1.0;
#endif /* DO_FPT */
static gsl_rng *rng;
static size_t size, reserved;
#ifdef DO_PHONEBOOK
static real_t *pbtimes, *pbvalues;
#else /* if !DO_PHONEBOOK */
#define pbtimes times
#define pbvalues values
static real_t *cinv, *gamma_, *g;
#endif /* !DO_PHONEBOOK */
static real_t *times, *values;
#ifdef DO_MAX
static bridge_t *queue;
static size_t top;
#endif /* DO_MAX */
static unsigned *bisects;

#ifdef DO_MAX
static int compare(void const *lhs_, void const *rhs_) {
	bridge_t const *lhs = lhs_, *rhs = rhs_;
	real_t value = MAX(lhs->lpos, lhs->rpos) - MAX(rhs->lpos, rhs->rpos);
	/* Sort from larger to smaller maxima */
	if (value > 0.0) return -1;
	if (value < 0.0) return  1;
	return 0;
}
#endif /* DO_MAX */

#ifndef DO_PHONEBOOK
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
	symrk1(cinv, 1.0/sigsq, g, size);
	scavec(cinv + size*(size+1)/2, -1.0/sigsq, g, size);
	cinv[(size+1)*(size+2)/2 - 1] = 1.0/sigsq;
}
#endif /* !DO_PHONEBOOK */

static real_t sample(real_t time, unsigned level) {
	assert(level > 0);

	if (size + 1 > reserved) {
		reserved *= 2;
#ifndef DO_PHONEBOOK
		cinv = realloc(cinv, reserved*(reserved+1)/2 * sizeof(*cinv));
		gamma_ = realloc(gamma_, reserved * sizeof(*gamma_));
		g = realloc(g, reserved * sizeof(*g));
#endif /* !DO_PHONEBOOK */
		times = realloc(times, reserved * sizeof(*times));
		values = realloc(values, reserved * sizeof(*values));
#ifdef DO_MAX
		queue = realloc(queue, 2*reserved * sizeof(*queue));
#endif /* DO_MAX */
	}

#ifdef DO_PHONEBOOK
	size_t i;
	times[size] = time;
	for (i = 0; /* true */; i++) {
		if (pbtimes[i] == time) break;
		assert(pbtimes[i] < time);
	}
	real_t value = values[size++] = pbvalues[i];
#else /* !DO_PHONEBOOK */
	extend(cinv, times, time, size);
	real_t var = 1.0 / cinv[(size+1)*(size+2)/2 - 1],
	       mean = -var * inner(values, cinv + size*(size+1)/2, size);
	real_t value = values[size++] =
	    var >= 0 ? mean + gsl_ran_gaussian_ziggurat(rng, sqrtr(var)) : NAN;
#endif /* !DO_PHONEBOOK */
	real_t pos = value + lindrift * time + fracdrift * powr(time, 2*hurst);

#ifndef DO_PHONEBOOK
	if slower(trace & TVARIANCE)
		printf("# variance %u %.17e\n", level, (double)var);
#endif /* !DO_PHONEBOOK */
	if slower(trace & TBISECTS)
		bisects[level-1]++;
	return pos;
}

#ifdef DO_FPT
static bool visitfpt(real_t *fpt, real_t ltime, real_t lpos, real_t rtime,
                     real_t rpos, unsigned level, real_t strip) {
	if slower(trace & TBRIDGES) {
		printf("# bridge %u %.17e %.17e %.17e %.17e\n",
		       level, (double)ltime, (double)lpos, (double)rtime,
		       (double)rpos);
	}
	if slower(trace & TSTRIP) {
		printf("# strip %.17e %.17e\n",
		       (double)barrier, (double)strip);
	}
	if (level == 0) {
		if (rpos < barrier)
			return false;
		*fpt = ltime + (rtime - ltime)*(barrier - lpos)/(rpos - lpos);
		return true;
	}
	if (MAX(lpos, rpos) < barrier - strip)
		return false;

	real_t mtime = (ltime + rtime)/2, mpos = sample(mtime, level);
	if (isnan(mpos)) {
		*fpt = mpos;
		return true;
	}
	strip *= stripfac;
	return visitfpt(fpt, ltime, lpos, mtime, mpos, level-1, strip) ||
	       visitfpt(fpt, mtime, mpos, rtime, rpos, level-1, strip);
}
#endif /* DO_FPT */

#ifdef DO_MAX
static void visitmax(real_t *maxtime, real_t *maxpos, real_t ltime, real_t lpos,
                     real_t rtime, real_t rpos, unsigned level, real_t strip) {
	if slower(trace & TBRIDGES) {
		printf("# bridge %u %.17e %.17e %.17e %.17e\n",
		       level, (double)ltime, (double)lpos, (double)rtime,
		       (double)rpos);
	}
	if slower(trace & TSTRIP) {
		printf("# strip %.17e %.17e\n",
		       (double)*maxpos, (double)strip);
	}
	if (MAX(lpos, rpos) < *maxpos - strip)
		return;
	if (lpos > *maxpos) {
		*maxtime = ltime;
		*maxpos = lpos;
	}
	if (rpos > *maxpos) {
		*maxtime = rtime;
		*maxpos = rpos;
	}
	if (level == 0)
		return;

	real_t mtime = (ltime + rtime)/2, mpos = sample(mtime, level);
	if (isnan(mpos)) {
		*maxtime = mpos;
		*maxpos = mpos;
	}
	bridge_t left  = {ltime, lpos, mtime, mpos},
	         right = {mtime, mpos, rtime, rpos};
	assert(top + 2 <= 2 * size);
	queue[top++] = left;
	queue[top++] = right;
}
#endif /* DO_MAX */

#ifdef DO_PHONEBOOK
static unsigned findlevel(real_t time) {
	real_t left = 0.0, right = 1.0, pbdt = pbtimes[0];
	for (size_t i = 0; i < size; i++) {
		assert(times[i] != time);
		if (times[i] < time && time - times[i] < time - left)
			left = times[i];
		if (times[i] > time && times[i] - time < right - time)
			right = times[i];
	}
	return (unsigned)lrintr(log2r((right - left)/pbdt));
}
#endif /* DO_PHONEBOOK */

int main(int argc, char **argv) {
	unsigned logn = 8, levels = 0, iters = 0;
	bool reseed = false;
	unsigned long seed = 0;

	int c;
	while ((c = getopt(argc, argv, "b:g:h:m:n:t:E:G:I:RS:")) != -1) {
		switch (c) {
#ifdef DO_FPT
		case 'b': barrier = atof(optarg); break;
#endif /* DO_FPT */
		case 'g': logn = atoi(optarg); break;
		case 'h': hurst = atof(optarg); break;
		case 'm': lindrift = -atof(optarg); break; /* NB sign */
		case 'n': fracdrift = -atof(optarg); break; /* NB sign */
		case 'E': epsilon = atof(optarg); break;
		case 'G': levels = atoi(optarg); break;
		case 'I': iters = atoi(optarg); break;
		case 'R': reseed = true; break;
		case 'S': seed = atol(optarg); break;

		case 't':
			for (size_t i = 0; i < countof(tflags); i++) {
				if (!strcmp(optarg, tflags[i].name))
					trace |= tflags[i].flag;
			}
			break;
		}
	}

#ifdef DO_PHONEBOOK
	size_t pbn = (size_t)1 << (logn + levels);
	real_t pbdt = 1.0 / (real_t)pbn;
#else /* !DO_PHONEBOOK */
#define pbn n
#define pbdt dt
#endif /* !DO_PHONEBOOK */
	size_t n = (size_t)1 << logn;
	real_t dt = 1.0 / (real_t)n;
	real_t strip = erfcinv(2*epsilon) * sqrtr(4 / powr(2, 2*hurst) - 1) *
	               powr(dt, hurst);
	stripfac = powr(2, -hurst);

	gsl_rng_env_setup();
	gsl_rng *seedrng;
	rng = gsl_rng_alloc(gsl_rng_default);
	if (reseed) {
		seedrng = gsl_rng_alloc(gsl_rng_default);
		gsl_rng_set(seedrng, seed);
	} else {
		gsl_rng_set(rng, seed);
	}

	printf("# Program: %s\n"
	       "# Hurst parameter: %.17e\n"
	       "# Linear drift: %.17e\n"
	       "# Fractional drift: %.17e\n"
#ifdef DO_FPT
	       "# Barrier height: %.17e\n"
#endif /* DO_FPT */
	       "# Log of grid size: %u\n"
	       "# Levels to descend: %u\n"
	       "# Error tolerance: %.17e\n"
	       "# Iterations: %u\n"
	       "# RNG seed: %lu\n"
	       "# Reseeding: %s\n"
	       "#\n",
#ifdef DO_FPT
	       "fpt"
#endif /* DO_FPT */
#ifdef DO_MAX
	       "max"
#endif /* DO_MAX */
#ifdef DO_PHONEBOOK
	       "-pb"
#endif /* DO_PHONEBOOK */
	       , (double)hurst, (double)lindrift, (double)fracdrift,
#ifdef DO_FPT
	       (double)barrier,
#endif /* DO_FPT */
	       logn, levels, (double)epsilon, iters, seed,
	       reseed ? "yes" : "no");
	fflush(stdout);

	/* Compute circulant eigenvalues */

	real_t *eigen = fftwr_alloc_real(pbn + 1);
	real_t prevexp, currexp, nextexp;
	currexp = powr(1.0 / pbn, 2 * hurst);
	nextexp = 0.0;
	for (size_t i = 0; i < pbn; i++) {
		prevexp = currexp;
		currexp = nextexp;
		nextexp = powr((real_t)(i+1) / (real_t)pbn, 2*hurst);
		eigen[i] = prevexp + nextexp - 2.0*currexp;
	}
	eigen[pbn] = 0.0;

	fftwr_plan eigenplan = fftwr_plan_r2r_1d(pbn + 1, eigen, eigen,
	                                         FFTW_REDFT00, FFTW_ESTIMATE);
	fftwr_execute(eigenplan);
	fftwr_destroy_plan(eigenplan);

	/* Compute inverse correlation matrices */

	reserved = n;
#ifdef DO_PHONEBOOK
	pbtimes = malloc(pbn * sizeof(*pbtimes));
	pbvalues = malloc(pbn * sizeof(*pbvalues));
#else /* if !DO_PHONEBOOK */
	cinv = malloc(n*(n+1)/2 * sizeof(*cinv));
	gamma_ = malloc(n * sizeof(*gamma_));
	g = malloc(n * sizeof(*g));
#endif /* !DO_PHONEBOOK */
	times = malloc(n * sizeof(*times));
	values = malloc(n * sizeof(*values));
#ifdef DO_MAX
	queue = malloc(2*n * sizeof(*queue));
#endif /* DO_MAX */

#ifndef DO_PHONEBOOK
	real_t **cinvs = NULL;
	if faster(levels > 0) {
		cinvs = malloc(n * sizeof(*cinvs));
		for (size_t i = 0; i < n; i++) {
			cinvs[i] = malloc((i+1)*(i+2)/2 * sizeof(*cinvs[i]));
			if (i > 0) { /* Avoid undefined behaviour */
				memcpy(cinvs[i], cinvs[i-1],
				       i*(i+1)/2 * sizeof(*cinvs[i]));
			}
			extend(cinvs[i], times, (i+1)*dt, i);
		}
	}
#endif /* !DO_PHONEBOOK */

	/* Iterate */

	real_t *noise = fftwr_alloc_real(2*pbn);
	fftwr_plan noiseplan = fftwr_plan_r2r_1d(2*pbn, noise, noise,
	                                         FFTW_HC2R, FFTW_ESTIMATE);

	if (trace & TBISECTS)
		bisects = malloc(levels * sizeof(*bisects));

	for (unsigned iter = 0; iter < iters; iter++) {
		/* Reseed */

		if (reseed)
			gsl_rng_set(rng, gsl_rng_get(seedrng));

		/* Generate noise */

		for (size_t i = 1; i < pbn; i++) {
			noise[i] = sqrtr(0.25 * eigen[i] / pbn) *
			           gsl_ran_gaussian_ziggurat(rng, 1.0);
			noise[2*pbn-i] = sqrtr(0.25 * eigen[i] / pbn) *
			                 gsl_ran_gaussian_ziggurat(rng, 1.0);
		}
		noise[0] = sqrtr(0.5 * eigen[0] / pbn) *
		           gsl_ran_gaussian_ziggurat(rng, 1.0);
		noise[pbn] = sqrtr(0.5 * eigen[n] / pbn) *
		             gsl_ran_gaussian_ziggurat(rng, 1.0);
		fftwr_execute(noiseplan);

		/* Integrate */

		real_t sum = 0.0;
		for (size_t i = 0; i < pbn; i++) {
			pbtimes[i] = (i + 1) * pbdt;
			pbvalues[i] = sum += noise[i];
		}
#ifdef DO_PHONEBOOK
		for (size_t i = 0; i < n; i++) {
			size_t j = (1 << levels) * (i+1) - 1;
			times[i] = pbtimes[j];
			values[i] = pbvalues[j];
		}
#endif /* DO_PHONEBOOK */
		for (size = 1; size < n; size++) {
#ifdef DO_FPT
			if (values[size-1] + lindrift * size*dt + fracdrift * powr(size*dt, 2*hurst) >= barrier)
				break;
#endif /* DO_FPT */
		}
#ifndef DO_PHONEBOOK
		if faster(levels > 0) {
			memcpy(cinv, cinvs[size-1],
			       size*(size+1)/2 * sizeof(*cinv));
		}
#endif /* !DO_PHONEBOOK */

		/* Find first passage or maximum */

		if slower(trace & TBISECTS)
			memset(bisects, 0, levels * sizeof(*bisects));

#ifdef DO_FPT
		real_t fpt = 1.0;
#endif /* DO_FPT */
		real_t prevtime = 0.0, prevpos = 0.0;
		for (size_t i = 0; i < n; i++) {
			real_t time = (i + 1) * dt,
			       pos = values[i] + lindrift * time +
			             fracdrift * powr(time, 2*hurst);
#ifdef DO_FPT
			if (visitfpt(&fpt, prevtime, prevpos, time, pos,
			             levels, strip))
				break;
#endif /* DO_FPT */
#ifdef DO_MAX
			queue[i].ltime = prevtime;
			queue[i].lpos  = prevpos;
			queue[i].rtime = time;
			queue[i].rpos  = pos;
#endif /* DO_MAX */
			prevtime = time;
			prevpos  = pos;
		}

#ifdef DO_MAX
		real_t maxtime = 0.0, maxpos = 0.0;
		size_t bottom = 0;
		top = size;
		real_t levelstrip = strip;
		for (unsigned level = levels+1; level > 0; level--) {
			size_t prevtop = top;
			qsort(queue + bottom, prevtop - bottom, sizeof(*queue),
			      compare);
			for (size_t i = bottom; i < prevtop; i++) {
				if (isnan(maxpos))
					break;
				visitmax(&maxtime, &maxpos,
				         queue[i].ltime, queue[i].lpos,
				         queue[i].rtime, queue[i].rpos,
				         level-1, levelstrip);
			}
			bottom = prevtop;
			levelstrip *= stripfac;
		}
#endif /* DO_MAX */

		if slower(trace & TBISECTS) {
			printf("# bisects");
			for (unsigned i = levels; i > 0; i--)
				printf(" %u", bisects[i-1]);
			printf("\n");
		}

#ifdef DO_PHONEBOOK
#ifdef DO_FPT
		real_t pbfpt = 1.0;
		prevtime = 0.0; prevpos = 0.0;
#endif /* DO_FPT */
#ifdef DO_MAX
		real_t pbmaxtime = 0.0, pbmaxpos = 0.0;
#endif /* DO_MAX */
		for (size_t i = 0; i < pbn; i++) {
			real_t time = (i + 1) * pbdt,
			       pos = pbvalues[i] + lindrift * time +
			             fracdrift * powr(time, 2*hurst);
#ifdef DO_FPT
			if (pos >= barrier) {
				pbfpt = prevtime + (time - prevtime) *
				                   (barrier - prevpos) /
				                   (pos - prevpos);
				break;
			}
			prevtime = time; prevpos = pos;
#endif /* DO_FPT */
#ifdef DO_MAX
			if (pbmaxpos < pos) {
				pbmaxtime = time;
				pbmaxpos = pos;
			}
#endif /* DO_MAX */
		}
#ifdef DO_FPT
		if (fpt != pbfpt) {
			assert(fpt > pbfpt);
			printf("# error %u %.17e %.17e\n",
			       findlevel(pbfpt), (double)pbfpt, (double)fpt);
		}
#endif /* DO_FPT */
#ifdef DO_MAX
		if (maxpos != pbmaxpos) {
			assert(maxpos < pbmaxpos);
			printf("# error %u %.17e %.17e %.17e %.17e\n",
			       findlevel(pbmaxtime), (double)pbmaxtime,
			       (double)pbmaxpos, (double)maxtime,
			       (double)maxpos);
		}
#endif /* DO_MAX */
#endif /* DO_PHONEBOOK */

#ifdef DO_FPT
		printf("%.17e\n", (double)fpt);
#endif /* DO_FPT */
#ifdef DO_MAX
		printf("%.17e %.17e\n", (double)maxtime, (double)maxpos);
#endif /* DO_MAX */
		fflush(stdout);
	}

	if (trace & TBISECTS)
		free(bisects);

	fftwr_free(noise); fftwr_destroy_plan(noiseplan);

#ifndef DO_PHONEBOOK
	if faster(levels > 0) {
		for (size_t i = 0; i < n; i++)
			free(cinvs[i]);
		free(cinvs);
	}
#endif /* !DO_PHONEBOOK */

#ifdef DO_PHONEBOOK
	free(pbtimes); free(pbvalues);
#else /* !DO_PHONEBOOK */
	free(cinv); free(gamma_); free(g);
#endif /* !DO_PHONEBOOK */
	free(times); free(values);
#ifdef DO_MAX
	free(queue);
#endif /* DO_MAX */

	fftwr_free(eigen);

	gsl_rng_free(rng);
	if (reseed)
		gsl_rng_free(seedrng);

	fftwr_cleanup();
	return EXIT_SUCCESS;
}
