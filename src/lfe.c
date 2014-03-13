/*
 Author: Simen Gaure
 Copyright: 2011, Simen Gaure
 Licence: Artistic 2.0
*/

#include "config.h"

#if defined(_WIN32) || defined(_WIN64) || defined(WIN64)
#define WIN
#endif


#ifdef WIN
#include <windows.h>
#else
#ifndef NOTHREADS
#include <semaphore.h>
#include <pthread.h>
#endif
#endif

#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include <string.h>
/* Need sprintf */
#include <stdio.h>  
#include <R.h>
#include <Rversion.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <R_ext/BLAS.h>

#if defined(R_VERSION) && R_VERSION >= R_Version(3, 0, 0)
typedef R_xlen_t mybigint_t;
#else
typedef int mybigint_t;
#endif
/* If the the number of G() terms times the number of observations
   exceeds the 2^31 4-byte integer limit, define mysize_t as size_t
   This will increase the memory usage, so we wait until it's needed.
*/

#ifdef HUGE_INT 
typedef R_xlen_t mysize_t;
#else
typedef unsigned int mysize_t;
#endif

#ifdef NOTHREADS
#define LOCK_T int*
#define LOCK(l)
#define UNLOCK(l)
#else

#ifdef WIN
#define LOCK_T HANDLE*
#define LOCK(l) WaitForSingleObject(l,INFINITE)
#define UNLOCK(l) ReleaseMutex(l)
#else
#define LOCK_T pthread_mutex_t*
#define LOCK(l) (void)pthread_mutex_lock(l)
#define UNLOCK(l) (void)pthread_mutex_unlock(l)
#endif
#endif


/* Trickery to check for interrupts when using threads */
static void chkIntFn(void *dummy) {
  R_CheckUserInterrupt();
}

/* this will call the above in a top-level context so it won't 
   longjmp-out of context */

int checkInterrupt() {
  return (R_ToplevelExec(chkIntFn, NULL) == 0);
}


/* More trickery, we can't printf in threads since the R API is not
   thread-safe.  So set up a message stack.  This is pushed in the threads
   and popped in the main thread. */

#define MSGLIM 256
static char *msgstack[MSGLIM];
static int msgptr;
/* Craft our own strdup, it's not supported everywhere */
static char *mystrdup(char *s) {
  char *sc = (char*)malloc(strlen(s)+1);
  if(sc != NULL) strcpy(sc,s);
  return(sc);
}
static void pushmsg(char *s, LOCK_T lock) {
#ifdef NOTHREADS  
  REprintf(s);
#else
  LOCK(lock);
  if(msgptr < MSGLIM) {
    msgstack[msgptr++] = mystrdup(s);
  }
  UNLOCK(lock);
#endif
}

static void printmsg(LOCK_T lock) {
#ifdef NOTHREADS
  return;
#else
  char *s;
  int i;
  LOCK(lock);
  for(i = 0; i < msgptr; i++) {
    s = msgstack[i];
    if(s != NULL) {
      REprintf(s);
      free(s);
    }
  }
  msgptr = 0;
  UNLOCK(lock);

#endif
}

typedef struct {
  int nlevels;
  /* group[i] is the level of observation i */
  int *group;
  /* invgpsize[j] is the 1/(size of level j) */
  double *invgpsize;  
  double *gpsize;
  int *gpl; /* group list */
  int *ii;  /* indices into gpl */
  double *x; /* optional interaction covariate */
} FACTOR;


/* 
   I've seen a couple of datasets which don't converge.
   Perhaps we could try a similar trick as in the Kaczmarz-iteration?
   In each factor, every level may be centred independently, we may
   consider each level of each factor a set of projections, and shuffle
   them randomly.  We would need some rewriting for this.  Would it help?
 */


/*
  Centre on all factors in succession.  Vector v. In place.
 */

static R_INLINE void centre(double *v, int N, 
			    FACTOR *factors[], int e, double *means) {

  for(int i=0; i < e; i++) {
    FACTOR *f = factors[i];
    const int *gp = f->group;
    int j=0;

    /* compute means */
    memset(means,0,sizeof(double)* f->nlevels);
    if(NULL != f->x) {
      for(j = 0; j < N; j++) {
	means[gp[j]-1] += v[j]*f->x[j];
      }
    } else {
      for(j = 0; j < N; j++) {
	means[gp[j]-1] += v[j];
      }
    }
    for(j = 0; j < f->nlevels; j++) {
      means[j] *= f->invgpsize[j];
    }

    if(NULL != f->x) {
      for(j = 0; j < N; j++) {
	v[j] -= means[gp[j]-1]*f->x[j];
      }
    } else {
      for(j = 0; j < N; j++) {
	v[j] -= means[gp[j]-1];
      }
    }
  }
}


/*
  Method of alternating projections.  Input v, output res.
*/

static int demean(double *v, int N, double *res,
		  FACTOR *infactors[],int e,double eps, double *means,
		  double *prev, double *prev2,
		  int *stop,int vecnum, LOCK_T lock, const int gkacc) {
  double delta;
  double norm2 = 0.0;
  int i;
  double prevdelta;
  double c;
  double neps;
  int iter;
  int okconv = 1;
  FACTOR *factors[2*e];
  time_t last, now;
  int lastiter;

  for(int i = 0; i < e; i++) factors[i] = infactors[i];
  /* make a copy to result vector, centre() is in-place */
  memcpy(res,v,N*sizeof(double));
  centre(res, N, factors, e, means);

  if(e <= 1) return(1);

  // Initialize
  memcpy(prev2, res, N*sizeof(double));

  norm2 = 1.0;
  for(i = 0; i < N; i++) norm2 += res[i]*res[i];
  norm2 = sqrt(norm2);
  delta = prevdelta = 2*norm2;
  neps = norm2*eps;
  iter = 0;
  last = time(NULL);
  lastiter = 0;

  double target=0;
  if(gkacc) target = 1e-4*neps;

  do {
    iter++;
    if(gkacc) memcpy(prev, res, N*sizeof(double));
    centre(res,N,factors,e,means);

    if(gkacc) {
      // Try the Gearhart-Koshy acceleration, AS 3.2
      // The result vector is an affine combination of the previous x and the centred vector Px
      // tPx + (1-t) x
      // where t = (x, x-Px)/|x - Px|^2   (t can be > 1).
      // It's quite much faster for some datasets (factor 2-3), perhaps more.
      // Though this seems to be related to lower precision.
      double ip = 0, nm=0, nm2=0;
      for(int i = 0; i < N; i++) {
	nm += prev[i]*prev[i];
	ip += prev[i]*(prev[i]-res[i]);
	nm2 += (prev[i]-res[i])*(prev[i] - res[i]);
      }
      if(nm2 > 1e-18*nm) {
	double t = ip/nm2;

	// By Lemma 3.9, t >= 0.5. If we're below zero we're numerically out of whack
	// Or does it mean we have converged? 
	if(t < 0.49) {
	  char buf[256];
	  if(nm2 > 1e-15*nm) {
	    sprintf(buf,"Demeaning of vec %d failed after %d iterations, t=%.1e, nm2=%.1enm\n",
		    vecnum,iter,t,nm2/nm);
	    pushmsg(buf,lock);
	    okconv = 0;
	  }
	  break;
	}
	
	for(int i = 0; i < N; i++) {
	  res[i] = t*res[i] + (1.0-t)*prev[i];
	}
      }
    }

    // make this a power of two, so we don't have to do integer division
#define IBATCH 128
    // Check convergence rate every now and then
    // For the purpose of computing time to convergence, we assume convergence is linear, 
    // i.e. that the decrease in norm since the previous iteration is a constant factor c.
    // To save some computing time we don't check every iteration
    if((iter & (IBATCH-1)) == 0) {
      // compute delta per iter
      delta = 0.0;
      for(int i = 0; i < N; i++) delta += (prev2[i]-res[i])*(prev2[i]-res[i]);
      memcpy(prev2,res,N*sizeof(double));
      // delta is the square norm improvement since last time
      // we normalize it to be per iteration
      // we divide it by the number of iterations to get an improvement per iteration
      delta = sqrt(delta/IBATCH);
      // then we compute how fast the improvement dimishes. We use this to predict when we're done
      // for e == 2 without gkacc it should diminish, so that c < 1, but for e > 2 it may
      // increase, in particular when we do acceleration. Then it is difficult to predict when we will be done
      c = pow(delta/prevdelta,1.0/IBATCH);
      prevdelta = delta;
    /* c is the convergence rate per iteration, at infinity they add up to 1/(1-c).
       This is true for e==2, but also in the ballpark for e>2 we guess.
       Only fail if it's after some rounds. It seems a bit unstable in the
       beginning.  This is to be expected due to Deutsch and Hundal's result 
       If we're doing acceleration, we test above for t < 0.5 for failure
    */

      if(!gkacc && c >= 1.0 && iter > 100) {
	char buf[256];
	sprintf(buf,"Demeaning of vec %d failed after %d iterations, c-1=%.0e, delta=%.1e\n",
		vecnum,iter,c-1.0,delta);
	pushmsg(buf,lock);
	okconv = 0;
	break;
      }

      if(c < 1.0) target = (1.0-c)*neps;
      // target should not be smaller than 3 bits of precision in each coordinate
      if(target < N*1e-15) target = N*1e-15;
      now = time(NULL);
      if(now - last >= 3600 && delta > 0.0) {
	int reqiter;
	double eta;
	char tu;
	char buf[256];
	reqiter = log(target/delta)/log(c);
	eta = 1.0*(now-last)*reqiter/(iter-lastiter);
	if(eta < 0) eta = NA_REAL; 
	tu = 's';
	if(eta > 3600.0) {
	  eta /= 3600.0;
	  tu='h';
	}
	if(gkacc&&0) {
	  sprintf(buf,"...centering vec %d iter %d, delta=%.1e(target %.1e)\n",
		  vecnum,iter,delta,target);

	} else {
	  sprintf(buf,"...centering vec %d iter %d, c=1-%.1e, delta=%.1e(target %.1e), ETA in %.1f%c\n",
		  vecnum,iter,1.0-c,delta,target,eta,tu);
	}
	pushmsg(buf,lock);
	lastiter = iter;
	last = now;
      }
    }

#ifdef NOTHREADS
    R_CheckUserInterrupt();
#else
    if(*stop != 0) {okconv = 0; break;}
#endif

  } while(delta > target);
  return(okconv);
}


/*
  Argument block for threads.
 */
typedef struct {
  int N;
  int K;
  double **v;
  double **res;
  FACTOR **factors;
  int e;
  double eps;
  int nowdoing;
  int badconv;
  time_t last;
  time_t start;
  int done;
  int quiet;
  int stop;
  int threadnum;
  double **means;
  double **tmp1, **tmp2;
  LOCK_T lock;
  int gkacc;
#ifndef NOTHREADS
#ifdef WIN
  /*  HANDLE *lock;*/
#else
  /*  pthread_mutex_t *lock;*/
  int running;
#ifdef HAVE_SEM
  sem_t finished;
#endif
#endif
#endif
} PTARG;


/*
  Thread routine for demeaning.
  Threads pick a vector to demean by looking
  in arg->nowdoing, a common counter.  Proteced
  by a mutex.
 */




/* The thread routine 
   Each thread centres an entire vector.
*/
#ifdef WIN
DWORD WINAPI demeanlist_thr(LPVOID varg) {
#else
static void *demeanlist_thr(void *varg) {
#endif
  PTARG *arg = (PTARG*) varg;
  double *means;
  double *tmp1, *tmp2;
  int vecnum;
  int okconv = 0;
  int myid;
  /* Find the preallocated means buffer */
  LOCK(arg->lock);
  myid = arg->threadnum++;
  UNLOCK(arg->lock);
  means = arg->means[myid];
  tmp1 = arg->tmp1[myid];
  tmp2 = arg->tmp2[myid];
  while(1) {
    time_t now;
    /* Find next vector */
    LOCK(arg->lock);
    vecnum = (arg->nowdoing)++;
    UNLOCK(arg->lock);
    if(vecnum >= arg->K) break;

    okconv = demean(arg->v[vecnum],arg->N,arg->res[vecnum],
		    arg->factors,arg->e,arg->eps,means,
		    tmp1,tmp2,
		    &arg->stop,vecnum+1,
		    arg->lock,arg->gkacc);

    LOCK(arg->lock);
    now = time(NULL);
    if(!okconv) {
      arg->badconv++;
    }
    (arg->done)++;
    if(arg->quiet > 0 && now > arg->last + arg->quiet && arg->K > 1) {
      char buf[256];
      sprintf(buf,"...finished centering vector %d of %d in %d seconds\n",
	     arg->done,arg->K,(int)(now-arg->start));
      UNLOCK(arg->lock); /* release lock, pushmsg takes it */
      pushmsg(buf,arg->lock);
      LOCK(arg->lock);
      arg->last = now;
    }
    UNLOCK(arg->lock);
  }
  /* signal we've finished */
#ifndef NOTHREADS
#ifndef WIN
  LOCK(arg->lock);
  /* yuck, main thread can't wait for zero on a semaphore, so do it this way.
   Last leaver turns out the light.*/
  arg->running--;
#ifdef HAVE_SEM
  if(arg->running == 0) sem_post(&arg->finished);
#endif
  UNLOCK(arg->lock);
#endif
#endif
  return 0;
}


/* main C entry-point for demeaning.  Called from the R-entrypoint */
static int demeanlist(double **vp, int N, int K, double **res,
		      FACTOR *factors[], int e, double eps,int cores,
		      int quiet, const int gkacc) {
  PTARG arg;
  int numthr = 1;
  int thr;
  int maxlev;
  int i;
#ifndef NOTHREADS
#ifdef WIN
  HANDLE *threads;
  DWORD *threadids;
  HANDLE lock;
#else
  pthread_t *threads;
  pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
#endif

#ifndef NOTHREADS
  /* Set up the message stack */
  msgptr = 0;
#endif

  numthr = cores;
  if(numthr > K) numthr = K;
  if(numthr < 1) numthr = 1;

#ifdef WIN
  lock = CreateMutex(NULL,FALSE,NULL);
  threads = (HANDLE*) R_alloc(numthr,sizeof(HANDLE));
  threadids = (DWORD*) R_alloc(numthr,sizeof(DWORD));
#else
  threads = (pthread_t*) R_alloc(numthr,sizeof(pthread_t));
#ifdef HAVE_SEM
  if(sem_init(&arg.finished,0,0) != 0) error("sem_init failed, errno=%d",errno);
#endif
  arg.running = numthr;
#endif
  arg.lock = &lock;
#endif

  maxlev = 0;
  for(i = 0; i < e; i++)
    if(maxlev < factors[i]->nlevels) maxlev = factors[i]->nlevels;

  arg.means = (double **) R_alloc(numthr,sizeof(double*));
  arg.tmp1 = (double **) R_alloc(numthr,sizeof(double*));
  arg.tmp2 = (double **) R_alloc(numthr,sizeof(double*));
  for(thr = 0; thr < numthr; thr++) {
    arg.means[thr] = (double*) R_alloc(maxlev,sizeof(double));
    if(e > 1) {
      arg.tmp1[thr] = (double*) R_alloc(N,sizeof(double));
      arg.tmp2[thr] = (double*) R_alloc(N,sizeof(double));
    }
  }

  arg.badconv = 0;
  arg.N = N;
  arg.v = vp;
  arg.res = res;
  arg.K = K;
  arg.nowdoing = 0;
  arg.factors = factors;
  arg.last = time(NULL);
  arg.start = arg.last;
  arg.done = 0;
  arg.e = e;
  arg.eps = eps;
  arg.quiet = quiet;
  arg.stop = 0;
  arg.threadnum = 0;
  arg.gkacc = gkacc;

#ifdef NOTHREADS
  demeanlist_thr((void*)&arg);
#else
  /* Do it in separate threads */
  for(thr = 0; thr < numthr; thr++) {
#ifdef WIN
    threads[thr] = CreateThread(NULL,0,demeanlist_thr,&arg,0,&threadids[thr]);
    if(0 == threads[thr]) error("Failed to create thread");
#else
    int stat = pthread_create(&threads[thr],NULL, demeanlist_thr,&arg);
    if(0 != stat) error("Failed to create thread, stat=%d",stat);
#endif
  }

  /* wait for completion */
  /* We want to check for interrupts regularly, and
     set a stop flag */
  while(1) {
    printmsg(arg.lock);
    if(arg.stop == 0 && checkInterrupt()) {
      REprintf("...stopping centering threads...\n");
      arg.stop=1;
    }

#ifdef WIN
    if(WaitForMultipleObjects(numthr,threads,TRUE,3000) != WAIT_TIMEOUT) {
      for(thr = 0; thr < numthr; thr++) {
	CloseHandle(threads[thr]);
      }
      CloseHandle(lock);
      break;
    }
#else
    {
#ifndef HAVE_SEM
      struct timespec atmo = {0,50000000};
      /* Kludge when no timedwait, i.e. MacOSX */
      
      if(arg.stop == 0) nanosleep(&atmo,NULL);
      if(arg.stop == 1 || arg.running == 0) {
#else
	struct timespec tmo = {time(NULL)+3,0};
	if(arg.stop == 1 || sem_timedwait(&arg.finished,&tmo) == 0) {
#endif
	  for(thr = 0; thr < numthr; thr++) {
	    (void)pthread_join(threads[thr], NULL);
	  }
#ifdef HAVE_SEM
	  sem_destroy(&arg.finished);
#endif
	  pthread_mutex_destroy(arg.lock);
	  break;
	}
      }
#endif
    }
#endif

    /* print remaining messages */
  printmsg(arg.lock);
  if(arg.stop == 1) error("centering interrupted");
  if(quiet > 0 && arg.start != arg.last && K > 1)
    REprintf("...%d vectors centred in %d seconds\n",K,(int)(time(NULL)-arg.start));

  return(arg.badconv);
}



/* Helper: Make lists of members of each level */

static void invertfactor(FACTOR *f,int N) {
  int nlev = f->nlevels;
  int *curoff;
  int i;
  f->ii = (int*) R_alloc(nlev+1,sizeof(int));
  f->gpl = (int*) R_alloc(N,sizeof(int));

  memset(f->ii,0,sizeof(int)*(nlev+1));

  /* find sizes of groups */
  for(i = 0; i < N; i++) {
    f->ii[f->group[i]]++;
  }

  /* cumulative */
  for(i = 1; i < nlev+1; i++) {
    f->ii[i] += f->ii[i-1];
  }

  curoff = Calloc(nlev+1,int);
  for(i = 0; i < N; i++) {
    int gp = f->group[i]-1;
    f->gpl[f->ii[gp]+curoff[gp]] = i;
    curoff[gp]++;
  }
  Free(curoff);
}

static FACTOR** makefactors(SEXP flist) {
  FACTOR **factors;
  int numfac = LENGTH(flist);
  int N=0;

  numfac = 0;
  for(int i = 0; i < LENGTH(flist); i++) {
    SEXP sf = VECTOR_ELT(flist,i);
    SEXP xattr = getAttrib(sf,install("x"));    
    if(isNull(xattr)) {
      numfac++;
      continue;
    }
    if(!isMatrix(xattr)) {
      numfac++;
      continue;
    }
    numfac += ncols(xattr);
  }

  factors = (FACTOR**) R_alloc(numfac+1,sizeof(FACTOR*));
  factors[numfac] = NULL;
  int truefac = 0;
  for(int i = 0; i < LENGTH(flist); i++) {
    int len;
    FACTOR *f;
    len = LENGTH(VECTOR_ELT(flist,i));
    if(i == 0) {
      N = len;
    } else if(len != N) {
      error("All factors must have the same length %d %d",len,N);
    }
    
    f = factors[truefac++] = (FACTOR*) R_alloc(1,sizeof(FACTOR));
    f->group = INTEGER(VECTOR_ELT(flist,i));
    f->nlevels = LENGTH(getAttrib(VECTOR_ELT(flist,i),R_LevelsSymbol));
    SEXP xattr = getAttrib(VECTOR_ELT(flist,i),install("x"));
    if(isNull(xattr)) {
      f->x = NULL;
    } else {
      if(isMatrix(xattr)) {
	if(nrows(xattr) != len) {
	  error("Factor interaction terms (%d) must have the same length (%d) as the factor",
		LENGTH(xattr),len);
	}
	truefac--;
	for(int j = 0; j < ncols(xattr); j++) {
	  FACTOR *g = factors[truefac++] = (FACTOR*) R_alloc(1,sizeof(FACTOR));
	  g->group = f->group;
	  g->nlevels = f->nlevels;
	  g->x = &REAL(xattr)[(mybigint_t)nrows(xattr)*j];
	}
      } else {
	if(LENGTH(xattr) != len) {
	  error("Factor interaction terms (%d) must have the same length (%d) as the factor",
		LENGTH(xattr),len);
	}
	f->x = REAL(xattr);
      }
    }
  }

  /* Make array for holding precomputed group levels 
     Now, what about entries which don't belong to a group
     I.e. NA entries, how is that handled, I wonder.
     seems to be negative.  Anyway, we fail on them. No, we don't */
  
  for(int i = 0; i < truefac; i++) {
    FACTOR *f = factors[i];
    f->gpsize = (double *)R_alloc(f->nlevels,sizeof(double));
    f->invgpsize = (double *)R_alloc(f->nlevels,sizeof(double));
    memset(f->gpsize,0,f->nlevels*sizeof(double));
    /* first count it */
    for(int j = 0; j < N; j++) {
	/* skip entries without a group, do we need that? */
	/* if(f->group[j] < 1) error("Factors can't have missing levels"); */
      if(f->group[j] > 0) {
	if(NULL == f->x)
	  f->gpsize[f->group[j]-1] += 1.0;
	else
	  f->gpsize[f->group[j]-1] += f->x[j]*f->x[j];
      } else {
	error("Factors can't have missing levels");
      }
      }
      /* then invert it, it's much faster to multiply than to divide */
      /* in the iterations */
    for(int j = 0; j < f->nlevels; j++) {
      f->invgpsize[j] = 1.0/f->gpsize[j];
    }
  }
  return(factors);
}


static double kaczmarz(FACTOR *factors[],int e, mysize_t N, double *R, double *x,
		       double eps, mysize_t *work, int *stop, LOCK_T lock) {
		       //		       double eps, double *newR, int *indices, int *stop, LOCK_T lock) {
  /* The factors define a matrix D, we will solve Dx = R
     There are N rows in D, each row a_i contains e non-zero entries, one
     for each factor, the level at that position.
     We iterate on k, start with x=0, an iteration consists
     of adding a multiple of row i (with i = k %% N), the multiplying coefficient
     is (R[i] - (a_i,x))/e 

     To get better memory locality, we create an (e X N)-matrix with the non-zero indices

     Now, as an afterthought we have made possible to interact a covariate
     with each factor, it's in factors[i]->x
     The theory is exactly the same, but the projection is different.
     I.e. a_i are not 0 or one, and we should not divide by e, but by
     ||a_i||^2.

  */

  double einv = 1.0/e;
  double norm2;
  double prevdiff,neweps;
  double c,sd;
  int iter=0;

  int ie;
  int newN;

  /* We should remove duplicates, at least when they're consecutive.
     Non-consecutives duplicates are just alternating projections and may
     not be too inefficient.  Besides, removing them is more work...
     If any factor is interacted with a covariate, we don't remove
     anything
  */

  int hasinteract = 0;
  for(int i = 0; i < e; i++)
    if(NULL != factors[i]->x) hasinteract = 1;
  mysize_t workpos = 0;

  /* Do the doubles first to keep alignment */
  double *newR = (double *) work;
  mysize_t *indices = (mysize_t *) &work[workpos += N*sizeof(double)/sizeof(mysize_t)];
  mysize_t *perm = (mysize_t *) &work[workpos += e*N];
  mysize_t *prev = (mysize_t *) &work[workpos += N];
  mysize_t *this = (mysize_t *) &work[workpos += e];

  int numlev = 0;
  for(int i = 0; i < e; i++) {
    numlev += factors[i]->nlevels;
    prev[i] = 0;
  }
  newN = 0;
  ie = 0;
  for(mysize_t i = 0; i < N; i++) {
    perm[i] = i;
    for(int j = 0; j < e; j++) {
      this[j] = factors[j]->group[i];
    }
    if(hasinteract || (memcmp(this,prev,e*sizeof(int)) != 0) ) {
      int nlev=0;
      /* not duplicate, store in indices */
      
      for(int j = 0; j < e; j++) {
	indices[ie+j] = this[j]-1 + nlev;
	nlev += factors[j]->nlevels;
      }
      newR[newN] = R[i];
      newN++;
      ie += e;
      memcpy(prev,this,e*sizeof(int));
    }
  }

  /* 
     At this point we should perhaps randomly shuffle the equations.
     We don't know when this ought to be done, but we've only seen it when
     there are many factors.
     We want to use unif_rand to be able to get reproducible results, 
     at least for single-threaded things, and keeping the random number
     generator in the same state when we return to R.
     The unif_rand isn't concurrency-safe, so protect with mutex.
     The Get/PutRNGstate is done in the main thread.
     Knuth-Fisher-Yates shuffle.
     The randomization done here should not have consequences for the result
     of the routine, only the speed.
  */
  LOCK(lock);
  if(e > 1) {
    for(mysize_t i = newN-1; i > 0; i--) {
      mysize_t j;
      /* Pick j between 0 and i inclusive */
      j = (mysize_t) floor((i+1) * unif_rand());
      if(j == i) continue;
      /* exchange newR[i] and newR[j]
	 as well as indices[i*e:i*e+e-1] and indices[j*e:j*e+e-1]
      */
      double dtmp = newR[j];
      newR[j] = newR[i];
      newR[i] = dtmp;
      mysize_t itmp = perm[j];
      perm[j] = perm[i];
      perm[i] = itmp;
      for(mysize_t k = 0; k < e; k++) {
	mysize_t itmp;
	itmp = indices[j*e+k];
	indices[j*e+k] = indices[i*e+k];
	indices[i*e+k] = itmp;
      }
    }
  }
  UNLOCK(lock);

  /* Then, do the Kaczmarz iterations */
  norm2 =0.0;
  for(mysize_t i = 0; i < newN; i++) norm2 += newR[i]*newR[i];
  norm2 = sqrt(norm2);
  prevdiff = 2*norm2;
  neweps = eps*norm2;

  do {
    mysize_t ie = 0; /* equals i*e; integer multiplication is slow, keep track instead */
    double diff = 0.0;
    if(hasinteract) {
      for(mysize_t i = 0; i < newN; i++,ie+=e) {
	mysize_t ip = perm[i];
	double upd = 0.0;
	double ai2 = 0.0;
	upd = newR[i];

	/* Subtract inner product */
	for(int j = 0; j < e; j++) {
	  mysize_t idx = indices[ie + j];
	  double *fx = factors[j]->x;
	  if(NULL == fx) {
	    upd -= x[idx];
	    ai2 += 1.0;
	  } else {
	    upd -= x[idx]*fx[ip];
	    ai2 += fx[ip]*fx[ip];
	  }
	}
	upd /= ai2;
	for(int j = 0; j < e; j++) {
	  mysize_t idx = indices[ie + j];
	  double *fx = factors[j]->x;
	  if(NULL == fx) {
	    x[idx] += upd;
	    diff += upd*upd;
	  } else {
	    x[idx] += upd*fx[ip];
	    diff += upd*upd*fx[ip]*fx[ip];
	  }
	}
      }
      diff *= einv;
    } else {
      for(mysize_t i = 0; i < newN; i++,ie+=e) {
	double upd = 0.0;
	upd = newR[i];
	
	/* Subtract inner product */
	for(int j = 0; j < e; j++) {
	  mysize_t idx = indices[ie + j];
	  upd -= x[idx];
	}
	upd *= einv;
      
	/* Update x */
	for(int j = 0; j < e; j++) {
	  mysize_t idx = indices[ie + j];      
	  x[idx] += upd;
	}
	/* Keep track of update */
	diff += upd*upd;
      }
    }
    iter++;
    /* Relax a bit */
    /*    R_CheckUserInterrupt();*/
    sd = sqrt(diff*e);
    c = sd/prevdiff;
    prevdiff = sd;
    if(c >= 1.0 && iter > 20) {
      char buf[256];
      sprintf(buf,"Kaczmarz failed in iter %d*%d, c=1-%.0e, delta=%.1e, eps=%.1e\n",iter,newN,1.0-c,sd,neweps);
      pushmsg(buf,lock);
      break;
    }
#ifdef NOTHREADS
    R_CheckUserInterrupt();
#else
    if(*stop != 0) return(0);
#endif
  } while(sd >= neweps*(1.0-c) && neweps > 1e-15);

  return(sd);
}

typedef struct {
  int nowdoing;
  
  double **source;
  double **target;
  FACTOR **factors;
  int e;
  int N;
  int threadnum;
  int stop;
  int numvec;
  double eps;
  double *tol;
  //  double **newR;
  //  int **indices;
  mysize_t **work;
  LOCK_T lock;
#ifndef NOTHREADS
#ifdef WIN
#else
  int running;
#ifdef HAVE_SEM
  sem_t finished;
#endif
#endif
#endif
} KARG;

#ifdef WIN
DWORD WINAPI kaczmarz_thr(LPVOID varg) {
#else
static void *kaczmarz_thr(void *varg) {
#endif
  KARG *arg = (KARG*) varg;
  int myid;
  int vecnum;
  /* Get a thread id */
  LOCK(arg->lock);
  myid = arg->threadnum++;
  UNLOCK(arg->lock);
  while(1) {
    LOCK(arg->lock);
    vecnum = arg->nowdoing++;
    UNLOCK(arg->lock);
    if(vecnum >= arg->numvec) break;
    (void) kaczmarz(arg->factors,arg->e,arg->N,
		    arg->source[vecnum],arg->target[vecnum],arg->eps,
		    arg->work[myid], &arg->stop, arg->lock);
  }
#ifndef NOTHREADS
#ifndef WIN
  LOCK(arg->lock);
  arg->running--;
#ifdef HAVE_SEM
  if(arg->running == 0) sem_post(&arg->finished);
#endif
  UNLOCK(arg->lock);
#endif
#endif
  return 0;
}

static SEXP R_kaczmarz(SEXP flist, SEXP vlist, SEXP Reps, SEXP initial, SEXP Rcores) {

  double eps = REAL(Reps)[0];
  /*  double *R = REAL(RR);*/
  FACTOR **factors;
  double *init = 0;
  int cores = INTEGER(Rcores)[0];
  int numfac;
  mysize_t N = 0;
  //  int i;
  int sumlev = 0;
  int listlen;
  int numvec;
  SEXP reslist;
  KARG arg;
  double **vectors, **target;
  int cnt;
  int numthr = 1;
#ifndef NOTHREADS
  int thr;
#ifdef WIN
  HANDLE *threads;
  DWORD *threadids;
  HANDLE lock;
#else
  pthread_t *threads;
  pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
#endif
#endif

#ifndef NOTHREADS
  /* Set up the message stack */
  msgptr = 0;
#endif
  
  if(!isNull(initial)) {
    init = REAL(initial);
  }
  PROTECT(flist = AS_LIST(flist));
  //  numfac = LENGTH(flist);

  factors = makefactors(flist);
  numfac = 0;
  for(FACTOR **f = factors; *f != NULL; f++) numfac++;
  N = LENGTH(VECTOR_ELT(flist,0));
  for(int i = 0; i < numfac; i++)
    sumlev += factors[i]->nlevels;

  if(!isNull(initial) && LENGTH(initial) != sumlev)
    error("Initial vector must have length %d, but is %d\n",sumlev, LENGTH(initial));

  /* Then the vectors */

  PROTECT(vlist = AS_LIST(vlist));
  listlen = LENGTH(vlist);
  PROTECT(reslist = NEW_LIST(listlen));

  /* First, count the number of vectors in total */
  numvec = 0;
  for(int i = 0; i < listlen; i++) {
    SEXP elt = VECTOR_ELT(vlist,i);
    if(!IS_NUMERIC(elt)) error("Entries must be numeric");
    /* Each entry in the list is either a vector or a matrix */
    if(!isMatrix(elt)) {
      if(LENGTH(elt) != N) 
	error("Vector length (%d) must be equal to factor length (%d)",LENGTH(elt),N);
      numvec++;
    } else {
      if(nrows(elt) != N)
	error("Vector length must be equal to factor length %d %d ",nrows(elt),N);
      numvec += ncols(elt);
    }
  }

  /* Allocate pointers to source vectors */
  vectors = (double **)R_alloc(numvec,sizeof(double*));
  /* Allocate pointers to result vectors */
  target = (double**) R_alloc(numvec,sizeof(double*));
  /* Loop through list again to set up result structure */
  cnt = 0;
  for(int i = 0; i < listlen; i++) {
    SEXP elt = VECTOR_ELT(vlist,i);
    if(!isMatrix(elt)) {
      /* It's a vector */
      SEXP resvec;
      vectors[cnt] = REAL(elt);
      PROTECT(resvec = allocVector(REALSXP,sumlev));
      target[cnt] = REAL(resvec);
      SET_VECTOR_ELT(reslist,i,resvec);
      UNPROTECT(1);
      cnt++;
    } else {
      /* It's a matrix */
      int cols = ncols(elt);
      SEXP mtx;
      /* Allocate a matrix */
      PROTECT(mtx = allocMatrix(REALSXP,sumlev,cols));
      SET_VECTOR_ELT(reslist,i,mtx);
      UNPROTECT(1);
      /* Set up pointers */
      for(int j = 0; j < cols; j++) {
	vectors[cnt] = REAL(elt) + j*N;
	target[cnt] = REAL(mtx) + j*sumlev;
	cnt++;
      }
    }
  }


  for(int cnt = 0; cnt < numvec; cnt++) {
    if(init != 0)
      for(int i = 0; i < sumlev; i++) target[cnt][i] = init[i];
    else
      for(int i = 0; i < sumlev; i++) target[cnt][i] = 0.0;
  }

  
  /* set up for threading */

  numthr = cores;
  if(numthr > numvec) numthr = numvec;
  if(numthr < 1) numthr = 1;
  GetRNGstate();
#ifndef NOTHREADS
#ifdef WIN
  lock = CreateMutex(NULL,FALSE,NULL);
  threads = (HANDLE*) R_alloc(numthr,sizeof(HANDLE));
  threadids = (DWORD*) R_alloc(numthr,sizeof(DWORD));
#else
  threads = (pthread_t*) R_alloc(numthr,sizeof(pthread_t));
  arg.running = numthr;
#ifdef HAVE_SEM
  if(sem_init(&arg.finished,0,0) != 0) error("sem_init failed, errno=%d",errno);
#endif
#endif
  arg.lock = &lock;
#endif

  arg.factors = factors;
  arg.source = vectors;
  arg.target = target;
  arg.nowdoing = 0;
  arg.threadnum = 0;
  arg.e = numfac;
  arg.eps = eps;
  arg.numvec = numvec;
  arg.N = N;
  arg.stop = 0;
  arg.work = (mysize_t**) R_alloc(numthr, sizeof(mysize_t*));

  // when allocating the work, we use mysize_t, but parts of it is accessed as double which may be
  // larger. So we allocate some more (8*sizeof(mysize_t) more), and adjust the address so it's aligned on a double
  // When using it, make sure we do all the doubles first.

#ifdef NOTHREADS
  arg.work[0] = (mysize_t*) R_alloc(numfac*N + N*sizeof(double)/sizeof(mysize_t) + N + 2*numfac+8, sizeof(mysize_t));
  uintptr_t amiss = (uintptr_t) arg.work[0] % sizeof(double);
  if(amiss != 0) arg.work[0] = (mysize_t*) ((uintptr_t)arg.work[0] + sizeof(double)-amiss);
  kaczmarz_thr((void*)&arg);
#else
  /* Do it in separate threads */
  for(thr = 0; thr < numthr; thr++) {
    // allocate some thread-specific storage, we can't use R_alloc in a thread
    arg.work[thr] = (mysize_t*) R_alloc(numfac*N + N*sizeof(double)/sizeof(mysize_t) + N + 2*numfac+8, sizeof(mysize_t));
    uintptr_t amiss = (uintptr_t) arg.work[thr] % sizeof(double);
    if(amiss != 0) arg.work[thr] = (mysize_t*) ((uintptr_t)arg.work[thr] + sizeof(double)-amiss);

#ifdef WIN
    threads[thr] = CreateThread(NULL,0,kaczmarz_thr,&arg,0,&threadids[thr]);
    if(0 == threads[thr]) error("Failed to create kaczmarz thread");
#else
    int stat = pthread_create(&threads[thr],NULL,kaczmarz_thr,&arg);
    if(0 != stat) error("Failed to create kaczmarz thread, stat=%d",stat);
#endif
  }
  /* wait for completion */
  /* We want to check for interrupts regularly, and
     set a stop flag */
  while(1) {
    printmsg(arg.lock);
    if(arg.stop == 0 && checkInterrupt()) {
      REprintf("...stopping Kaczmarz threads...\n");
      arg.stop=1;
    }

#ifdef WIN
    if(WaitForMultipleObjects(numthr,threads,TRUE,3000) != WAIT_TIMEOUT) {
      for(thr = 0; thr < numthr; thr++) {
	CloseHandle(threads[thr]);
      }
      CloseHandle(lock);
      break;
    }
#else
    {
#ifndef HAVE_SEM
      struct timespec atmo = {0,50000000};
      /* Kludge in MacOSX because no timedwait */
      
      if(arg.stop == 0) nanosleep(&atmo,NULL);
      if(arg.stop == 1 || arg.running == 0) {
#else
      struct timespec tmo = {time(NULL)+3,0};
      if(arg.stop == 1 || sem_timedwait(&arg.finished,&tmo) == 0) {
#endif
	for(thr = 0; thr < numthr; thr++) {
	  (void)pthread_join(threads[thr], NULL);
	}
#ifdef HAVE_SEM
	sem_destroy(&arg.finished);
#endif
	pthread_mutex_destroy(arg.lock);
	break;
      }
    }
#endif
  }
#endif
  PutRNGstate();
  /* Print any remaining messages */
  printmsg(arg.lock);
  if(arg.stop == 1) error("Kaczmarz interrupted");
  UNPROTECT(3);
  return(reslist);
}

/*
 Now, for the R-interface.  We only export a demeanlist function
 It should be called from R with 
 demeanlist(matrix,list(f1,f2,...),double eps)
 and preferably update matrix inplace to save copying memory
 unless that's completely un-R'ish. Hmm, it is.
 We use the .C-interface (do we?), or do we use .Call?
 We should take the factor list and store it in our
 format (unless it's already suitable... let's see)
 and call our demeanlist and return
*/

static SEXP R_demeanlist(SEXP vlist, SEXP flist, SEXP Ricpt, SEXP Reps,
			 SEXP scores, SEXP quiet, SEXP gkacc, SEXP Rmeans) {
  int numvec;
  int numfac;
  int cnt;
  double **vectors;
  double **target;
  double eps;
  SEXP reslist;
  FACTOR **factors;
  int cores;
  int i;
  int N=0;
  int icpt;
  SEXP badconv;
  int listlen;
  int domeans=0;

  domeans = LOGICAL(Rmeans)[0];
  icpt = INTEGER(Ricpt)[0] - 1; /* convert from 1-based to zero-based */
  eps = REAL(Reps)[0];
  cores = INTEGER(scores)[0];


  PROTECT(flist = AS_LIST(flist));
  //  numfac = LENGTH(flist);
  factors = makefactors(flist);
  numfac = 0;
  for(FACTOR **f = factors; *f != NULL; f++) numfac++;

  N = LENGTH(VECTOR_ELT(flist,0));
  /* Then the vectors */

  PROTECT(vlist = AS_LIST(vlist));
  listlen = LENGTH(vlist);
  PROTECT(reslist = NEW_LIST(listlen));

  /* First, count the number of vectors in total */
  numvec = 0;
  for(i = 0; i < listlen; i++) {
    SEXP elt = VECTOR_ELT(vlist,i);
    if(!IS_NUMERIC(elt)) error("Entries must be numeric");
    /* Each entry in the list is either a vector or a matrix */
    if(!isMatrix(elt)) {
      if(LENGTH(elt) != N) 
	error("Factors and vectors must have the same length %d != %d",N,LENGTH(elt));
      numvec++;
    } else {
      if(nrows(elt) != N)
	error("Factors and vectors must have the same length %d != %d",N,nrows(elt));
      numvec += ncols(elt);
      if(icpt >= 0 && icpt < ncols(elt)) numvec--;
    }
  }

  /* Allocate pointers to source vectors */
  vectors = (double **)R_alloc(numvec,sizeof(double*));
  /* Allocate pointers to result vectors */
  target = (double**) R_alloc(numvec,sizeof(double*));
  /* Loop through list again to set up result structure */
  cnt = 0;
  for(i = 0; i < listlen; i++) {
    SEXP elt = VECTOR_ELT(vlist,i);
    if(!isMatrix(elt)) {
      /* It's a vector */
      SEXP resvec;
      vectors[cnt] = REAL(elt);
      PROTECT(resvec = allocVector(REALSXP,LENGTH(elt)));
      DUPLICATE_ATTRIB(resvec, elt);
      target[cnt] = REAL(resvec);
      SET_VECTOR_ELT(reslist,i,resvec);
      UNPROTECT(1);
      cnt++;
    } else {
      /* It's a matrix */
      int cols = ncols(elt);
      int rcols = cols;
      int j;
      SEXP mtx;
      /* Allocate a matrix */
      if(icpt >= 0 && icpt < cols) rcols--;
      PROTECT(mtx = allocMatrix(REALSXP,N,rcols));
      SET_VECTOR_ELT(reslist,i,mtx);
      UNPROTECT(1);
      /* Set up pointers */
      rcols = 0;
      for(j = 0; j < cols; j++) {
	if(j == icpt) continue;
	vectors[cnt] = REAL(elt) + j * (mybigint_t)N;
	target[cnt] = REAL(mtx) + (rcols++) * (mybigint_t)N;
	cnt++;
      }
    }
  }

  /* Then do stuff */
  PROTECT(badconv = allocVector(INTSXP,1));
  INTEGER(badconv)[0] = demeanlist(vectors,N,numvec,target,factors,numfac,
				   eps,cores,INTEGER(quiet)[0],
				   INTEGER(gkacc)[0]);

  if(INTEGER(badconv)[0] > 0)
    warning("%d vectors failed to centre to tolerance %.1le",INTEGER(badconv)[0],eps);
  setAttrib(reslist,install("badconv"),badconv);
  SET_NAMES(reslist, GET_NAMES(vlist));

  if(domeans) {
    for(int i = 0; i < numvec; i++) {
      double *srcvec = vectors[i];
      double *dstvec = target[i];
      for(int j = 0; j < N; j++) {
	dstvec[j] = (srcvec[j]-dstvec[j]);
      }
    }
  }
  /* unprotect the reslist */
  UNPROTECT(4);
  return(reslist);
}

// an .External version of demeanlist
static SEXP RE_demeanlist(SEXP args) {
  args = CDR(args); /* 'name' */
  /* Now, loop through the args and set up arguments for the 
     .Call demeanlist */
  SEXP flist=R_NilValue, Ricpt=R_NilValue, Reps=R_NilValue,
    scores=R_NilValue, quiet=R_NilValue, gkacc=R_NilValue, Rmeans=R_NilValue;
  // The number of vectors/matrices to be centered
  if(length(args) < 7) error("Passed %d arguments to RE_demeanlist", length(args));
  int numvec = length(args) - 7;
  int curvec = 0;
  SEXP vlist,lnames;
  
  PROTECT(vlist = allocVector(VECSXP,numvec));
  PROTECT(lnames = allocVector(STRSXP, numvec));
  for(int i = 0; args != R_NilValue; i++, args = CDR(args)) {
    SEXP tag = PRINTNAME(TAG(args));
    const char *name =
      isNull(TAG(args)) ? "" : CHAR(tag);
    SEXP el = CAR(args);
    if(strcmp(name,"fl") == 0) {
      flist = el;
    } else if(strcmp(name,"icpt") == 0) {
      Ricpt = el;
    } else if(strcmp(name,"eps") == 0) {
      Reps = el;
    } else if(strcmp(name,"threads") == 0) {
      scores = el;
    } else if(strcmp(name,"progress") == 0) {
      quiet = el;
    } else if(strcmp(name,"accel") == 0) {
      gkacc = el;
    } else if(strcmp(name,"means") == 0) {
      Rmeans = el;
    } else {
      //      Rprintf("set vec %d to %s, length=%d, tag-len %d\n",curvec,name,length(el), length(tag));
      SET_VECTOR_ELT(vlist,curvec,el);
      SET_STRING_ELT(lnames,curvec,tag);
      curvec++;
    }
  }

  SET_NAMES(vlist, lnames);
  SEXP ret = R_demeanlist(vlist, flist, Ricpt, Reps, scores, quiet, gkacc, Rmeans);
  UNPROTECT(2);
  return ret;
}




static SEXP R_scalecols(SEXP mat, SEXP vec) {
  if(!isMatrix(mat)) error("first argument should be a matrix");
  mybigint_t col = ncols(mat), row = nrows(mat);
  if(row != LENGTH(vec)) error("length of vector %d is different from number of rows %d",LENGTH(vec),row);
  double *cmat = REAL(mat);
  double *cvec = REAL(AS_NUMERIC(vec));
  for(mybigint_t j = 0; j < col; j++) {
    double *cc = &cmat[j*row];
    for(mybigint_t i = 0; i < row; i++)
      cc[i] *= cvec[i];
  }
  return R_NilValue;
}


/* 
Then for finding connection components 
From R we take a list of factors, we return
a factor of the same length with the connection
components
*/


/*
  This one is a bit tricky.  We could do it quite elegantly recursively,
  but we would suffer from a deep call-stack.  Hence, we make our own
  stack and push/pop in a loop
 */
static int Components(int **vertices, FACTOR **factors, int K) {
  int stacklevel = 0;
  int *stack;
  int curfac, curvert, curcomp,candvert;
  int startfac=0,startvert=0;
  int ii=0,i;
  int numvert=0;

  /* How big a stack ? */
  /* The number of vertices */

  for(i = 0; i < K; i++) numvert += factors[i]->nlevels;
  /* Never used in threads, so use R's Calloc */
  stack = Calloc(numvert*4,int);
#define PUSH(x) stack[stacklevel++] = x
#define POP(x) x = stack[--stacklevel]
#define PUSHALL {PUSH(startvert); PUSH(startfac); PUSH(curfac); PUSH(ii);}
#define POPALL {POP(ii); POP(curfac); POP(startfac); POP(startvert);}

  curcomp = 1;
  candvert = 0;
  do {
    curvert = candvert;
    curfac = 0;
    /* Find the entire component */

    while(1) {
      /* At the top here, curfac,curvert is a candidate for marking off as 
	 a vertex in curcomp
	 For each iteration:  

	 If it's not marked, mark it, push on stack, go to first datapoint
	 If it's already marked, go to next datapoint for the vertex (incidence matrix)
 	 If data-points are exhausted, go to next factor, start over with datapoints.
	 If factors are exhausted, pop the stack
	 If final stack-frame, we're done with component.
      */
      

      if(vertices[curfac][curvert] == 0) {
	/* Mark new vertex, find incidence list */

	vertices[curfac][curvert] = curcomp;
	PUSHALL;
	startvert = curvert;
	startfac = curfac;
	curfac = (startfac+1)%K;
	ii = factors[startfac]->ii[startvert];
      } else {
	/* Been there, try next in group */
	ii++;
      }
      if(ii >= factors[startfac]->ii[startvert+1]) {
	/* No more, move to next factor */

	curfac = (curfac + 1) % K;
	if(curfac == startfac) {
	  /* This is where we began, pop */
	  /* No more neighbours, go back to previous */

	  POPALL;
	  /* Get outta here */
	  if(0 == stacklevel) break;
	} else {
	  /* start over with data-points */
	  ii = factors[startfac]->ii[startvert];	
	}
      }
      curvert = factors[curfac]->group[factors[startfac]->gpl[ii]]-1;
    }
    /* Find next component */
    while(candvert < factors[0]->nlevels && vertices[0][candvert] != 0) candvert++;
    curcomp++;
  } while(candvert < factors[0]->nlevels);
  Free(stack);
  return(curcomp-1);
#undef PUSH
#undef POP
#undef PUSHALL
#undef POPALL
}
 
// Algorithm from:
// "A Note on the Determination of Connectedness in an N-Way Cross Classification"
// D.L. Weeks and D.R. Williams, Technometrics, vol 6 no 3, August 1964
// There probably exists faster algorithms, this one is quite slow
static void wwcomp(FACTOR *factors[], int numfac, int N, int *newlevels) {
   int level = 0;
   int chead = 0;
   int *newlist = Calloc(N,int);
   int *oldlist = Calloc(N,int);
   int oldstart = 0;
   // For cache-efficiency, make a numfac x N matrix of the factors
   int *facmat = Calloc(numfac*N, int);

   for(mysize_t i = 0; i < N; i++) {
     int *obsp = &facmat[i*numfac];
     newlevels[i] = 0;
     oldlist[i] = i;
     for(int j = 0; j < numfac; j++) {
       obsp[j] = factors[j]->group[i];
     }
   }
   while(oldstart < N) {
     int newidx;
     // set component number for first node in component
     // find the next we haven't checked, it's oldlist[oldstart]
     // increase oldstart by one
     level++;
     chead = oldlist[oldstart++];
     // put it as the first element in our newlist
     newlist[0] = chead;
     newlevels[chead] = level;
     newidx = 1;
     // loop over the list of newly added nodes, including the head
     // Note that we may increase newidx during the loop
     for(int i = 0; i < newidx; i++) {
       mysize_t newnode = newlist[i];
       int *newp = &facmat[newnode*numfac];
       // search for observations with distance 1 from newnode, mark them with level
       for(int jidx = oldstart; jidx < N; jidx++) { 
	 mysize_t trynode = oldlist[jidx];
	 int *tryp = &facmat[trynode*numfac];
	 int dist = 0;
	 // compute distance
	 for(int fi = 0; fi < numfac && dist < 2; fi++) 
	   dist += (newp[fi] != tryp[fi]);
	 //dist += (factors[fi]->group[newnode] != factors[fi]->group[trynode]);
	 // if close, set its level, add it to the list, move the start node
	 // to the empty place in oldlist. 
	 if(dist < 2) {
	   newlevels[trynode] = level;
	   newlist[newidx++] = trynode;
	   oldlist[jidx] = oldlist[oldstart++];
	 }
       }
     }
   }
   Free(facmat);
   Free(newlist);
   Free(oldlist);
 }

/*
R entry-point for conncomp.  Takes a list of factors as input.
 */

 static SEXP R_wwcomp(SEXP flist) {
   int numfac, N;
   FACTOR **factors;
   SEXP result;

  
  numfac = LENGTH(flist);
  if(numfac < 2) error("At least two factors must be specified");

  N = LENGTH(VECTOR_ELT(flist,0));
  for(int i = 0; i < numfac; i++) {
    if(N != LENGTH(VECTOR_ELT(flist,i))) 
      error("Factors must have the same length");
  }

  factors = (FACTOR**) R_alloc(numfac,sizeof(FACTOR*));
  for(int i = 0; i < numfac; i++) {
    FACTOR *f;

    factors[i] = (FACTOR*) R_alloc(1,sizeof(FACTOR));
    f = factors[i];
    f->group = INTEGER(VECTOR_ELT(flist,i));
  }
  PROTECT(result = allocVector(INTSXP,N));
  int *fac = INTEGER(result);
  wwcomp(factors, numfac, N, fac);
  // Now it's time to order the levels by decreasing size, so let's compute the sizes
  int levels = 0;
  for(int i = 0; i < N; i++) if(fac[i] > levels) levels = fac[i];
  double *levsize = (double*) R_alloc(levels, sizeof(double));
  int *index = (int*) R_alloc(levels, sizeof(int));
  for(int i = 0; i < levels; i++) {
    levsize[i] = 0.0;
    index[i] = i;
  }
  for(int i = 0; i < N; i++) levsize[fac[i]-1] = levsize[fac[i]-1]+1;
  revsort(levsize,index,levels);
  int *rindex = (int*) R_alloc(levels, sizeof(int));
  for(int i = 0; i < levels; i++) rindex[index[i]] = i;
  for(int i = 0; i < N; i++) {
    fac[i] = rindex[fac[i]-1]+1;
  }

  
  UNPROTECT(1);
  return result;
}

static SEXP R_conncomp(SEXP flist) {
  int numfac;
  int i;
  
  int N;
  FACTOR **factors;
  int *group;
  int **vertices;
  SEXP result;
  int *resgroup;
  int comps;
  int *levtrl;
  double *gpsiz;
  int *idx;

  
  numfac = LENGTH(flist);
  if(numfac < 2) error("At least two factors must be specified");

  N = LENGTH(VECTOR_ELT(flist,0));
  for(i = 0; i < numfac; i++) {
    if(N != LENGTH(VECTOR_ELT(flist,i))) 
      error("Factors must have the same length");
  }

  factors = (FACTOR**) R_alloc(numfac,sizeof(FACTOR*));
  PROTECT(flist = AS_LIST(flist));
  for(i = 0; i < numfac; i++) {
    FACTOR *f;

    factors[i] = (FACTOR*) R_alloc(1,sizeof(FACTOR));
    f = factors[i];
    f->group = INTEGER(VECTOR_ELT(flist,i));
    f->nlevels = LENGTH(getAttrib(VECTOR_ELT(flist,i),R_LevelsSymbol));
    invertfactor(f,N);
  }

  /* Create vertices */
  vertices = (int**) R_alloc(numfac,sizeof(int*));
  /* Create arrays for them */
  for(i = 0; i < numfac; i++) {
    vertices[i] = (int*) R_alloc(factors[i]->nlevels,sizeof(int));
    /* Assign no component to them*/
    memset(vertices[i],0,sizeof(int)*factors[i]->nlevels);
  }

  /* Do the stuff */
  comps = Components(vertices,factors,numfac);

  /* allocate result structure */
  PROTECT(result = allocVector(INTSXP,N));
  resgroup = INTEGER(result);
  group = factors[0]->group;
  for(i = 0; i < N; i++) {
    resgroup[i] = vertices[0][group[i]-1];
  }
  /* the levels should be ordered by decreasing size. How do we do this? 
     Hmm, we should have a look at revsort supplied by R.
     There must be an easier way, I'm clumsy today.   
  */

  gpsiz = Calloc(comps,double);
  idx = Calloc(comps,int);
  
  for(i = 0; i < comps; i++) idx[i] = i;
  for(i = 0; i < N; i++) {
    gpsiz[resgroup[i]-1]++;
  }

  revsort(gpsiz,idx,comps);
  Free(gpsiz);
  levtrl = Calloc(comps,int);
  for(i = 0; i < comps; i++) levtrl[idx[i]] = i+1;
  Free(idx);

  for(i = 0; i < N; i++) {
    resgroup[i] = levtrl[resgroup[i]-1];
  }
  Free(levtrl);
  
  UNPROTECT(2);
  return(result);
}

// copy-free dimnames<-

static SEXP R_setdimnames(SEXP obj, SEXP nm) {
  setAttrib(obj, R_DimNamesSymbol, nm);
  return(R_NilValue);
}

// copy-free dsyrk
static SEXP R_dsyrk(SEXP inbeta, SEXP inC, SEXP inalpha, SEXP inA) {
  double beta = REAL(AS_NUMERIC(inbeta))[0];
  double alpha = REAL(AS_NUMERIC(inalpha))[0];
  if(!isMatrix(inC)) error("C must be a matrix");
  if(!isMatrix(inA)) error("A must be a matrix");

  if(ncols(inC) != nrows(inC)) {
    error("C must be a square matrix, it is %d x %d",nrows(inC), ncols(inC));
  }
  int N = nrows(inC);
  double *C = REAL(inC);
  if(ncols(inA) != ncols(inC)) {
    error("A (%d x %d) must have the same number of columns as C (%d x %d)",nrows(inA),ncols(inA),nrows(inC),nrows(inC));
  }
  int K = nrows(inA);
  double *A = REAL(inA);
  F77_CALL(dsyrk)("U","T",&N, &K, &alpha, A, &K, &beta, C, &N);
  // fill in the lower triangular part
  for(mybigint_t row=0; row < N; row++) {
    for(mybigint_t col=0; col < row; col++) {
      C[col*N + row] = C[row*N+col];
    }
  }
  return R_NilValue;
}


// perform C <- beta*C + alpha * (A'B + B'A)
 static SEXP R_dsyr2k(SEXP inbeta, SEXP inC, SEXP inalpha, SEXP inA, SEXP inB) {
  double beta = REAL(AS_NUMERIC(inbeta))[0];
  double alpha = REAL(AS_NUMERIC(inalpha))[0];
  if(!isMatrix(inC)) error("C must be a matrix");
  if(!isMatrix(inA)) error("A must be a matrix");
  if(!isMatrix(inB)) error("B must be a matrix");

  if(ncols(inC) != nrows(inC)) {
    error("C must be a square matrix, it is %d x %d",nrows(inC), ncols(inC));
  }
  int N = nrows(inC);
  double *C = REAL(inC);
  if(ncols(inA) != ncols(inC)) {
    error("A (%d x %d) must have the same number of columns as C (%d x %d)",nrows(inA),ncols(inA),nrows(inC),nrows(inC));
  }
  if(nrows(inA) != nrows(inB)) {
    error("A (%d x %d) must have the same number of rows as B (%d x %d)",nrows(inA), ncols(inA), nrows(inB), ncols(inB));
  }
  int K = nrows(inA);
  
  double *A = REAL(inA);
  double *B = REAL(inB);
  F77_CALL(dsyr2k)("U","T",&N, &K, &alpha, A, &K, B, &K, &beta, C, &N);
  // fill in the lower triangular part
  for(mybigint_t row=0; row < N; row++) {
    for(mybigint_t col=0; col < row; col++) {
      C[col*N + row] = C[row*N+col];
    }
  }
  return R_NilValue;
}

static R_CallMethodDef callMethods[] = {
  {"conncomp", (DL_FUNC) &R_conncomp, 1},
  {"wwcomp", (DL_FUNC) &R_wwcomp, 1},
  {"demeanlist", (DL_FUNC) &R_demeanlist, 8},
  {"kaczmarz", (DL_FUNC) &R_kaczmarz, 5},
  {"setdimnames", (DL_FUNC) &R_setdimnames, 2},
  {"scalecols", (DL_FUNC) &R_scalecols, 2},
  {"dsyrk", (DL_FUNC) &R_dsyrk, 4},
  {"dsyr2k", (DL_FUNC) &R_dsyr2k, 5},
  {NULL, NULL, 0}
};
static R_ExternalMethodDef externalMethods[] = {
  {"edemeanlist", (DL_FUNC) &RE_demeanlist, -1},
  {NULL, NULL, 0}
};

void attribute_visible R_init_lfe(DllInfo *info) {
  /* register our routines */
  (void)R_registerRoutines(info,NULL,callMethods,NULL,externalMethods);
}

