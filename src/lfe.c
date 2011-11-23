/*
 Author: Simen Gaure
 Copyright: 2011, Simen Gaure
 Licence: Artistic 2.0
*/

#define CSTACK_DEFNS

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
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#ifndef NOTHREADS
#if !defined(WIN) && !defined(HAVE_UINTPTR_T) && !defined(uintptr_t)
typedef unsigned long uintptr_t;
#endif
extern uintptr_t R_CStackLimit; /* C stack limit */
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


typedef struct {
  int nlevels;
  /* group[i] is the level of observation i */
  int *group;
  /* invgpsize[j] is the 1/(size of level j) */
  double *invgpsize;  
  int *gpl; /* group list */
  int *ii;  /* indices into gpl */
} FACTOR;

#if __GNUC__ > 3 && __GNUC_MINOR__ > 3
/*#pragma GCC optimize ("O3", "loop-block", "unroll-loops")*/
/*#pragma GCC optimize ("O3", "unroll-loops")*/
#endif

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

static R_INLINE double centre(double *v, int N, 
		   FACTOR *factors[], int e, double *means) {

  int i;
  double sum = 0.0;
  
  for(i=0; i < e; i++) {
    FACTOR *f = factors[i];
    const int *gp = f->group;
    int j=0;

    /* compute means */
    memset(means,0,sizeof(double)* f->nlevels);
    for(j = 0; j < N; j++) {
      means[gp[j]-1] += v[j];
    }

    for(j = 0; j < f->nlevels; j++) {
      means[j] *= f->invgpsize[j];
      sum += means[j]*means[j];
    }

    for(j = 0; j < N; j++) {
      /*      sum +=  means[gp[j]-1]*means[gp[j]-1];*/
      v[j] -= means[gp[j]-1];
    }
  }

  return sqrt(sum);
}


static R_INLINE double slcentre(double *v, int N, 
		   FACTOR *factors[], int e, double *means) {
  int done[e];
  double w,best;
  int i,prev;
  double sum=0.0;
  int iter;
  FACTOR *f;
  if(e <= 2) {
    return centre(v,N,factors,e,means);
  }

  /*
    Factors with few number of levels should be centred on more often,
    since they typically will "mix" more.  Just a hunch.
    Use Sainte-Lagues method, inverse number of levels is vote. But don't allow
    the same factor in two consecutive draws. */
  w = 0;
  prev = -1;
  iter = 1;
  for(i = 0; i < e; i++) done[i] = 0;
  sum = 0.0;
  while(1) {
    int alldone;
    int ibest;
    int  j;
    const int *gp;
    /* Pick the highest score, different from previous */
    best = 0.0;
    for(i = 0; i < e; i++) {
      double sc;
      if(i == prev) continue;
      /* score is V/(2*s+1) where V is 1/nlevels */
      sc = 1.0/(sqrt(1.0*factors[i]->nlevels)*(2.0*done[i]+1.0));
      if(sc > best) {
	ibest = i;
	best = sc;
      }
    }
    done[ibest]++;
    prev = ibest;
    f = factors[ibest];
    gp = f->group;
    /* Centre on f */

    /* compute means */
    memset(means,0,sizeof(double)* f->nlevels);
    for(j = 0; j < N; j++) {
      means[gp[j]-1] += v[j];
    }

    for(j = 0; j < f->nlevels; j++) {
      means[j] *= f->invgpsize[j];
      sum += means[j]*means[j];
    }

    for(j = 0; j < N; j++) {
      /*      sum +=  means[gp[j]-1]*means[gp[j]-1];*/
      v[j] -= means[gp[j]-1];
    }
    /* Continue until everyone has got a seat */
    alldone = 1;
    for(i = 0; i < e; i++) alldone = alldone && (done[i] > 0);
    if(alldone) break;
  }
  return sqrt(sum);
}

/*
  Method of alternating projections.  Input v, output res.
*/

static int demean(double *v, int N, double *res,
		  FACTOR *factors[],int e,double eps, double *means,
		  int *stop,int vecnum, LOCK_T lock) {
  double delta;
  double norm2 = 0.0;
  int i;
  double prevdelta;
  double c;
  double neps;
  int iter;
  int okconv = 1;
  time_t last, now;
  int lastiter;
  /* make a copy to result vector, centre() is in-place */

  memcpy(res,v,N*sizeof(double));

  if(1 == e) {
    centre(res,N,factors,e,means);
    return(1);
  }

  for(i = 0; i < N; i++) norm2 += res[i]*res[i];
  norm2 = sqrt(norm2);
  prevdelta = 2*norm2;
  neps = norm2*eps;
  iter = 0;
  last = time(NULL);
  lastiter = 0;
  do {
    iter++;
    delta = centre(res,N,factors,e,means);
    c = delta/prevdelta;
    /* c is the convergence rate, at infinity they add up to 1/(1-c).
       This is true for e==2, but also in the ballpark for e>2 we guess.
    */
    if(c >= 1.0) {
      /* Seems warning() has problems with threads, do printf instead */
      LOCK(lock);
      REprintf("Demeaning of vec %d failed after %d iterations, 1-c=%.0e, delta=%.1e\n",
	       vecnum,iter,1.0-c,delta);
      UNLOCK(lock);
      okconv = 0;
      break;
    }
    now = time(NULL);
    if(now - last >= 3600 && delta > 0.0) {
      int reqiter;
      double eta;
      char tu;
      reqiter = log(neps*(1.0-c)/delta)/log(c);
      eta = 1.0*(now-last)*reqiter/(iter-lastiter);
      tu = 's';
      if(eta > 3600.0) {
	eta /= 3600.0;
	tu='h';
      }
      LOCK(lock);
      REprintf("...centering vec %d iter %d, 1-c=%.1e, delta=%.1e(%.1e), ETA in %.1f%c\n",
	       vecnum,iter,1.0-c,delta,neps*(1.0-c),eta,tu);
      UNLOCK(lock);
      lastiter = iter;
      last = now;
    }
#ifdef NOTHREADS
    R_CheckUserInterrupt();
#else
    if(*stop != 0) return(0);
#endif
    prevdelta = delta;
  } while(delta > neps*(1.0-c));
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
  LOCK_T lock;
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
  int vecnum;
  int okconv = 0;
  int myid;
  /* Find the preallocated means buffer */
  LOCK(arg->lock);
  myid = arg->threadnum++;
  UNLOCK(arg->lock);
  means = arg->means[myid];
  while(1) {
    time_t now;
    /* Find next vector */
    LOCK(arg->lock);
    vecnum = (arg->nowdoing)++;
    UNLOCK(arg->lock);
    if(vecnum >= arg->K) break;

    okconv = demean(arg->v[vecnum],arg->N,arg->res[vecnum],arg->factors,arg->e,arg->eps,means,&arg->stop,vecnum+1,arg->lock);
    LOCK(arg->lock);
    now = time(NULL);
    if(!okconv) {
      arg->badconv++;
    }
    (arg->done)++;
    if(arg->quiet > 0 && now > arg->last + arg->quiet && arg->K > 1) {
      REprintf("...finished centering vector %d of %d in %d seconds\n",
	     arg->done,arg->K,(int)(now-arg->start));
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
		      FACTOR *factors[], int e, double eps,int cores, int quiet) {
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
  for(thr = 0; thr < numthr; thr++) {
    arg.means[thr] = (double*) R_alloc(maxlev,sizeof(double));
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

#ifdef NOTHREADS
  demeanlist_thr((void*)&arg);
#else
  /* Do it in separate threads */
  for(thr = 0; thr < numthr; thr++) {
#ifdef WIN
    threads[thr] = CreateThread(NULL,0,demeanlist_thr,&arg,0,&threadids[thr]);
    if(0 == threads[thr]) error("Failed to create thread");
#else
    int stat = pthread_create(&threads[thr],NULL,demeanlist_thr,&arg);
    if(0 != stat) error("Failed to create thread, stat=%d",stat);
#endif
  }

  /* wait for completion */
  /* We want to check for interrupts regularly, and
     set a stop flag */
  while(1) {
    /* Serialize R-calls with a lock */
    LOCK(arg.lock);
    if(arg.stop == 0 && checkInterrupt()) {
      REprintf("...stopping centering threads...\n");
      arg.stop=1;
    }
    UNLOCK(arg.lock);
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
  curoff = Calloc(nlev+1,int);

  memset(f->ii,0,sizeof(int)*(nlev+1));

  /* find sizes of groups */
  for(i = 0; i < N; i++) {
    f->ii[f->group[i]]++;
  }

  /* cumulative */
  for(i = 1; i < nlev+1; i++) {
    f->ii[i] += f->ii[i-1];
  }

  for(i = 0; i < N; i++) {
    int gp = f->group[i]-1;
    f->gpl[f->ii[gp]+curoff[gp]] = i;
    curoff[gp]++;
  }
  Free(curoff);
}


static double kaczmarz(FACTOR *factors[],int e,int N, double *R,double *x,
		       double eps, double *newR, int *indices, int *stop, LOCK_T lock) {
  /* The factors define a matrix D, we will solve Dx = R
     There are N rows in D, each row a_i contains e non-zero entries, one
     for each factor, the level at that position.
     We iterate on k, start with x=0, an iteration consists
     of adding a multiple of row i (with i = k %% N), the multiplying coefficient
     is (R[i] - (a_i,x))/e 

     To get better memory locality, we create an (e X N)-matrix with the non-zero indices
  */

  double einv = 1.0/e;
  double norm2;
  double prevdiff,neweps;
  double c,sd;
  int iter=0;
  int i;
  int ie;
  int *prev,*this;
  int newN;

  /* We should remove duplicates, at least when they're consecutive.
     Non-consecutives duplicates are just alternating projections and may
     not be too inefficient.  Besides, removing them is more work...
  */

  prev = (int*) malloc(e*sizeof(int));
  this = (int*) malloc(e*sizeof(int));
  newN = 0;
  ie = 0;
  for(i = 0; i < N; i++) {
    int j;
    for(j = 0; j < e; j++) {
      this[j] = factors[j]->group[i];
    }
    if(memcmp(this,prev,e*sizeof(int)) != 0) {
      int nlev=0;
      /* not duplicate, store in indices */
      
      for(j = 0; j < e; j++) {
	indices[ie+j] = this[j]-1 + nlev;
	nlev += factors[j]->nlevels;
      }
      newR[newN] = R[i];
      newN++;
      ie += e;
      memcpy(prev,this,e*sizeof(int));
    }
  }
  free(this);
  free(prev);

  /* 
     At this point we should perhaps randomly shuffle the equations.
     We don't know when this ought to be done, but we've only seen it when
     there are many factors.
     The unif_rand isn't thread-safe, so protect with mutex.
     The Get/PutRNGstate is done in the main thread.
  */
  LOCK(lock);
  if(e > 2) {
    for(i = newN-1; i > 0; i--) {
      double dtmp;
      int k,j;
      /* Pick j between 0 and i inclusive */
      j = (int) floor((i+1) * unif_rand());
      if(j == i) continue;
      /* exchange newR[i] and newR[j]
	 as well as indices[i*e:i*e+e-1] and indices[j*e:j*e+e-1]
      */
      dtmp = newR[j];
      newR[j] = newR[i];
      newR[i] = dtmp;
      for(k = 0; k < e; k++) {
	int itmp;
	itmp = indices[j*e+k];
	indices[j*e+k] = indices[i*e+k];
	indices[i*e+k] = itmp;
      }
    }
  }
  UNLOCK(lock);

  /* Then, do the Kaczmarz iterations */
  norm2 =0.0;
  for(i = 0; i < newN; i++) norm2 += newR[i]*newR[i];
  norm2 = sqrt(norm2);
  prevdiff = 2*norm2;
  neweps = eps*norm2;
  do {
    int ie = 0; /* equals i*e; integer multiplication is slow, keep track instead */
    double diff = 0.0;
    for(i = 0; i < newN; i++,ie+=e) {
      double upd = 0.0;
      int j;
      upd = newR[i];

      /* Subtract inner product */
      for(j = 0; j < e; j++) {
	int idx = indices[ie + j];
	upd -= x[idx];
      }
      upd *= einv;

      /* Update x */
      for(j = 0; j < e; j++) {
	int idx = indices[ie + j];      
	x[idx] += upd;
      }

      /* Keep track of update */
      diff += upd*upd;
    }

    /* Relax a bit */
    /*    R_CheckUserInterrupt();*/
    sd = sqrt(diff*e);
    c = sd/prevdiff;
    prevdiff = sd;
    iter++;
    if(c >= 1.0) {
      REprintf("Kaczmarz failed in iter %d*%d, 1-c=%.0e, delta=%.1e, eps=%.1e\n",iter,newN,1.0-c,sd,neweps);
      break;
    }
#ifdef NOTHREADS
    R_CheckUserInterrupt();
#else
    if(*stop != 0) return(0);
#endif
  } while(sd >= neweps*(1.0-c)) ;

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
  double **newR;
  int **indices;
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
		    arg->newR[myid],arg->indices[myid], &arg->stop, arg->lock);
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
  int N = 0;
  int i;
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

  
  if(!isNull(initial)) {
    init = REAL(initial);
  }
  PROTECT(flist = AS_LIST(flist));
  numfac = LENGTH(flist);

  factors = (FACTOR**) R_alloc(numfac,sizeof(FACTOR*));
  for(i = 0; i < numfac; i++) {
    int len;
    FACTOR *f;
    len = LENGTH(VECTOR_ELT(flist,i));
    if(N == 0) N = len;
    if(len != N) {
       error("Factors must have the same length %d %d",len,N);
    }

    factors[i] = (FACTOR*) R_alloc(1,sizeof(FACTOR));
    f = factors[i];
    f->group = INTEGER(VECTOR_ELT(flist,i));
    f->nlevels = LENGTH(getAttrib(VECTOR_ELT(flist,i),R_LevelsSymbol));
    sumlev += f->nlevels;
  }

  if(!isNull(initial) && LENGTH(initial) != sumlev)
    error("Initial vector must have length %d, but is %d\n",sumlev, LENGTH(initial));

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
	error("Vector length must be equal to factor length %d %d",LENGTH(elt),N);
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
  for(i = 0; i < listlen; i++) {
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

      int j;
      SEXP mtx;
      /* Allocate a matrix */
      PROTECT(mtx = allocMatrix(REALSXP,sumlev,cols));
      SET_VECTOR_ELT(reslist,i,mtx);
      UNPROTECT(1);
      /* Set up pointers */
      for(j = 0; j < cols; j++) {
	vectors[cnt] = REAL(elt) + j*N;
	target[cnt] = REAL(mtx) + j*sumlev;
	cnt++;
      }
    }
  }


  for(cnt = 0; cnt < numvec; cnt++) {
    if(init != 0)
      for(i = 0; i < sumlev; i++) target[cnt][i] = init[i];
    else
      for(i = 0; i < sumlev; i++) target[cnt][i] = 0.0;
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
  arg.newR = (double**) R_alloc(numthr,sizeof(double*));
  arg.indices = (int**) R_alloc(numthr,sizeof(int*));

#ifdef NOTHREADS
  arg.newR[0] = (double*) R_alloc(N,sizeof(double));
  arg.indices[0] = (int*) R_alloc(numfac*N,sizeof(int));
  kaczmarz_thr((void*)&arg);
#else
  /* Do it in separate threads */
  for(thr = 0; thr < numthr; thr++) {
    int stat;
    arg.newR[thr] = (double*) R_alloc(N,sizeof(double));
    arg.indices[thr] = (int*) R_alloc(numfac*N,sizeof(int));
#ifdef WIN
    threads[thr] = CreateThread(NULL,0,kaczmarz_thr,&arg,0,&threadids[thr]);
    if(0 == threads[thr]) error("Failed to create kaczmarz thread");
#else
    stat = pthread_create(&threads[thr],NULL,kaczmarz_thr,&arg);
    if(0 != stat) error("Failed to create kaczmarz thread, stat=%d",stat);
#endif
  }
  /* wait for completion */
  /* We want to check for interrupts regularly, and
     set a stop flag */
  while(1) {
    if(arg.stop == 0 && checkInterrupt()) {
      REprintf("...stopping Kaczmarz threads...\n");
      arg.stop=1;
    }

#ifdef WIN
    /*    WaitForMultipleObjects(numthr,threads,TRUE,INFINITE);*/
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

static SEXP R_demeanlist(SEXP vlist, SEXP flist, SEXP Ricpt, SEXP Reps, SEXP scores, SEXP quiet) {
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

#ifndef NOTHREADS
  R_CStackLimit = (uintptr_t)-1;
#endif

  icpt = INTEGER(Ricpt)[0] - 1; /* convert from 1-based to zero-based */
  eps = REAL(Reps)[0];
  cores = INTEGER(scores)[0];


  PROTECT(flist = AS_LIST(flist));
  numfac = LENGTH(flist);

  factors = (FACTOR**) R_alloc(numfac,sizeof(FACTOR*));
  for(i = 0; i < numfac; i++) {
    int j;
    int len;
    FACTOR *f;
    len = LENGTH(VECTOR_ELT(flist,i));
    if(i == 0) {
      N = len;
    } else if(len != N) {
      error("All factors must have the same length %d %d",len,N);
    }

    factors[i] = (FACTOR*) R_alloc(1,sizeof(FACTOR));
    f = factors[i];
    f->group = INTEGER(VECTOR_ELT(flist,i));
    f->nlevels = LENGTH(getAttrib(VECTOR_ELT(flist,i),R_LevelsSymbol));
    /* Make array for holding precomputed group levels */
    /* Now, what about entries which don't belong to a group */
    /* I.e. NA entries, how is that handled, I wonder. */
    /* seems to be negative.  Anyway, we fail on them. No, we don't */
    f->invgpsize = (double *)R_alloc(f->nlevels,sizeof(double));
    memset(f->invgpsize,0,f->nlevels*sizeof(double));
    /* first count it */
    for(j = 0; j < N; j++) {
      /* skip entries without a group, do we need that? */
      /* if(f->group[j] < 1) error("Factors can't have missing levels"); */
      if(f->group[j] > 0) {
	f->invgpsize[f->group[j]-1] += 1.0;
      } else {
	error("Factors can't have missing levels");
      }
    }
    /* then invert it, it's much faster to multiply than to divide */
    /* in the iterations */
    for(j = 0; j < f->nlevels; j++) {
      f->invgpsize[j] = 1.0/f->invgpsize[j];
    }
    /*invertfactor(f,N);*/
  }


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
	vectors[cnt] = REAL(elt) + j*N;
	target[cnt] = REAL(mtx) + (rcols++)*N;
	cnt++;
      }
    }
  }

  /* Then do stuff */
  PROTECT(badconv = allocVector(INTSXP,1));
  INTEGER(badconv)[0] = demeanlist(vectors,N,numvec,target,factors,numfac,eps,cores,INTEGER(quiet)[0]);

  if(INTEGER(badconv)[0] > 0)
    warning("%d vectors failed to centre to tolerance %.1le",INTEGER(badconv)[0],eps);
  setAttrib(reslist,install("badconv"),badconv);
  /* unprotect the reslist */
  UNPROTECT(4);
  return(reslist);
}


/* 
Then for finding connection components 
From R we take a list of factors, we return
a factor of the same length with the connection
components
*/


/*
  This one is a bit tricky.  We could do it quite elegantly recursively,
  but we would suffer from a deep call-stack.  Thus, we make our own
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
}

/*
R entry-point for conncomp.  Takes a list of factors as input.
 */

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


static R_CallMethodDef callMethods[] = {
  {"conncomp", (DL_FUNC) &R_conncomp, 1},
  {"demeanlist", (DL_FUNC) &R_demeanlist, 6},
  {"kaczmarz", (DL_FUNC) &R_kaczmarz, 5},
  {NULL, NULL, 0}
};

void attribute_visible R_init_lfe(DllInfo *info) {
  /* register our routines */
  (void)R_registerRoutines(info,NULL,callMethods,NULL,NULL);
}

