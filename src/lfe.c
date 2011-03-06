/*
 Author: Simen Gaure
 Copyright: 2011, Simen Gaure
 Licence: Artistic 2.0
*/

#if defined(_WIN32) || defined(_WIN64) || defined(WIN64)
#define WIN
#endif

#ifdef WIN
#include <windows.h>
#else
#ifndef NOTHREADS
#include <pthread.h>

#ifdef _POSIX_BARRIERS
#define COOP
#endif

/* It's slower, too much synchronization */
#undef COOP

#endif
#endif

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

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

static R_INLINE double gmean(double *v,int start, int end, FACTOR *f,double *means) {
  /*
    Do the actual work, single vector, single factor, 
    output in preallocated vector means
  */
  int i;
  int *group = f->group;
  int nlev = f->nlevels;
  
  double sum = 0.0;
  /* zero the means */
  memset(means,0,sizeof(double)*nlev);

  /* add entries into right group */

  for(i = start; i < end; i++) {
    means[group[i]-1] += v[i];
  }

  /* divide by group size */
  for(i = 0; i < nlev; i++) {
    means[i] *= f->invgpsize[i];
    sum += means[i]*means[i];
  }
  return sum;
}

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
    sum += gmean(v,0,N,f,means);

    for(j = 0; j < N; j++) {
      v[j] -= means[gp[j]-1];
    }
  }

  return sqrt(sum);
}

/*
  Method of alternating projections.  Input v, output res.
*/

static int demean(double *v, int N, double *res,
		   FACTOR *factors[],int e,double eps, double *means) {
  double delta;
  double norm2 = 0.0;
  int i;
  double prevdelta;
  double c;
  double neps;
  int iter;
  int okconv = 1;
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
  do {
    iter++;
    delta = centre(res,N,factors,e,means);
    c = delta/prevdelta;
    /* c is the convergence rate, at infinity they add up to 1/(1-c).
       This is true for e==2, but also in the ballpark for e>2 we guess.
    */
    if(c >= 1.0) {
      /* Seems warning() has problems with threads, do printf instead */
      fprintf(stderr,"Demeaning failed after %d iterations, c-1=%.0e, delta=%.1e\n",iter,c-1.0,delta);
      okconv = 0;
      break;
    }
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
#ifdef COOP
  double **means;
  double sum;
  int threads;
  int tmp;
  double *gmeans;
  pthread_barrier_t *barrier;
#endif
#ifndef NOTHREADS
#ifdef WIN
  HANDLE *lock;
#else
  pthread_mutex_t *lock;
#endif
#endif
} PTARG;


/*
  Thread routine for demeaning.
  Threads pick a vector to demean by looking
  in arg->nowdoing, a common counter.  Proteced
  by a mutex.
 */

#ifdef COOP
#define BARRIER(b) pthread_barrier_wait(b)
#else
#define BARRIER(b)
#endif

#ifdef NOTHREADS
#define LOCK(l)
#define UNLOCK(l)
#else
#ifdef WIN
#define LOCK(l) WaitForSingleObject(l,INFINITE)
#define UNLOCK(l) ReleaseMutex(l)
#else
#define LOCK(l) (void)pthread_mutex_lock(l)
#define UNLOCK(l) (void)pthread_mutex_unlock(l)
#endif
#endif


/* Cooperative centering 
   Each vector is centered by all the threads.
   Each thread gets a portion of the vector.
   This is faster with few vectors, but there's 
   synchronization which may cause it to slow down.
   Indeed it's slow.  Deactivated.

   Haven't bothered to figure out how to do barriers in windows yet,
   so skip the cooperation threading there.
*/

#ifdef COOP

#ifdef WIN
DWORD WINAPI demeanlist_coop(LPVOID varg) {
#else
static void *demeanlist_coop(void *varg) {
#endif
  PTARG *arg = (PTARG*) varg;
  int maxlev;
  int i;
  int myid;
  double *means;
  int N = arg->N;
  int e = arg->e;
  int vecstart,veclen,vecend;
  /* preallocate means */
  maxlev = 0;
  for(i = 0; i < e; i++)
    if(maxlev < arg->factors[i]->nlevels) maxlev = arg->factors[i]->nlevels;

  /* Figure out which slice of the vectors I'm going to compute on */
  /* First, number the threads and store the means pointer */

  LOCK(arg->lock);
  myid = (arg->tmp)++;
  means = Calloc(maxlev,double);
  arg->means[myid] = means;
  UNLOCK(arg->lock);

  BARRIER(arg->barrier);
  veclen = N / arg->threads;
  vecstart = veclen*myid;
  /* Stuff the rest into the last */
  if(myid == arg->threads-1) {
    veclen = N - vecstart;
  }

  vecend = vecstart + veclen;
  for(i = 0; i < arg->K; i++) {
    time_t start = time(NULL);
    double sd;
    double c;
    double prevdiff;
    double neps;
    double *v = arg->v[i];
    double *t = arg->res[i];
    int k;
    /* Centre this vector */
    /* copy to target */
    neps = 0;
    for(k = 0; k < N; k++) {
      neps += v[k]*v[k];
    }
    prevdiff = 2.0*sqrt(neps);
    neps = arg->eps*sqrt(neps);
    if(myid == 0) {
      memcpy(t,v,N*sizeof(double));	
      arg->sum = 0;
    }
    /* All wait for the copy to the common area */
    BARRIER(arg->barrier);
    
    do {
      int j;
      /* Compute the means */
      /* Zero the square sums */
      if(BARRIER(arg->barrier) == PTHREAD_BARRIER_SERIAL_THREAD) {
	arg->sum = 0;
      }
      BARRIER(arg->barrier);

      for(j = 0; j < e; j++) {
	FACTOR *f = arg->factors[j];
	int *gp = f->group;
	int idx;

	gmean(t,vecstart,vecend,f,means);
	/* All wait for the means */
	if(BARRIER(arg->barrier) == PTHREAD_BARRIER_SERIAL_THREAD) {
	  int k,l;
	  double sum = 0;
	  /* sum the means from all threads, we're in only one now */
	  memset(arg->gmeans,0, f->nlevels*sizeof(double));
	  for(l = 0; l < arg->threads; l++) {
	    for(k = 0; k < f->nlevels; k++) {
	      arg->gmeans[k] += arg->means[l][k];
	    }
	  }
	  /* add to the square sum */
	  for(k = 0; k < f->nlevels; k++) {
	    sum += arg->gmeans[k] * arg->gmeans[k];
	  }
	  /* store it */
	  arg->sum += sum;
	}
	/* all wait for the means */
	BARRIER(arg->barrier);
	/* subtract the means, own slice */
	for(idx = vecstart; idx < vecend; idx++) {
	  t[idx] -= arg->gmeans[gp[idx]-1];
	}
      }
      BARRIER(arg->barrier);
      sd = sqrt(arg->sum);
      c = sd/prevdiff;
      if(c >= 1.0) {
	if(myid == 0) {
	  arg->badconv++;
	}
	break;
      }
      prevdiff = sd;
    } while(sd > neps*(1.0-c));
    now = time(NULL);
    arg->done++;
    if(now > arg->last + 300) {
      printf("...finished centering vector %d of %d in %d seconds\n",
	      i+1,arg->K,now-arg->start);
      arg->last = now;
    }
  }
  LOCK(arg->lock);
  Free(means);
  UNLOCK(arg->lock);
  return 0;
}
#define demeanlist_thr demeanlist_coop
#else
#define demeanlist_thr demeanlist_single
#endif


/* The thread routine 
   Each thread centres an entire vector.
*/
#ifdef WIN
DWORD WINAPI demeanlist_single(LPVOID varg) {
#else
static void *demeanlist_single(void *varg) {
#endif
  PTARG *arg = (PTARG*) varg;
  double *means;
  int i;
  int maxlev=0;
  int vecnum;
  int okconv = 0;

  /* preallocate means */
  for(i = 0; i < arg->e; i++)
    if(maxlev < arg->factors[i]->nlevels) maxlev = arg->factors[i]->nlevels;
  means = (double*) malloc(maxlev*sizeof(double));
  while(1) {
    time_t start,now;
    /* Find next vector */
    LOCK(arg->lock);
    vecnum = (arg->nowdoing)++;
    UNLOCK(arg->lock);
    if(vecnum >= arg->K) break;

    start = time(NULL);
    okconv = demean(arg->v[vecnum],arg->N,arg->res[vecnum],arg->factors,arg->e,arg->eps,means);
    now = time(NULL);
    LOCK(arg->lock);
    if(!okconv) {
      arg->badconv++;
    }
    (arg->done)++;
    if(!arg->quiet && now > arg->last + 300 && arg->K > 1) {
      printf("...finished centering vector %d of %d in %d seconds\n",
	     arg->done,arg->K,(int)(now-arg->start));
      arg->last = now;
    }
    UNLOCK(arg->lock);
  }
  free(means);
  return 0;
}


/* main C entry-point for demeaning.  Called from the R-entrypoint */
static int demeanlist(double **vp, int N, int K, double **res,
		      FACTOR *factors[], int e, double eps,int cores, int quiet) {
  PTARG arg;
#ifndef NOTHREADS
  int numthr = 1;
  int i;
  int thr;
  int maxlev;
#ifdef WIN
  HANDLE *threads;
  DWORD *threadids;
  HANDLE lock;
#else
  pthread_t *threads;
  pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
#endif
#ifdef COOP
  pthread_barrier_t barrier;
#endif

  maxlev = 0;
  for(i = 0; i < e; i++)
    if(maxlev < factors[i]->nlevels) maxlev = factors[i]->nlevels;


  numthr = cores;
  if(numthr > K) numthr = K;
  if(numthr < 1) numthr = 1;
#ifdef WIN
  lock = CreateMutex(NULL,FALSE,NULL);
  threads = (HANDLE*) R_alloc(numthr,sizeof(HANDLE));
  threadids = (DWORD*) R_alloc(numthr,sizeof(DWORD));
#else
  threads = (pthread_t*) R_alloc(numthr,sizeof(pthread_t));
#endif

  arg.lock = &lock;
#ifdef COOP
  arg.barrier = &barrier;
  pthread_barrier_init(arg.barrier,NULL,numthr);
#endif
#endif

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
#ifdef COOP
  arg.means = (double**) R_alloc(numthr,sizeof(double*));
  arg.gmeans = (double*) R_alloc(maxlev,sizeof(double));
  arg.threads = numthr;
  arg.tmp = 0;
#endif

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

#ifdef WIN
  WaitForMultipleObjects(numthr,threads,TRUE,INFINITE);
  CloseHandle(lock);
  for(thr = 0; thr < numthr; thr++) {
    CloseHandle(threads[thr]);
  }
#else
  for(thr = 0; thr < numthr; thr++) {
    (void)pthread_join(threads[thr], NULL);
  }
#ifdef COOP
  pthread_barrier_destroy(arg.barrier);
#endif
#endif
#endif
  if(!quiet && arg.start != arg.last && K > 1)
    printf("...%d vectors centred in %d seconds\n",K,(int)(time(NULL)-arg.start));
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
		       double eps, double *newR, int *indices) {
  /* The factors define a matrix D, we will solve Dx = R
     There are N rows in D, each row a_i contains e non-zero entries, one
     for each factor, the level at that position.
     We iterate on k, start with x=0, an iteration consists
     of adding a multiple of row i (with i = k %% N), the multiplying coefficient
     is (R[i] - (a_i,x))/e 

     To get better memory locality, we create an (e X N)-matrix with the non-zero indices
  */

  /*  int *indices = (int*) R_alloc(e*N,sizeof(int));*/
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
      fprintf(stderr,"Kaczmarz failed in iter %d*%d, c-1=%.0e, sd=%.1e, eps=%.1e\n",iter,newN,c-1.0,sd,neweps);
      /*      warning("Kaczmarz failed in iter %d*%d, c-1=%.0le, sd=%.1le, eps=%.1le\n",iter,newN,c-1.0,sd,neweps);*/
      break;
    }
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
  int numvec;
  double eps;
  double *tol;
  double **newR;
  int **indices;
#ifndef NOTHREADS
#ifdef WIN
  HANDLE *lock;
#else
  pthread_mutex_t *lock;
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
		    arg->newR[myid],arg->indices[myid]);
  }
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
#ifndef NOTHREADS
  int numthr = 1;
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
#ifdef WIN
  lock = CreateMutex(NULL,FALSE,NULL);
  threads = (HANDLE*) R_alloc(numthr,sizeof(HANDLE));
  threadids = (DWORD*) R_alloc(numthr,sizeof(DWORD));
#else
  threads = (pthread_t*) R_alloc(numthr,sizeof(pthread_t));
#endif


  arg.factors = factors;
  arg.source = vectors;
  arg.target = target;
  arg.nowdoing = 0;
  arg.threadnum = 0;
  arg.lock = &lock;
  arg.e = numfac;
  arg.eps = eps;
  arg.numvec = numvec;
  arg.N = N;
  arg.newR = (double**) R_alloc(numthr,sizeof(double*));
  arg.indices = (int**) R_alloc(numthr,sizeof(int*));

#ifdef NOTHREADS
  arg.newR[0] = R_alloc(N,sizeof(double));
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

#ifdef WIN
  WaitForMultipleObjects(numthr,threads,TRUE,INFINITE);
  CloseHandle(lock);
  for(thr = 0; thr < numthr; thr++) {
    CloseHandle(threads[thr]);
  }
#else
  for(thr = 0; thr < numthr; thr++) {
    (void)pthread_join(threads[thr], NULL);
  }
#endif
#endif
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

