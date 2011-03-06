/*
 Author: Simen Gaure
 Copyright: 2011, Simen Gaure
 Licence: Artistic 2.0
*/

#if defined(_WIN32) || defined(_WIN64) || defined(WIN64)
#define WIN
#endif

#include <stdlib.h>
#ifdef WIN
#include <windows.h>
#else
#ifndef NOTHREADS
#include <pthread.h>
#endif
#endif

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


static double gmean(double *v,int N, FACTOR *f,double *means) {
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

  for(i = 0; i < N; i++) {
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

static double centre(double *v, int N, 
		   FACTOR *factors[], int e, double *means) {

  int i;
  double sum = 0.0;
  
  for(i=0; i < e; i++) {
    FACTOR *f = factors[i];
    const int *gp = f->group;
    int j=0;
    /* compute means */
    sum += gmean(v,N,f,means);

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
      REprintf("Demeaning failed after %d iterations, c-1=%.0le, delta=%.1le\n",iter,c-1.0,delta);
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


/* The thread routine */
#ifdef WIN
DWORD WINAPI demeanlist_thr(LPVOID varg) {
#else
static void *demeanlist_thr(void *varg) {
#endif
  PTARG *arg = (PTARG*) varg;
  double *means;
  int i;
  int maxlev;
  int vecnum;
  int okconv = 0;
  /* preallocate means */
  maxlev = 0;
  for(i = 0; i < arg->e; i++)
    if(maxlev < arg->factors[i]->nlevels) maxlev = arg->factors[i]->nlevels;
  means = Calloc(maxlev,double);
  while(1) {
    LOCK(arg->lock);
    vecnum = (arg->nowdoing)++;
    UNLOCK(arg->lock);
    if(vecnum >= arg->K) break;
    okconv = demean(arg->v[vecnum],arg->N,arg->res[vecnum],arg->factors,arg->e,arg->eps,means);
    if(!okconv) {
      LOCK(arg->lock);
      arg->badconv++;
      UNLOCK(arg->lock);
    }
  }

  Free(means);
  return 0;
}


/* main C entry-point for demeaning.  Called from the R-entrypoint */
static int demeanlist(double **vp, int N, int K, double **res,
		       FACTOR *factors[], int e, double eps,int cores) {
  PTARG arg;

#ifndef NOTHREADS
  int thr;
  int numthr;

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

#ifdef WIN
  lock = CreateMutex(NULL,FALSE,NULL);
  threads = (HANDLE*) R_alloc(numthr,sizeof(HANDLE));
  threadids = (DWORD*) R_alloc(numthr,sizeof(DWORD));
#else
  threads = (pthread_t*) R_alloc(numthr,sizeof(pthread_t));
#endif

  arg.lock = &lock;
#endif

  arg.badconv = 0;
  arg.N = N;
  arg.v = vp;
  arg.res = res;
  arg.K = K;
  arg.nowdoing = 0;
  arg.factors = factors;
  arg.e = e;
  arg.eps = eps;

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
#endif
#endif

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


static double kaczmarz(FACTOR *factors[],int e,int N, double *R,double *x, double eps) {
  /* The factors define a matrix D, we will solve Dx = R
     There are N rows in D, each row a_i contains e non-zero entries, one
     for each factor, the level at that position.
     We iterate on k, start with x=0, an iteration consists
     of adding a multiple of row i (with i = k %% N), the multiplying coefficient
     is (R[i] - (a_i,x))/e 

     To get better memory locality, we create an (e X N)-matrix with the non-zero indices
  */

  int *indices = (int*) R_alloc(e*N,sizeof(int));
  double einv = 1.0/e;
  double norm2;
  double prevdiff,neweps;
  double c,sd;
  int iter=0;
  int i;
  int ie;
  int *prev,*this;
  double *newR;
  int newN;

  /* We should remove duplicates, at least when they're consecutive.
     Non-consecutives duplicates are just alternating projections and may
     not be too inefficient.  Besides, removing them is more work...
  */

  newR = (double*) R_alloc(N,sizeof(double));
  prev = Calloc(e,int);
  this = Calloc(e,int);
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
  Free(this);
  Free(prev);


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
    R_CheckUserInterrupt();
    sd = sqrt(diff*e);
    c = sd/prevdiff;
    prevdiff = sd;
    iter++;
    if(c >= 1.0) {
      REprintf("Kaczmarz failed in iter %d*%d, c-1=%.0le, sd=%.1le\n",iter,newN,c-1.0,sd);
      warning("Kaczmarz failed in iter %d*%d, c-1=%.0le, sd=%.1le\n",iter,newN,c-1.0,sd);
      break;
    }
  } while(sd >= neweps*(1.0-c)) ;
  return(sd);
}


static SEXP R_kaczmarz(SEXP flist, SEXP RR, SEXP Reps) {

  double eps = REAL(Reps)[0];
  double *R = REAL(RR);
  FACTOR **factors;
  SEXP result;
  int numfac;
  int N = LENGTH(RR);
  int i;
  int sumlev = 0;
  SEXP tol;
  /*  double norm2 = 0.0;*/
  PROTECT(flist = AS_LIST(flist));
  numfac = LENGTH(flist);

  factors = (FACTOR**) R_alloc(numfac,sizeof(FACTOR*));
  for(i = 0; i < numfac; i++) {
    int len;
    FACTOR *f;
    len = LENGTH(VECTOR_ELT(flist,i));
    if(len != N) {
       error("Factors and vector must have the same length %d %d",len,N);
    }

    factors[i] = (FACTOR*) R_alloc(1,sizeof(FACTOR));
    f = factors[i];
    f->group = INTEGER(VECTOR_ELT(flist,i));
    f->nlevels = LENGTH(getAttrib(VECTOR_ELT(flist,i),R_LevelsSymbol));
    sumlev += f->nlevels;
  }

  PROTECT(result = allocVector(REALSXP,sumlev));
  for(i = 0; i < sumlev; i++) REAL(result)[i] = 0.0;
  PROTECT(tol = allocVector(REALSXP,1));
  REAL(tol)[0] = kaczmarz(factors,numfac,N,R,REAL(result),eps);
  setAttrib(result,install("tolerance"),tol);
  UNPROTECT(3);
  return(result);
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

static SEXP R_demeanlist(SEXP vlist, SEXP flist, SEXP Ricpt, SEXP Reps, SEXP scores) {
  int numvec;
  int numfac;
  double **vectors;
  double **target;
  double eps;
  SEXP resvec;
  SEXP reslist;
  FACTOR **factors;
  int cores;
  int i;
  int N=0, len;
  int mtinput = 0;
  int icpt;
  SEXP badconv;
  icpt = INTEGER(Ricpt)[0] - 1; /* convert from 1-based to zero-based */
  eps = REAL(Reps)[0];
  cores = INTEGER(scores)[0];
  if(!IS_NUMERIC(vlist)) {
    PROTECT(vlist = AS_LIST(vlist));
    numvec = LENGTH(vlist);

    vectors = (double **)R_alloc(numvec,sizeof(double*));
    /* Create the vector pointers */
    for(i = 0; i < numvec; i++) {
      len = LENGTH(VECTOR_ELT(vlist,i));
      if(i > 0 && len != N) {
	error("Vectors must have equal length");
      }
      vectors[i] = REAL(VECTOR_ELT(vlist,i));
      N = len;
    }
  } else {
    int ncol;
    mtinput = 1;
    PROTECT(vlist = vlist);
    N = nrows(vlist);
    ncol = ncols(vlist);

    vectors = (double **)R_alloc(ncol,sizeof(double*));
    /* Create the vector pointers, skip intercept */
    numvec = 0;
    for(i = 0; i < ncol; i++) {
      /* to the i'th column */
      if(i == icpt) continue;
      vectors[numvec++] = REAL(vlist) + i*N;
    }
  } 

  PROTECT(flist = AS_LIST(flist));
  numfac = LENGTH(flist);

  factors = (FACTOR**) R_alloc(numfac,sizeof(FACTOR*));
  for(i = 0; i < numfac; i++) {
    int j;
    FACTOR *f;
    len = LENGTH(VECTOR_ELT(flist,i));
    if(len != N) {
       error("Factors and vectors must have the same length %d %d",len,N);
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
      if(f->group[j] > 0)
	f->invgpsize[f->group[j]-1] += 1.0;
    }
    /* then invert it, it's much faster to multiply than to divide */
    /* in the iterations */
    for(j = 0; j < f->nlevels; j++) {
      f->invgpsize[j] = 1.0/f->invgpsize[j];
    }
    /*invertfactor(f,N);*/
  }

  /* allocate result structure */

  target = (double**) R_alloc(numvec,sizeof(double*));
  if(mtinput) {
    double *mat;
    PROTECT(reslist = allocMatrix(REALSXP,N,numvec));
    mat = REAL(reslist);
    for(i = 0; i < numvec; i++) {
      target[i] = mat + i*N;
    }
  } else {
    PROTECT(reslist = NEW_LIST(numvec));

    for(i = 0; i < numvec; i++) {
      PROTECT(resvec = allocVector(REALSXP,N));  /* need protect? */
      target[i] = REAL(resvec);
      SET_VECTOR_ELT(reslist,i,resvec);
      /* as member of list it's protected */
      UNPROTECT(1);
    }
  }  
  PROTECT(badconv = allocVector(INTSXP,1));
  INTEGER(badconv)[0] = demeanlist(vectors,N,numvec,target,factors,numfac,eps,cores);
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
  
  UNPROTECT(1);
  /* unprotect flist  */
  UNPROTECT(1);
  return(result);
}


static R_CallMethodDef callMethods[] = {
  {"conncomp", (DL_FUNC) &R_conncomp, 1},
  {"demeanlist", (DL_FUNC) &R_demeanlist, 5},
  {"kaczmarz", (DL_FUNC) &R_kaczmarz, 3},
  {NULL, NULL, 0}
};

void attribute_visible R_init_lfe(DllInfo *info) {
  /* register our routines */
  (void)R_registerRoutines(info,NULL,callMethods,NULL,NULL);
}

