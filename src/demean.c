/*
  $Id: demean.c 1671 2015-03-23 13:04:42Z sgaure $
*/
#include "lfe.h"
/* Need sprintf */
#include <stdio.h>  

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
  double *weights;
  int *scale;
  LOCK_T lock;
  int gkacc;
#ifndef NOTHREADS
#ifndef WIN
  int running;
#ifdef HAVE_SEM
  sem_t finished;
#endif
#endif
#endif
} PTARG;



/*
  Centre on all factors in succession.  Vector v. In place.
 */

static R_INLINE void centre(double *v, int N, 
			    FACTOR *factors[], int e, double *means,
			    double *weights) {

  const int hw = (weights != NULL);
  for(int i=0; i < e; i++) {
    FACTOR *f = factors[i];
    const int *gp = f->group;
    const int hx = (f->x != NULL);
    int j=0;

    /* compute means */
    memset(means,0,sizeof(double)* f->nlevels);
    for(j = 0; j < N; j++) {
      double w = hw ? (hx ? f->x[j]*weights[j] : weights[j]) : (hx ? f->x[j] : 1.0);
      if(gp[j] > 0) means[gp[j]-1] += v[j]*w;
    }

    for(j = 0; j < f->nlevels; j++) {
      means[j] *= f->invgpsize[j];
    }

    /* subtract means */
    for(j = 0; j < N; j++) {
      double w = hw ? (hx ? f->x[j]*weights[j] : weights[j]) : (hx ? f->x[j] : 1.0);
      if(gp[j] > 0) v[j] -= means[gp[j]-1]*w;
    }
  }
}


/*
  Method of alternating projections.  Input v, output res.
*/

static int demean(double *v, int N, double *res, double *weights,int *scale,
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
  // zero out NaNs
  for(int i = 0; i < N; i++) {
    if(isnan(res[i])) res[i]=0.0;
  }
  if(weights != NULL && scale[0]) for(int i = 0; i < N; i++) res[i] = res[i]*weights[i];

  centre(res, N, factors, e, means, weights);

  if(e <= 1 || factors[0]->oneiter) {
    if(weights != NULL && scale[1]) for(int i = 0; i < N; i++) res[i] = res[i]/weights[i];
    return(1);
  }
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
  double prevdp=0.0;
  do {
    iter++;
    if(gkacc) memcpy(prev, res, N*sizeof(double));
    centre(res,N,factors,e,means, weights);

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

    delta = 0.0;
    for(int i = 0; i < N; i++) delta += (prev2[i]-res[i])*(prev2[i]-res[i]);

    // make this a power of two, so we don't have to do integer division
    // Check convergence rate every now and then
    // For the purpose of computing time to convergence, we assume convergence is linear, 
    // i.e. that the decrease in norm since the previous iteration is a constant factor c.
      // compute delta per iter
    memcpy(prev2,res,N*sizeof(double));
      // delta is the square norm improvement since last time
      // we normalize it to be per iteration
      // we divide it by the number of iterations to get an improvement per iteration
    delta = sqrt(delta);
      // then we compute how fast the improvement dimishes. We use this to predict when we're done
      // for e == 2 without gkacc it should diminish, so that c < 1, but for e > 2 it may
      // increase, in particular when we do acceleration. Then it is difficult to predict when we will be done
    c = delta/prevdelta;
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
      char buf[256];
      // Now, use a c which is based on the medium c since the last print
      if(prevdp == 0.0)
	c = delta/prevdelta;
      else
	c = pow(delta/prevdp, 1.0/(iter-lastiter));
      reqiter = log(target/delta)/log(c);
      eta = 1.0*(now-last)*reqiter/(iter-lastiter);
      if(eta < 0) eta = NA_REAL; 
      if(gkacc&&0) {
	sprintf(buf,"...centering vec %d iter %d, delta=%.1e(target %.1e)\n",
		vecnum,iter,delta,target);
	
      } else {
	time_t arriv = now + (time_t) eta;
	char timbuf[50];
	struct tm tmarriv;
#ifdef WIN
	localtime_s(&tmarriv, &arriv);
#else
	localtime_r(&arriv,&tmarriv);
#endif
	strftime(timbuf, sizeof(timbuf), "%c", &tmarriv);
	//	  ctime_r(&arriv, timbuf);
	sprintf(buf,"...centering vec %d i:%d c:%.1e d:%.1e(t:%.1e) ETA:%s\n",
		vecnum,iter,1.0-c,delta,target,timbuf);
	
      }
      pushmsg(buf,lock);
      lastiter = iter;
      prevdp = delta;
      last = now;
    }
    
    prevdelta = delta;    

#ifdef NOTHREADS
    R_CheckUserInterrupt();
#else
    int stp;
    LOCK(lock); stp = *stop; UNLOCK(lock);
    if(stp != 0) {okconv = 0; break;}
#endif

  } while(delta > target);
  if(weights != NULL && scale[1]) for(int i = 0; i < N; i++) res[i] = res[i]/weights[i];
  return(okconv);
}



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
#ifdef HAVE_THREADNAME
    char thrname[16];
    snprintf(thrname,16, "Ct %d/%d",vecnum+1, arg->K);
    STNAME(thrname);
#endif
    okconv = demean(arg->v[vecnum],arg->N,arg->res[vecnum],arg->weights,arg->scale,
		    arg->factors,arg->e,arg->eps,means,
		    tmp1,tmp2,
		    &arg->stop,vecnum+1,
		    arg->lock,arg->gkacc);
    now = time(NULL);
    LOCK(arg->lock);
    if(!okconv) {
      arg->badconv++;
    }
    (arg->done)++;
    if(arg->quiet > 0 && now > arg->last + arg->quiet && arg->K > 1) {
      char buf[256];
      sprintf(buf,"...finished centering vector %d of %d in %d seconds\n",
	     arg->done,arg->K,(int)(now-arg->start));
      arg->last = now;
      UNLOCK(arg->lock); 
      pushmsg(buf,arg->lock);
    } else {
      UNLOCK(arg->lock);
    }
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
static int demeanlist(double **vp, int N, int K, double **res, double *weights,int *scale,
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
  initmsg();
#endif

  numthr = cores;
  if(numthr > K) numthr = K;
  if(numthr < 1) numthr = 1;

#ifdef WIN
  lock = CreateMutex(NULL,FALSE,NULL);
  if(lock == NULL) {
    error("Couldn't create mutex (error=%d)", GetLastError());
  }
  threads = (HANDLE*) R_alloc(numthr,sizeof(HANDLE));
  threadids = (DWORD*) R_alloc(numthr,sizeof(DWORD));
  arg.lock = lock;
#else
  threads = (pthread_t*) R_alloc(numthr,sizeof(pthread_t));
  arg.lock = &lock;
#ifdef HAVE_SEM
  if(sem_init(&arg.finished,0,0) != 0) error("sem_init failed, errno=%d",errno);
#endif
  arg.running = numthr;
#endif
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

  arg.weights = weights;
  arg.scale = scale;
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
      break;
    }
#else
#ifndef HAVE_SEM
    struct timespec atmo = {0,50000000};
    /* Kludge when no timedwait, i.e. MacOSX */
    int done;
    if(arg.stop == 0) nanosleep(&atmo,NULL);
    LOCK(arg.lock);
    done = (arg.running == 0);
    UNLOCK(arg.lock);
    if(arg.stop == 1 || done) {
#else
    struct timespec tmo = {time(NULL)+3,0};
    if(arg.stop == 1 || sem_timedwait(&arg.finished,&tmo) == 0) {
#endif
      for(thr = 0; thr < numthr; thr++) {
	(void)pthread_join(threads[thr], NULL);
      }
      break;
    }
#endif
  }
#endif

    /* print remaining messages */
  printmsg(arg.lock);
  /* Release synchronization gear */
#ifdef WIN
  CloseHandle(arg.lock);
#else
#ifdef HAVE_SEM
  sem_destroy(&arg.finished);
#endif
  pthread_mutex_destroy(arg.lock);
#endif


  if(arg.stop == 1) error("centering interrupted");
  if(quiet > 0 && arg.start != arg.last && K > 1)
    REprintf("...%d vectors centred in %d seconds\n",K,(int)(time(NULL)-arg.start));

  return(arg.badconv);
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

SEXP MY_demeanlist(SEXP vlist, SEXP flist, SEXP Ricpt, SEXP Reps,
		   SEXP scores, SEXP quiet, SEXP gkacc, SEXP Rmeans,
		   SEXP Rweights, SEXP Rscale) {
  int numvec;
  int numfac;
  int cnt;
  double **vectors;
  double **target;
  double *weights = NULL;
  double eps;
  SEXP reslist;
  FACTOR **factors;
  int cores;
  int i;
  int N=0;
  int icpt;
  int *vicpt;
  int icptlen;
  SEXP badconv;
  int listlen;
  int domeans=0;
  int protectcount=0;
  int scale[2];
  // Find the length of the data
  // We are never called with length(flist) == 0
  PROTECT(flist = AS_LIST(flist));  protectcount++;
  if(LENGTH(flist) == 0) {
    warning("demeanlist called with length(fl)==0, internal error?");
    N = 0;
  } else 
    N = LENGTH(VECTOR_ELT(flist,0));

  if(LENGTH(Rscale) == 1) {
    scale[0] = LOGICAL(Rscale)[0];
    scale[1] = scale[0];
  } else if(LENGTH(Rscale) > 1) {
    scale[0] = LOGICAL(Rscale)[0];
    scale[1] = LOGICAL(Rscale)[1];
  } else {
    error("scale must have length > 2");
  }

  domeans = LOGICAL(Rmeans)[0];
  vicpt = INTEGER(Ricpt);
  icptlen = LENGTH(Ricpt);
  icpt = vicpt[0] - 1; /* convert from 1-based to zero-based */
  eps = REAL(Reps)[0];
  cores = INTEGER(scores)[0];
  if(!isNull(Rweights)) {
    if(LENGTH(Rweights) != N) error("Length of weights (%d) must equal length of data (%d)",
				    LENGTH(Rweights),N);
    weights = REAL(PROTECT(coerceVector(Rweights, REALSXP)));
    protectcount++;
  }

  //  numfac = LENGTH(flist);
  factors = makefactors(flist, 1, weights);
  numfac = 0;
  for(FACTOR **f = factors; *f != NULL; f++) numfac++;

  int isdf = inherits(vlist,"data.frame");
  if(isdf) {
    icptlen = 1;
    icpt = -1;
  }

  PROTECT(vlist = AS_LIST(vlist)); protectcount++;
  listlen = LENGTH(vlist); 
  if(icptlen != 1 && icptlen != listlen)
    error("Length of icpt (%d) should be 1 or the number of arguments (%d)", icptlen, listlen);

  PROTECT(reslist = NEW_LIST(listlen)); protectcount++;
  SET_NAMES(reslist, GET_NAMES(vlist));
  if(isdf) {
    setAttrib(reslist, R_RowNamesSymbol, getAttrib(vlist, R_RowNamesSymbol));
    classgets(reslist,df_string);
    //    classgets(reslist,PROTECT(mkString("data.frame"))); UNPROTECT(1);
  }

  /* First, count the number of vectors in total */
  numvec = 0;
  for(int i = 0; i < listlen; i++) {
    SEXP elt = VECTOR_ELT(vlist,i);
    /* Each entry in the list is either a vector or a matrix */
    if(isNull(elt)) continue;
    if(!isMatrix(elt)) {
      if(LENGTH(elt) != N) 
	error("mtx[%d]: Factors and vectors must have the same length %d != %d",i, N,LENGTH(elt));
      numvec++;
    } else {
      if(nrows(elt) != N)
	error("mtx[%d]: Factors and vectors must have the same length %d != %d",i, N,nrows(elt));
      numvec += ncols(elt);
      if(icptlen != 1) icpt = vicpt[i]-1;
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
    if(isNull(elt)) continue;
    if(!isReal(elt)) {
      PROTECT(elt = coerceVector(elt, REALSXP)); protectcount++;
    }

    if(!isMatrix(elt)) {
      /* It's a vector */
      SEXP resvec;
      vectors[cnt] = REAL(elt);
      PROTECT(resvec = allocVector(REALSXP,LENGTH(elt)));
      SET_VECTOR_ELT(reslist,i,resvec);
      UNPROTECT(1);
      SET_NAMES(resvec, GET_NAMES(elt));
      target[cnt] = REAL(resvec);
      cnt++;
    } else {
      /* It's a matrix */
      int cols = ncols(elt);
      int rcols = cols;
      int j;
      SEXP mtx;
      /* Allocate a matrix */
      if(icptlen != 1) icpt = vicpt[i]-1;
      if(icpt >= 0 && icpt < cols) rcols--;
      PROTECT(mtx = allocMatrix(REALSXP,N,rcols));
      SET_VECTOR_ELT(reslist,i,mtx);
      UNPROTECT(1);
      // copy colnames 
      if(cols == rcols) {
	SET_DIMNAMES(mtx, GET_DIMNAMES(elt));
      }
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
  setAttrib(reslist,install("badconv"),badconv);
  UNPROTECT(1);
  INTEGER(badconv)[0] = demeanlist(vectors,N,numvec,target,weights,scale,
				   factors,numfac,
				   eps,cores,INTEGER(quiet)[0],
				   INTEGER(gkacc)[0]);

  if(INTEGER(badconv)[0] > 0)
    warning("%d vectors failed to centre to tolerance %.1le",INTEGER(badconv)[0],eps);


  if(domeans) {
    for(int i = 0; i < numvec; i++) {
      double *srcvec = vectors[i];
      double *dstvec = target[i];
      for(int j = 0; j < N; j++) {
	dstvec[j] = (srcvec[j]-dstvec[j]);
      }
    }
  }
  /* unprotect things */
  UNPROTECT(protectcount);
  return(reslist);
}

