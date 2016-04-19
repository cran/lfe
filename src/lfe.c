/*
 $Id: lfe.c 1991 2016-04-14 15:56:45Z sgaure $
*/
#include "lfe.h"
SEXP MY_threads(SEXP rt) {
  if(LENGTH(rt) < 1) return R_NilValue;
  LFE_GLOBAL_THREADS = INTEGER(rt)[0];
  return R_NilValue;
}
static R_CallMethodDef callMethods[] = {
  {"conncomp", (DL_FUNC) &MY_conncomp, 1},
  {"wwcomp", (DL_FUNC) &MY_wwcomp, 1},
  {"demeanlist", (DL_FUNC) &MY_demeanlist, 10},
  {"kaczmarz", (DL_FUNC) &MY_kaczmarz, 5},
  {"setdimnames", (DL_FUNC) &MY_setdimnames, 2},
  {"scalecols", (DL_FUNC) &MY_scalecols, 2},
  {"pdaxpy", (DL_FUNC) &MY_pdaxpy, 3},
  {"sandwich", (DL_FUNC) &MY_sandwich, 3},
  {"piproduct", (DL_FUNC) &MY_piproduct, 2},
  {"dsyrk", (DL_FUNC) &MY_dsyrk, 4},
  {"dsyr2k", (DL_FUNC) &MY_dsyr2k, 5},
  {"address", (DL_FUNC) &MY_address, 1},
  {"named", (DL_FUNC) &MY_named, 2},
  {"rowsum", (DL_FUNC) &Crowsum, 3},
  //  {"setnamed", (DL_FUNC) &MY_setnamed, 2},
  //  {"ppf", (DL_FUNC) &MY_ppf, 2},
  //  {"threads", (DL_FUNC) &MY_threads, 1},
  {NULL, NULL, 0}
};

void attribute_visible R_init_lfe(DllInfo *info) {
  /* register our routines */
  (void) R_registerRoutines(info,NULL,callMethods,NULL,NULL);
  (void) R_PreserveObject(df_string=mkString("data.frame"));
  LFE_GLOBAL_THREADS=1;
}

void attribute_visible R_unload_lfe(DllInfo *info) {
  if(info != NULL){}; //avoid pedantic warning about unused parameter
  (void) R_ReleaseObject(df_string);
}
