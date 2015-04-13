/*
 $Id: lfe.c 1671 2015-03-23 13:04:42Z sgaure $
*/
#include "lfe.h"

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

  {NULL, NULL, 0}
};

void attribute_visible R_init_lfe(DllInfo *info) {
  /* register our routines */
  (void) R_registerRoutines(info,NULL,callMethods,NULL,NULL);
  (void) R_PreserveObject(df_string=mkString("data.frame"));
}

void attribute_visible R_unload_lfe(DllInfo *info) {
  if(info != NULL){}; //avoid pedantic warning about unused parameter
  (void) R_ReleaseObject(df_string);
}
