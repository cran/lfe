/*
 Author: Simen Gaure
 Copyright: 2011, Simen Gaure
 Licence: Artistic 2.0
*/
#include "lfe.h"

static R_CallMethodDef callMethods[] = {
  {"conncomp", (DL_FUNC) &R_conncomp, 1},
  {"wwcomp", (DL_FUNC) &R_wwcomp, 1},
  {"demeanlist", (DL_FUNC) &R_demeanlist, 8},
  {"kaczmarz", (DL_FUNC) &R_kaczmarz, 5},
  {"setdimnames", (DL_FUNC) &R_setdimnames, 2},
  {"scalecols", (DL_FUNC) &R_scalecols, 2},
  {"pdaxpy", (DL_FUNC) &R_pdaxpy, 3},
  {"sandwich", (DL_FUNC) &R_sandwich, 3},
  {"piproduct", (DL_FUNC) &R_piproduct, 2},
  {"dsyrk", (DL_FUNC) &R_dsyrk, 4},
  {"dsyr2k", (DL_FUNC) &R_dsyr2k, 5},
  {"address", (DL_FUNC) &R_address, 1},

  {NULL, NULL, 0}
};

void attribute_visible R_init_lfe(DllInfo *info) {
  /* register our routines */
  R_PreserveObject(df_string=mkString("data.frame"));
  (void)R_registerRoutines(info,NULL,callMethods,NULL,NULL);
}

void attribute_visible R_unload_lfe(DllInfo *info) {
  info=info; //avoid pedantic warning about unused parameter
  R_ReleaseObject(df_string);
}
