/*
 Author: Simen Gaure
 Copyright: 2011, Simen Gaure
 Licence: Artistic 2.0
*/
#include "lfe.h"

SEXP R_address(SEXP x);
static R_CallMethodDef callMethods[] = {
  {"conncomp", (DL_FUNC) &R_conncomp, 1},
  {"wwcomp", (DL_FUNC) &R_wwcomp, 1},
  {"demeanlist", (DL_FUNC) &R_demeanlist, 8},
  {"kaczmarz", (DL_FUNC) &R_kaczmarz, 5},
  {"setdimnames", (DL_FUNC) &R_setdimnames, 2},
  {"scalecols", (DL_FUNC) &R_scalecols, 2},
  {"pdaxpy", (DL_FUNC) &R_pdaxpy, 3},
  {"piproduct", (DL_FUNC) &R_piproduct, 2},
  {"dsyrk", (DL_FUNC) &R_dsyrk, 4},
  {"dsyr2k", (DL_FUNC) &R_dsyr2k, 5},
  {"address", (DL_FUNC) &R_address, 1},

  {NULL, NULL, 0}
};
static R_ExternalMethodDef externalMethods[] = {
  {"edemeanlist", (DL_FUNC) &RE_demeanlist, -1},
  {NULL, NULL, 0}
};

void attribute_visible R_init_lfe(DllInfo *info) {
  str_DF = mkString("data.frame");
  /* register our routines */
  (void)R_registerRoutines(info,NULL,callMethods,NULL,externalMethods);
}

