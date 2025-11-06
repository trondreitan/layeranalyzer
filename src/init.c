//#ifdef __linux__

//#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP layeranalyzer(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"layeranalyzer",                 (DL_FUNC) &layeranalyzer,                 25},
  {NULL, NULL, 0}
};

void R_init_layeranalyzer(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

// #endif // __linux__
