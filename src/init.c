/* #ifdef __linux__ // not necessary? */

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP layeranalyzer(SEXP /*input*/,
			  SEXP /*num_MCMC*/,
			  SEXP /*Burnin*/,
			  SEXP /*Spacing*/,
			  SEXP /*NumTemp*/,
			  SEXP /*do_model_likelihood*/,
			  SEXP /*do_maximum_likelihood*/,
			  SEXP /*maximum_likelihood_numstart*/,
			  SEXP /*SilentMode*/,
			  SEXP /*TalkativeBurnin*/,
			  SEXP /*TalkativeLikelihood*/,
			  SEXP /*IdStrategy*/,
			  SEXP /*UseStationarySdev*/,
			  SEXP /*TempGround*/,
			  SEXP /*UseHalfLives*/,
			  SEXP /*ReturnMCMC*/, 
			  SEXP /*causal*/,
			  SEXP /*causal_symmetric*/,
			  SEXP /*corr*/,
			  SEXP /*smooth_specs*/,
			  SEXP /*realization_specs*/,
			  SEXP /*ReturnResiduals*/,
			  SEXP /*mode*/,
			  SEXP /*loglik_laxness*/,
			  SEXP /*input_param_values*/,
			  SEXP /*input_mcmc*/,
			  SEXP /*external_input*/,
			  SEXP /*distance_matrix*/);


static const R_CMethodDef CEntries[] = {
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"layeranalyzer",                 (DL_FUNC) &layeranalyzer,28},
  {NULL, NULL, 0}
};

void R_init_layeranalyzer(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

// #endif // __linux__
