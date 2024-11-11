
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void C_logLn_cal1_wei(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                             void *, void *, void *, void *);
extern void C_logLn_obs_pwc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                            void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void C_logL_cal_pwc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                           void *, void *, void *, void *, void *, void *, void *, void *);
extern void C_ppwc(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"C_logLn_cal1_wei", (DL_FUNC) &C_logLn_cal1_wei, 14},
  {"C_logLn_obs_pwc",  (DL_FUNC) &C_logLn_obs_pwc,  20},
  {"C_logL_cal_pwc",   (DL_FUNC) &C_logL_cal_pwc,   18},
  {"C_ppwc",           (DL_FUNC) &C_ppwc,            8},
  {NULL, NULL, 0}
};

void R_init_TwoPhaseScreenTest(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
