#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void solveMyLasso(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void solveMyLassoF(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"solveMyLasso",  (DL_FUNC) &solveMyLasso,   9},
  {"solveMyLassoF", (DL_FUNC) &solveMyLassoF, 10},
  {NULL, NULL, 0}
};

void R_init_MGSDA(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}