#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void ks_matrix_R(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void matrix_density_R(void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"ks_matrix_R",      (DL_FUNC) &ks_matrix_R,      10},
    {"matrix_density_R", (DL_FUNC) &matrix_density_R,  7},
    {NULL, NULL, 0}
};

void R_init_psSubpathway(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
