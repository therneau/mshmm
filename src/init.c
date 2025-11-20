#include "R.h"
#include "R_ext/Rdynload.h"

/* Interface to expm package. */
typedef enum {Ward_2, Ward_1, Ward_buggy_octave} precond_type;
void (*expm)(double *x, int n, double *z, precond_type precond_kind);
void (*matexp_MH09)(double *x, int n, const int p, double *ret);
void R_init_hmm(DllInfo *dll)
{
    expm = (void (*)) R_GetCCallable("expm", "expm");
    matexp_MH09= (void (*)) R_GetCCallable("matexp_MH09", "matexp_MH09");
}

