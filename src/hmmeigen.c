/*
** This is a modification of the routine La_rg, from src/modules/lapack/Lapack.c
**  in the R source code.  The primary difference is that we return both the
**  right and left eigenvectors, as well as the eigenvalues.
*/
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <stdbool.h>

static SEXP unscramble(const double* imaginary, int n,
		       const double* vecs)
{
    int i, j;
    SEXP s = allocMatrix(CPLXSXP, n, n);
    size_t N = n;
    for (j = 0; j < n; j++) {
	if (imaginary[j] != 0) {
	    int j1 = j + 1;
	    for (i = 0; i < n; i++) {
		COMPLEX(s)[i+N*j].r = COMPLEX(s)[i+N*j1].r = vecs[i + j * N];
		COMPLEX(s)[i+N*j1].i = -(COMPLEX(s)[i+N*j].i = vecs[i + j1 * N]);
	    }
	    j = j1;
	} else {
	    for (i = 0; i < n; i++) {
		COMPLEX(s)[i+N*j].r = vecs[i + j * N];
		COMPLEX(s)[i+N*j].i = 0.0;
	    }
	}
    }
    return s;
}

SEXP hmmeigen(SEXP x)
{
    bool complexValues;
    int i, n, lwork, info, *xdims;
    double *work, *wR, *wI, *left, *right, *xvals, tmp;
    char jobVL[2] = "V", jobVR[2] = "V";
    static const char *outnames[]= {"values", "left", "right", ""};
    SEXP rval, val;

    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    n = xdims[0];
    if (n != xdims[1])
	Rf_error("'x' must be a square numeric matrix");

    /* work on a copy of x */
    if (!isReal(x)) {
	x = coerceVector(x, REALSXP);
	xvals = REAL(x);
    } else {
	xvals = (double *) R_alloc(n * (size_t)n, sizeof(double));
	Memcpy(xvals, REAL(x), (size_t) n * n);
    }
    PROTECT(x);

    /* space for left and right eigenvectors */
    right = (double *) R_alloc(n * (size_t)n, sizeof(double));
    left  = (double *) R_alloc(n * (size_t)n, sizeof(double));

    wR = (double *) R_alloc(n, sizeof(double));
    wI = (double *) R_alloc(n, sizeof(double));
    /* ask for optimal size of work array */
    lwork = -1;
    F77_CALL(dgeev)(jobVL, jobVR, &n, xvals, &n, wR, wI,
		    left, &n, right, &n, &tmp, &lwork, &info FCONE FCONE);
    if (info != 0)
	Rf_error("error code %d from Lapack routine dgeev", info);
    lwork = (int) tmp;
    work = (double *) R_alloc(lwork, sizeof(double));
    F77_CALL(dgeev)(jobVL, jobVR, &n, xvals, &n, wR, wI,
		    left, &n, right, &n, work, &lwork, &info FCONE FCONE);
    if (info != 0)
	Rf_error("error code %d from Lapack routine dgeev", info);

    complexValues = false;
    for (i = 0; i < n; i++)
	/* This test used to be !=0 for R < 2.3.0.  This is OK for 0+0i */
	/* tmt: I can't find R_AccuracyInfo in any .h files, sub in a value */
	/*if (fabs(wI[i]) >  10 * R_AccuracyInfo.eps * fabs(wR[i])) { */
	if (fabs(wI[i]) >  2.2e-15 * fabs(wR[i])) {
	    complexValues = true;
	    break;
	}

    PROTECT(rval = mkNamed(VECSXP, outnames));
    if (complexValues) {
	val = SET_VECTOR_ELT(rval, 0, allocVector(CPLXSXP, n));
	for (i = 0; i < n; i++) {
	    COMPLEX(val)[i].r = wR[i];
	    COMPLEX(val)[i].i = wI[i];
	}
	SET_VECTOR_ELT(rval, 1, unscramble(wI, n, left));
	SET_VECTOR_ELT(rval, 2, unscramble(wI, n, right));
    } else {
	val = SET_VECTOR_ELT(rval, 0, allocVector(REALSXP, n));
	for (i = 0; i < n; i++) REAL(val)[i] = wR[i];

 	val = SET_VECTOR_ELT(rval, 1, allocMatrix(REALSXP, n, n));
	for (i = 0; i < (n * n); i++) REAL(val)[i] = left[i];

	val = SET_VECTOR_ELT(rval, 2, allocMatrix(REALSXP, n, n));
	for (i = 0; i < (n * n); i++) REAL(val)[i] = right[i];
    }
    UNPROTECT(2);
    return(rval);
}
