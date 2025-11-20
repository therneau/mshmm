/*
** Compute the matrix exponential of a set of upper triangular matrices
**  nstate = integer number of states
**  eta    = a matrix with n rows (one per subject) and one column for
**    each linear predictor
**  map  = a vector with length = number of linear predictors, containing the
**     element in the rate matrix it plugs into
**  time = compute the matrix exponential exp(Rt), each row has it's own value
**     of t.
**  nterm = the order of approx from 6-13 for Higham09, if <6 the Ward method
**    is used
**  return value an array P(nstate, nstate, n) containing the matrix
**    exponentials
*/
#include <math.h>
#include "R.h"
#include "Rinternals.h"

typedef enum {Ward_2, Ward_1, Ward_buggy_octave} precond_type;
extern void (*expm)(double *x, int n, double *z, precond_type precond_kind);
extern void (*matexp_MH09)(double *x, int n, const int p, double *ret);

SEXP upper(SEXP nstate2, SEXP eta2, SEXP time2, SEXP map2, SEXP eps2,
	   SEXP nterm2) {
    int nterm =0;  
    int i, ii, j,k, irow;
    int nc, n, npar;
    
    static const char *outnames[] = {"P", "ties", ""};

    SEXP P2, rval;
    double *P, *time;
    int *map, *tied_eigen;
    double *eta;
    int n2;    /*n2 is how many we are doing */
    double *A, *R, *d;
    double temp;
    double eps;

    nterm  = asInteger(nterm2);
    nc     = asInteger(nstate2);  /* number of rows/cols in rate matrix R */
    eta    = REAL(eta2);
    map    = INTEGER(map2);
    time   = REAL(time2);
    eps    = REAL(eps2)[0];

    npar = LENGTH(map2);   /* number of columns */
    n    = nrows(eta2);   /* rows in the eta matrix */
    n2   = LENGTH(time2);

    PROTECT(rval = mkNamed(VECSXP, outnames));
    P2 = SET_VECTOR_ELT(rval, 0, allocVector(REALSXP, nc*nc*n2));
    P  = REAL(P2);
    tied_eigen = INTEGER(SET_VECTOR_ELT(rval, 1, allocVector(INTSXP, 1)));
    /*
    ** scratch space
    */
    A = (double *) R_alloc(nc*nc, sizeof(double));
    R = (double *) R_alloc(nc*nc, sizeof(double));
    d = (double *) R_alloc(nc, sizeof(double));
 
    for (i=0; i<nc*nc; i++) {  /* R doesn't zero memory */
	R[i] =0.0;  
    }
    for (i=0; i<nc*nc*n2; i++) P[i] = 0.0;

    tied_eigen[0] = 0;  /* for my curiosity -- how often does this happen? */
    for (irow=0; irow<n2; irow++) {
	/* populate the R matrix */
	for (i=0; i<npar; i++) R[map[i]-1] = exp(eta[irow + i*n]);
	/* 
	**  fill in the diagonal of R such that row sums are 0 
	**  since this is upper triangular and row sums have to be zero,
	**   the lower right corner element must be zero
	*/
	k =0;
	for (i =0; i<(nc-1); i++) {
	    temp =0;
	    for (j=i+1; j<nc; j++) temp += R[i + j*nc];
	    R[i + i*nc] = -temp;
	    d[i] = temp;
	    if (temp== 0.){ /* backsolve algorithm fails with 2 zeros though */
		for (j=0; j<i; j++) if (d[j] ==0) k=1;
	    } else{	
		for (j=0; j<i; j++) /* use the same rule as the R code */
		    if (fabs(temp- d[j]) < eps) k=1;
	    }
	}	

	/* compute exp(R).  
	** If there are tied eigenvalued default to the expm package
	*/
	if (k>0) {
	    tied_eigen[0]++;
	    for (i=0; i<(nc*nc); i++) R[i] *= time[irow];
	    if (nterm <6) expm(R, nc, P, Ward_2) ;
	    else matexp_MH09(R, nc, nterm, P); 
	    }
	else {
	    /* 
	    **   For each column of R, find x such that Rx = kx
	    **   The eigenvalue k is R[i,i], A contains the eigenvectors x
	    **  Remember that R is in column order, so the i,j element is in
	    **   location i + j*nc
	    */
	    ii =0; /* contains i * nc */
	    for (i=0; i<nc; i++) { /* computations for column i */
		d[i] = R[i +ii];    /* the i,i diagonal element = eigenvalue*/
		A[i +ii] = 1.0;
		for (j=(i-1); j >=0; j--) {  /* fill in the rest */
		    temp =0;
		    for (k=j+1; k<=i; k++) temp += R[j + k*nc]* A[k +ii];
		    A[j +ii] = temp/(d[i]- R[j + j*nc]);
		}
		ii += nc;
	    }
	    
	    /* Now compute exp(R) = A exp(D) A-inverse 
	    **  we actually never need an explicit A-inverse, since A is
	    **  triangluar
	    */
	    ii =0; /* contains i * nc */
	    for (i=0; i<nc; i++) d[i] = exp(time[irow]* d[i]);
	    for (i=0; i<nc; i++) {
		P[i +ii] = d[i];
		for (j=i+1; j<nc; j++) { 
		    temp =0;
		    for (k=i; k<j; k++) temp += P[i+ k*nc]* A[k + j*nc];
		    P[i + j*nc] = (A[i + j*nc]*d[j] - temp)/A[j + j*nc];
		} 
	    ii += nc;
	    }
	}
	P += nc*nc;
    }
    UNPROTECT(1);
    return(rval);
}
