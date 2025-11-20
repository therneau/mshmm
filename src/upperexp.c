/*
** Compute the eigenvectors for the upper triangular matrix R
**  If R = ADA^{-1} where D is a diagonal matrix of eigenvalues, then
**  S = exp(R) = A exp(D) A^{-1}
*/
#include <math.h>
#include "R.h"
#include "Rinternals.h"

SEXP upperexp(SEXP R2) {
    int i,j,k;
    int nc, ii;
    
    SEXP S2;
    double *R, *A, *S;
    double *dd, temp;

    nc = ncols(R2);   /* number of columns */
    R = REAL(R2);

    /* Make the output matrix a copy of R, so as to inherit
    **   the dimnames and etc
    */
    PROTECT(S2 = duplicate(R2));  /* exp(R)*/
    S = REAL(S2);
    dd = (double *) R_alloc(nc, sizeof(double));
    A =  (double *) R_alloc(nc*nc, sizeof(double)); 

    /*
    **	Compute the eigenvectors
    **   For each column of R, find x such that Rx = kx
    **   The eigenvalue k is R[i,i], x is a column of A
    **  Remember that R is in column order, so the i,j element is in
    **   location i + j*nc
    */
    ii =0; /* contains i * nc */
    for (i=0; i<nc; i++) { /* computations for column i */
	dd[i] = R[i +ii];    /* the i,i diagonal element */
	A[i +ii] = 1.0;
	for (j=(i-1); j >=0; j--) {  /* fill in the rest */
	    temp =0;
	    for (k=j+1; k<=i; k++) temp += R[j + k*nc]* A[k +ii];
	    A[j +ii] = temp/(dd[i]- R[j + j*nc]);
	    }
	ii += nc;
	}

    /* 
    ** The solution satisfies S= AD(A-inverse) where S is the solution and
    **  D is a diagonal matrix, or equivalently SA = AD.
    **  The elements of D are exp(diagonal of R).
    ** S will also be upper triangular, and we can solve for it using
    **  nearly the same code as above.  The prior block had Rz = kz with
    **  z an unknown column of A and k the eigenvalue.
    **  This has SA =x, with x a scaled column of A and S unknown.
    ** Imagine S and A are 4x4 and we are solving for the second row
    **  of S.  The equations are (remember that S[2,1]= A[2,3] = A[2,4] =0)
    **
    **    0*A[1,2] + S[2,2]A[2,2] + S[2,3] 0     + S[2,4] 0     = A[2,2] D[2]
    **    0*A[1,3] + S[2,2]A[2,3] + S[2,3]A[3,3] + S[2,4] 0     = A[2,3] D[3]
    **    0*A[1,4] + S[2,2]A[2,4] + S[2,3]A[3,4] + S[2,4]A[4,4] = A[2,4] D[4]
    **  
    */
    for (i=0; i<nc; i++) dd[i] = exp(dd[i]);  /* the scale */

    ii =0; /* contains i * nc */
    for (i=0; i<nc; i++) { /* computations for row i of S*/
	S[i +ii] = dd[i];
	for (j=i+1; j<nc; j++) {  /* fill in the rest of the row*/
	    temp =0;
	    for (k=i; k<j; k++) temp += S[i+ k*nc]* A[k + j*nc];
	    S[i + j*nc] = (A[i + j*nc]*dd[j] - temp)/A[j + j*nc];
	    }
	ii += nc;
	}

    unprotect(1);
    return(S2);
    }

	    
    
