/* C header for worst_VaR.c ***************************************************/

#ifndef worst_VaR_H
#define worst_VaR_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/* For RA and extensions specifically */
int num_opp_order_cols(double **X, int N, int d); /* count number of oppositely ordered cols in X */
void RA_aux(double **X, int N, int d, const char *method, const char *err,
	    int maxiter, double eps, int m_row_sums_size,
            double *individual_err, double *m_row_sums,
            int *num_opp_ordered, int *count); /* Steps 4 and 5 for the RA */
SEXP RA_aux_(SEXP X, SEXP method, SEXP err, SEXP maxiter, SEXP eps); /* interface to R */

#endif

