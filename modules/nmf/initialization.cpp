#include <iostream>
#include <fstream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <list>
#include <math.h>
#include <queue>
#include <string>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include "initialization.h"
#include "nmf_wo_init.h"


extern "C" {
	#include "random_svd/low_rank_svd_algorithms_gsl.h"
}
// A ~ W^TH
int init::nndsvd(gsl_matrix* A, gsl_matrix* W, gsl_matrix* H, int seed) {
	int m = A->size1;
	int n = A->size2;
	int k = W->size1;

	double a = 0; //avg 
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			a += A->data[i * A->tda + j];	
		}
	}
	a /= (m*n);
	
	
	gsl_matrix *U, *S, *V;
	
	randomized_low_rank_svd3(A, k, 3, 1, &U, &S, &V, seed);
	
	gsl_vector_view x = gsl_matrix_column(U, 0);
	gsl_vector_view y = gsl_matrix_column(V, 0);
	double sigma = sqrt(gsl_matrix_get(S, 0, 0));
	gsl_vector_scale(&x.vector, sigma);
	gsl_vector_scale(&y.vector, sigma);
	gsl_matrix_set_row(W, 0, &x.vector);
	gsl_matrix_set_row(H, 0, &y.vector);

	for (int i = 1; i < k; i++) {
		gsl_vector* xp = gsl_vector_calloc(m);
		gsl_vector* xn = gsl_vector_calloc(m);
		gsl_vector* yp = gsl_vector_calloc(n);
		gsl_vector* yn = gsl_vector_calloc(n);

		for (int j = 0; j < m; j++) {
			double val = gsl_matrix_get(U, j, i);
			if (val >= 0) {
				gsl_vector_set(xp, j, val);
			} else {	
				gsl_vector_set(xn, j, val*-1);
			}
		}		

		for (int j = 0; j < n; j++) {
			double val = gsl_matrix_get(V, j, i);
			if (val >= 0) {
				gsl_vector_set(yp, j, val);
			} else {	
				gsl_vector_set(yn, j, val*-1);
			}
		}		
	
		double xpNorm = gsl_blas_dnrm2(xp);
		double xnNorm = gsl_blas_dnrm2(xn);
		double ypNorm = gsl_blas_dnrm2(yp);
		double ynNorm = gsl_blas_dnrm2(yn);

		double mp = xpNorm * ypNorm;
		double mn = xnNorm * ynNorm;


		if (mp > mn) {
			gsl_vector_scale(xp, 1/xpNorm * sqrt(gsl_matrix_get(S,i,i)*mp));
			gsl_vector_scale(yp, 1/ypNorm * sqrt(gsl_matrix_get(S,i,i)*mp));	
			gsl_matrix_set_row(W, i, xp);
			gsl_matrix_set_row(H, i, yp);	
		} else {
			gsl_vector_scale(xn, 1/xnNorm * sqrt(gsl_matrix_get(S,i,i)*mn));
			gsl_vector_scale(yn, 1/ynNorm * sqrt(gsl_matrix_get(S,i,i)*mn));
			gsl_matrix_set_row(W, i, xn);
			gsl_matrix_set_row(H, i, yn);
		}

		gsl_vector_free(xp);
		gsl_vector_free(xn);
		gsl_vector_free(yp);
		gsl_vector_free(yn);
	}
	
	gsl_matrix_free(U);
	gsl_matrix_free(V);
	gsl_matrix_free(S);	
	
	for (int i=0; i < k; i++) {
		for (int j=0; j < m; j++) {
			double *val = &(W->data[i * W->tda + j]);
			if (*val <= 0) {
				*val = a;
			}			
		}
	}

	for (int i=0; i < k; i++) {
		for (int j=0; j < n; j++) {
			double *val = &(H->data[i * H->tda + j]);
			if (*val <= 0) {
				*val = a;
			}			
		}
	}

	return 0;
}
// Allow negative entries
int init::nndsvd_general(gsl_matrix* A, gsl_matrix* W, gsl_matrix* H, int seed) {
	int m = A->size1;
	int n = A->size2;
	int k = W->size1;

	double a = 0; //avg 
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			a += A->data[i * A->tda + j];	
		}
	}
	a /= (m*n);
	
	
	gsl_matrix *U, *S, *V;
	
	randomized_low_rank_svd3(A, k, 3, 1, &U, &S, &V, seed);
	
	gsl_vector_view x = gsl_matrix_column(U, 0);
	gsl_vector_view y = gsl_matrix_column(V, 0);
	double sigma = sqrt(gsl_matrix_get(S, 0, 0));
	gsl_vector_scale(&x.vector, sigma);
	gsl_vector_scale(&y.vector, sigma);
	gsl_matrix_set_row(W, 0, &x.vector);
	gsl_matrix_set_row(H, 0, &y.vector);

	for (int i = 1; i < k; i++) {
		gsl_vector* xp = gsl_vector_calloc(m);
		gsl_vector* xn = gsl_vector_calloc(m);
		gsl_vector* yp = gsl_vector_calloc(n);
		gsl_vector* yn = gsl_vector_calloc(n);

		for (int j = 0; j < m; j++) {
			double val = gsl_matrix_get(U, j, i);
			if (val >= 0) {
				gsl_vector_set(xp, j, val);
			} else {	
				gsl_vector_set(xn, j, val*-1);
			}
		}		

		for (int j = 0; j < n; j++) {
			double val = gsl_matrix_get(V, j, i);
			if (val >= 0) {
				gsl_vector_set(yp, j, val);
			} else {	
				gsl_vector_set(yn, j, val*-1);
			}
		}		
	
		double xpNorm = gsl_blas_dnrm2(xp);
		double xnNorm = gsl_blas_dnrm2(xn);
		double ypNorm = gsl_blas_dnrm2(yp);
		double ynNorm = gsl_blas_dnrm2(yn);

		double mp = xpNorm * ypNorm;
		double mn = xnNorm * ynNorm;


		if (mp > mn) {
			gsl_vector_scale(xp, 1/xpNorm * sqrt(gsl_matrix_get(S,i,i)*mp));
			gsl_vector_scale(yp, 1/ypNorm * sqrt(gsl_matrix_get(S,i,i)*mp));	
			gsl_matrix_set_row(W, i, xp);
			gsl_matrix_set_row(H, i, yp);	
		} else {
			gsl_vector_scale(xn, 1/xnNorm * sqrt(gsl_matrix_get(S,i,i)*mn));
			gsl_vector_scale(yn, 1/ynNorm * sqrt(gsl_matrix_get(S,i,i)*mn));
			gsl_matrix_set_row(W, i, xn);
			gsl_matrix_set_row(H, i, yn);
		}

		gsl_vector_free(xp);
		gsl_vector_free(xn);
		gsl_vector_free(yp);
		gsl_vector_free(yn);
	}
	
	gsl_matrix_free(U);
	gsl_matrix_free(V);
	gsl_matrix_free(S);	

	return 0;
}


// A ~ WH^T
int init::nndsvd2(gsl_matrix* A, gsl_matrix* W, gsl_matrix* H, int seed) {
	int m = A->size1;
	int n = A->size2;
	int k = W->size2;

	// double a = 0; //avg 
	// for (int i = 0; i < m; i++) {
	// 	for (int j = 0; j < n; j++) {
	// 		a += A->data[i * A->tda + j];	
	// 	}
	// }
	// a /= (m*n);
	
	
	gsl_matrix *U, *S, *V;
	
	randomized_low_rank_svd3(A, k, 3, 1, &U, &S, &V, seed);
	
	gsl_vector_view x = gsl_matrix_column(U, 0);
	gsl_vector_view y = gsl_matrix_column(V, 0);
	double sigma = sqrt(gsl_matrix_get(S, 0, 0));
	gsl_vector_scale(&x.vector, sigma);
	gsl_vector_scale(&y.vector, sigma);
	gsl_matrix_set_col(W, 0, &x.vector);
	gsl_matrix_set_col(H, 0, &y.vector);

	for (int i = 1; i < k; i++) {
		gsl_vector* xp = gsl_vector_calloc(m);
		gsl_vector* xn = gsl_vector_calloc(m);
		gsl_vector* yp = gsl_vector_calloc(n);
		gsl_vector* yn = gsl_vector_calloc(n);

		for (int j = 0; j < m; j++) {
			double val = gsl_matrix_get(U, j, i);
			if (val >= 0) {
				gsl_vector_set(xp, j, val);
			} else {	
				gsl_vector_set(xn, j, val*-1);
			}
		}		

		for (int j = 0; j < n; j++) {
			double val = gsl_matrix_get(V, j, i);
			if (val >= 0) {
				gsl_vector_set(yp, j, val);
			} else {	
				gsl_vector_set(yn, j, val*-1);
			}
		}		
	
		double xpNorm = gsl_blas_dnrm2(xp);
		double xnNorm = gsl_blas_dnrm2(xn);
		double ypNorm = gsl_blas_dnrm2(yp);
		double ynNorm = gsl_blas_dnrm2(yn);

		double mp = xpNorm * ypNorm;
		double mn = xnNorm * ynNorm;


		if (mp > mn) {
			gsl_vector_scale(xp, 1/xpNorm * sqrt(gsl_matrix_get(S,i,i)*mp));
			gsl_vector_scale(yp, 1/ypNorm * sqrt(gsl_matrix_get(S,i,i)*mp));	
			gsl_matrix_set_col(W, i, xp);
			gsl_matrix_set_col(H, i, yp);	
		} else {
			gsl_vector_scale(xn, 1/xnNorm * sqrt(gsl_matrix_get(S,i,i)*mn));
			gsl_vector_scale(yn, 1/ynNorm * sqrt(gsl_matrix_get(S,i,i)*mn));
			gsl_matrix_set_col(W, i, xn);
			gsl_matrix_set_col(H, i, yn);
		}

		gsl_vector_free(xp);
		gsl_vector_free(xn);
		gsl_vector_free(yp);
		gsl_vector_free(yn);
	}
	
	gsl_matrix_free(U);
	gsl_matrix_free(V);
	gsl_matrix_free(S);	
	
	// for (int i=0; i < k; i++) {
	// 	for (int j=0; j < m; j++) {
	// 		double *val = &(W->data[i * W->tda + j]);
	// 		if (*val <= 0) {
	// 			*val = a;
	// 		}			
	// 	}
	// }

	// for (int i=0; i < k; i++) {
	// 	for (int j=0; j < n; j++) {
	// 		double *val = &(H->data[i * H->tda + j]);
	// 		if (*val <= 0) {
	// 			*val = a;
	// 		}			
	// 	}
	// }

	return 0;
}

// A ~ UV^TWH^T 
// don't use. Not stable
int init::nndsvd(gsl_matrix* A, gsl_matrix* U, gsl_matrix* V, gsl_matrix* W, gsl_matrix* H, int seed) {
	int k1 = U->size2;
	int p = W->size1;
	int k2 = W->size2;
	int m = A->size1;
	int n = A->size2;

	gsl_matrix* X = gsl_matrix_alloc(p, m); // VU^T
	gsl_matrix* Y = gsl_matrix_alloc(p, n); // WH^T
	double a = 0; //avg 
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			a += A->data[i * A->tda + j];	
		}
	}
	a /= (m*n);

	nndsvd(A, X, Y, seed); // A = X^TY
	NMFWI nmf_initializer(p, 5, seed, false, 1);
	nmf_initializer.fit(A, X, Y);

	nndsvd2(X, V, U, seed);
	nndsvd2(Y, W, H, seed);
	NMFWI nmf_X_initializer(k1, 300, seed, false, 0.1);
	NMFWI nmf_Y_initializer(k2, 300, seed, false, 0.1);

	gsl_matrix* UT = gsl_matrix_alloc(k1, m);
	gsl_matrix* VT = gsl_matrix_alloc(k1, p);
	gsl_matrix_transpose_memcpy(UT, U);
	gsl_matrix_transpose_memcpy(VT, V);
	nmf_X_initializer.fit(X, VT, UT);
	gsl_matrix_transpose_memcpy(U, UT);
	gsl_matrix_transpose_memcpy(V, VT);

	gsl_matrix* WT = gsl_matrix_alloc(k2, p);
	gsl_matrix* HT = gsl_matrix_alloc(k2, n);
	gsl_matrix_transpose_memcpy(WT, W);
	gsl_matrix_transpose_memcpy(HT, H);
	nmf_Y_initializer.fit(Y, WT, HT);
	gsl_matrix_transpose_memcpy(W, WT);
	gsl_matrix_transpose_memcpy(H, HT);

	gsl_matrix_free(X);
	gsl_matrix_free(Y);
	gsl_matrix_free(UT);
	gsl_matrix_free(VT);
	gsl_matrix_free(WT);
	gsl_matrix_free(HT);

	gsl_matrix* u_mat[4];
	u_mat[0] = U;
	u_mat[1] = V;
	u_mat[2] = W;
	u_mat[3] = H;
	
	// the below causes the initial loss to spike. Don't use.
	// // avoid divide by zero
	// for(int k = 0; k < 4; ++k){
	// 	gsl_matrix* cur_mat = u_mat[k];
	// 	for (int i=0; i < cur_mat->size1; i++) {
	// 		for (int j=0; j < cur_mat->size2; j++) {
	// 			double *val = &(cur_mat->data[i * cur_mat->tda + j]);
	// 			if (*val <= 0) {
	// 				*val = a;
	// 			}			
	// 		}
	// 	}
	// }

	// for (int k = 0; k < k1; ++k) {
	// 	gsl_vector_view u_k = gsl_matrix_column(U, k);
	// 	gsl_vector_view v_k = gsl_matrix_column(V, k);
	// 	double norm = gsl_blas_dnrm2(&u_k.vector);
	// 	gsl_vector_scale(&u_k.vector, 1/norm);
	// 	gsl_vector_scale(&v_k.vector, norm);
	// }
	// for (int k = 0; k < k2; ++k) {
	// 	gsl_vector_view w_k = gsl_matrix_column(W, k);
	// 	gsl_vector_view h_k = gsl_matrix_column(H, k);
	// 	double norm = gsl_blas_dnrm2(&h_k.vector);
	// 	gsl_vector_scale(&h_k.vector, 1/norm);
	// 	gsl_vector_scale(&w_k.vector, norm);
	// }
	return 0;
}

int init::random(gsl_matrix* W, gsl_matrix* H, int seed) {
	const gsl_rng_type* T;
	gsl_rng* ri;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	ri = gsl_rng_alloc(T);
	gsl_rng_set(ri, seed);
	int row = W->size1;
	int col = W->size2;
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			double *val = &(W->data[i * W->tda + j]);
			*val = gsl_rng_uniform(ri);
		}
	}
	row = H->size1;
	col = H->size2;
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			double *val = &(H->data[i * H->tda + j]);
			*val = gsl_rng_uniform(ri);
		}
	}

	gsl_rng_free(ri);
	return 0;
}

int init::random(vector<gsl_matrix*>& U, int seed) {
	const gsl_rng_type* T;
	gsl_rng* ri;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	ri = gsl_rng_alloc(T);
	gsl_rng_set(ri, seed);
	for (vector<gsl_matrix*>::iterator it = U.begin(); it != U.end(); ++it) {
    	gsl_matrix* W = *it;
		int row = W->size1;
		int col = W->size2;
		for (int i = 0; i < row; ++i) {
			for (int j = 0; j < col; ++j) {
				double *val = &(W->data[i * W->tda + j]);
				*val = gsl_rng_uniform(ri);
			}
		}
	}
	gsl_rng_free(ri);
	return 0;
}



