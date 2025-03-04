#include <gsl/gsl_blas.h>
#include <algorithm>

#include "similarity_mat.h"

SIMMAT::SIMMAT() {};


SIMMAT::~SIMMAT() {};


int 
SIMMAT::getWeightedJaccard(const gsl_matrix* Confidence, CBLAS_UPLO_t row_uplo, CBLAS_UPLO_t col_uplo, gsl_matrix* row_s, gsl_matrix* col_s){
	int row = Confidence->size1;
	int col = Confidence->size2;

	// fill up row_s	
	for(int i = 0; i < row; ++i){
		gsl_matrix_set(row_s, i, i, 1);
		if(row_uplo == CblasUpper){
			for(int j = i+1; j < row; ++j){
				double min_sum = 0, max_sum = 0;
				for(int k = 0; k < col; ++k){
					double x = Confidence->data[i * Confidence->tda + k];
					double y = Confidence->data[j * Confidence->tda + k];
					min_sum += min(x,y);
					max_sum += max(x,y);
				}
				if(max_sum == 0){
					gsl_matrix_set(row_s, i, j, 0);
				} else {
					gsl_matrix_set(row_s, i, j, min_sum/max_sum);
				}
			}
		} else {
			for(int j = 0; j < i; ++j){
				double min_sum = 0, max_sum = 0;
				for(int k = 0; k < col; ++k){
					double x = Confidence->data[i * Confidence->tda + k];
					double y = Confidence->data[j * Confidence->tda + k];
					min_sum += min(x,y);
					max_sum += max(x,y);
				}
				if(max_sum == 0){
					gsl_matrix_set(row_s, i, j, 0);
				} else {
					gsl_matrix_set(row_s, i, j, min_sum/max_sum);
				}
			}
		}
	}

	// fill up col_s	
	for(int i = 0; i < col; ++i){
		gsl_matrix_set(col_s, i, i, 1);
		if(col_uplo == CblasUpper){
			for(int j = i+1; j < col; ++j){
				double min_sum = 0, max_sum = 0;
				for(int k = 0; k < row; ++k){
					double x = Confidence->data[k * Confidence->tda + i];
					double y = Confidence->data[k * Confidence->tda + j];
					min_sum += min(x,y);
					max_sum += max(x,y);
				}
				if(max_sum == 0){
					gsl_matrix_set(col_s, i, j, 0);
				} else {
					gsl_matrix_set(col_s, i, j, min_sum/max_sum);
				}
			}
		} else {
			for(int j = 0; j < i; ++j){
				double min_sum = 0, max_sum = 0;
				for(int k = 0; k < row; ++k){
					double x = Confidence->data[k * Confidence->tda + i];
					double y = Confidence->data[k * Confidence->tda + j];
					min_sum += min(x,y);
					max_sum += max(x,y);
				}
				if(max_sum == 0){
					gsl_matrix_set(col_s, i, j, 0);
				} else {
					gsl_matrix_set(col_s, i, j, min_sum/max_sum);
				}
			}
		}
	}
	return 0;
}

int 
SIMMAT::getWeightedJaccard(const gsl_matrix* Confidence, CBLAS_UPLO_t uplo, gsl_matrix* outsim, int row_0_col_1){
	int row = Confidence->size1;
	int col = Confidence->size2;

	if(row_0_col_1==0){
		// fill up row_s	
		for(int i = 0; i < row; ++i){
			gsl_matrix_set(outsim, i, i, 1);
			if(uplo == CblasUpper){
				for(int j = i+1; j < row; ++j){
					double min_sum = 0, max_sum = 0;
					for(int k = 0; k < col; ++k){
						double x = Confidence->data[i * Confidence->tda + k];
						double y = Confidence->data[j * Confidence->tda + k];
						min_sum += min(x,y);
						max_sum += max(x,y);
					}
					if(max_sum == 0){
						gsl_matrix_set(outsim, i, j, 0);
					} else {
						gsl_matrix_set(outsim, i, j, min_sum/max_sum);
					}
				}
			} else {
				for(int j = 0; j < i; ++j){
					double min_sum = 0, max_sum = 0;
					for(int k = 0; k < col; ++k){
						double x = Confidence->data[i * Confidence->tda + k];
						double y = Confidence->data[j * Confidence->tda + k];
						min_sum += min(x,y);
						max_sum += max(x,y);
					}
					if(max_sum == 0){
						gsl_matrix_set(outsim, i, j, 0);
					} else {
						gsl_matrix_set(outsim, i, j, min_sum/max_sum);
					}
				}
			}
		}
	} else if (row_0_col_1==1) {
		// fill up col_s	
		for(int i = 0; i < col; ++i){
			gsl_matrix_set(outsim, i, i, 1);
			if(uplo == CblasUpper){
				for(int j = i+1; j < col; ++j){
					double min_sum = 0, max_sum = 0;
					for(int k = 0; k < row; ++k){
						double x = Confidence->data[k * Confidence->tda + i];
						double y = Confidence->data[k * Confidence->tda + j];
						min_sum += min(x,y);
						max_sum += max(x,y);
					}
					if(max_sum == 0){
						gsl_matrix_set(outsim, i, j, 0);
					} else {
						gsl_matrix_set(outsim, i, j, min_sum/max_sum);
					}
				}
			} else {
				for(int j = 0; j < i; ++j){
					double min_sum = 0, max_sum = 0;
					for(int k = 0; k < row; ++k){
						double x = Confidence->data[k * Confidence->tda + i];
						double y = Confidence->data[k * Confidence->tda + j];
						min_sum += min(x,y);
						max_sum += max(x,y);
					}
					if(max_sum == 0){
						gsl_matrix_set(outsim, i, j, 0);
					} else {
						gsl_matrix_set(outsim, i, j, min_sum/max_sum);
					}
				}
			}
		}
	} else {
		printf("Please specify 0 (for similarity among rows) or 1 (cols)");
		return -1;
	}
	return 0;
}

int 
SIMMAT::getGraphLaplacian(const gsl_matrix* symMatrix, gsl_matrix* laplacian, CBLAS_UPLO_t uplo){
	int n = symMatrix->size1;
    gsl_matrix_memcpy(laplacian, symMatrix);
	gsl_matrix_scale(laplacian, -1);
	if (uplo == CblasUpper){
		for (int i = 0; i < n; ++i) {
			double degree = 0;
			for (int k = 0; k < i; ++k){
				degree += gsl_matrix_get(symMatrix, k, i);
			}
			for (int j = i+1; j < n; ++j) {
				degree += gsl_matrix_get(symMatrix, i, j);
			}
			gsl_matrix_set(laplacian, i, i, degree);
		}
	} else {
		for (int i = 0; i < n; ++i) {
			double degree = 0;
			for (int k = 0; k < i; ++k){
				degree += gsl_matrix_get(symMatrix, i, k);
			}
			for (int j = i+1; j < n; ++j) {
				degree += gsl_matrix_get(symMatrix, j, i);
			}
			gsl_matrix_set(laplacian, i, i, degree);
		}
	}
	return 0;
}

int 
SIMMAT::getGLFromPrior(const gsl_matrix* Confidence, CBLAS_UPLO_t row_uplo, CBLAS_UPLO_t col_uplo, gsl_matrix* row_s, gsl_matrix* col_s){
	getWeightedJaccard(Confidence, row_uplo, col_uplo, row_s, col_s);
	gsl_matrix* row_tmp = gsl_matrix_alloc(row_s->size1,row_s->size2);
	gsl_matrix* col_tmp = gsl_matrix_alloc(col_s->size1,col_s->size2);
	getGraphLaplacian(row_s, row_tmp, row_uplo);
	getGraphLaplacian(col_s, col_tmp, col_uplo);
	gsl_matrix_memcpy(row_s, row_tmp);
	gsl_matrix_memcpy(col_s, col_tmp);
	gsl_matrix_free(row_tmp);
	gsl_matrix_free(col_tmp);
	return 0;
}

int 
SIMMAT::getGLFromPrior(const gsl_matrix* Confidence, CBLAS_UPLO_t uplo, gsl_matrix* outsim, int row_0_col_1){
	getWeightedJaccard(Confidence, uplo, outsim, row_0_col_1);
	gsl_matrix* outsim_tmp = gsl_matrix_alloc(outsim->size1,outsim->size2);
	getGraphLaplacian(outsim, outsim_tmp, uplo);
	gsl_matrix_memcpy(outsim, outsim_tmp);
	gsl_matrix_free(outsim_tmp);
	return 0;
}

// mat should follow CblasLower
double
SIMMAT::simMatSpectralNorm(const gsl_matrix* mat){
	int n = mat->size1;

	gsl_matrix* matrix = gsl_matrix_alloc(n, n);
	gsl_vector* eigenvalues = gsl_vector_alloc(n);
	gsl_matrix_memcpy(matrix, mat);
    gsl_eigen_symm_workspace* w = gsl_eigen_symm_alloc(n);
	gsl_vector_view diag = gsl_matrix_diagonal(matrix); // preconditioning. The QR step just halts.
    gsl_vector_add_constant(&diag.vector, 1e-4);
    gsl_eigen_symm(matrix, eigenvalues, w);

    double ret = gsl_vector_max(eigenvalues);
    
    gsl_eigen_symm_free(w);
    gsl_matrix_free(matrix);
    gsl_vector_free(eigenvalues);
    return ret;
}