#ifndef _similarityMat_
#define _similarityMat_
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>

using namespace std;

enum similarity_type {w_jaccard};

class SIMMAT {
	public:
		SIMMAT();
		~SIMMAT();
		int getWeightedJaccard(const gsl_matrix*, CBLAS_UPLO_t, CBLAS_UPLO_t, gsl_matrix*, gsl_matrix*); // get both row similarity and column similarity
		int getWeightedJaccard(const gsl_matrix*, CBLAS_UPLO_t, gsl_matrix*, int); // 0: row; 1: col; Uplo specifies the triangle to store the symmetric matrix
		int getGraphLaplacian(const gsl_matrix*, gsl_matrix*, CBLAS_UPLO_t);
		int getGLFromPrior(const gsl_matrix*, CBLAS_UPLO_t, CBLAS_UPLO_t, gsl_matrix*, gsl_matrix*); // format the same as getWeightedJaccard
		int getGLFromPrior(const gsl_matrix*, CBLAS_UPLO_t, gsl_matrix*, int);
		double simMatSpectralNorm(const gsl_matrix*);

};
#endif
