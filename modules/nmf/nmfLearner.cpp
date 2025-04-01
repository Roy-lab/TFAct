#include <pthread.h>
#include <semaphore.h>
#include <gsl/gsl_blas.h>

#include "nmfLearner.h"
#include "utils.h"

NMFLearner::NMFLearner(int k, int UV_dim, int WH_dim, init_method initMethod, 
int maxIter, int k1_maxIter, int k2_maxIter, int seed, bool verb, double termTol, bool perm_up,
nmfType nt) 
{
	n_components = k;
	init = initMethod;
	max_iter = maxIter;
	random_state = seed;
	verbose = verb;
	tol = termTol;
	numNcaRep = -1;
	permuted_update = perm_up;
	k1 = UV_dim;
	k2 = WH_dim;
	UV_maxIter = k1_maxIter;
	WH_maxIter = k2_maxIter;
	nmft = nt;
	lambda_U = 0;
	lambda_V = 0;
}

// 2 factors
NMFLearner::NMFLearner(int k, init_method initMethod, 
int maxIter, int seed, bool verb, double termTol,
nmfType nt) 
{
	n_components = k;
	init = initMethod;
	max_iter = maxIter;
	random_state = seed;
	verbose = verb;
	tol = termTol;
	numNcaRep = -1;
	nmft = nt;
}

NMFLearner::~NMFLearner() {
	if (numNcaRep > 0){
		if ((nmft == recur) || (nmft == nmqf)|| (nmft == ibpg)){
			for (int i = 0; i < numNcaRep; ++i) {
				gsl_matrix_free(allA[i]);
				gsl_matrix_free(allP[i]);
				gsl_matrix_free(allU[i]);
				gsl_matrix_free(allV[i]);
				gsl_matrix_free(allW[i]);
				gsl_matrix_free(allH[i]);
			}
			delete[] allA;
			delete[] allP;
			delete[] allU;
			delete[] allV;
			delete[] allW;
			delete[] allH;
			for (int i = 0; i < numNcaRep; ++i) {
				delete nmfArray[i];
			}
			delete[] nmfArray;
		} else if ((nmft == hals_prior) || (nmft == ibpg_2f)){
			for (int i = 0; i < numNcaRep; ++i) {
				gsl_matrix_free(allA[i]);
				gsl_matrix_free(allP[i]);
			}
			delete[] allA;
			delete[] allP;
			for (int i = 0; i < numNcaRep; ++i) {
				delete nmfArray[i];
			}
			delete[] nmfArray;
		} else if ((nmft == ibpg_2fct) || (nmft == KL_NMF_2fct)){
			map<int, gsl_matrix*>::iterator it;
			for (int i = 0; i < numNcaRep; ++i) {
				for(it = (allA_CT[i]).begin(); it != (allA_CT[i]).end(); ++it){
					gsl_matrix_free(it->second);
					gsl_matrix_free(allP_CT[i][it->first]);
				}
			}
			for (int i = 0; i < numNcaRep; ++i) {
				delete nmfArray[i];
			}
			delete[] nmfArray;
		}

	}
}

int
NMFLearner::setPrior(gsl_matrix* prior, double lam1)
{
	priorMat=prior;
	lambda_1 = lam1;
	return 0;
}

// cell type specific prior
int
NMFLearner::setPrior(map<int, pair<gsl_matrix*, double>>* prior, double lam1)
{
	p_net = prior;
	lambda_1 = lam1;
	return 0;
}

int
NMFLearner::setLambdaUxUGraphL(gsl_matrix* GraphL, double lam1, double norm)
{
	lambdaUxUGraphL=GraphL;
	lambda_U = lam1;
	UGraphLNorm = norm; // precomputed Spectral norm of GraphL
	return 0;
}

int
NMFLearner::setLambdaVxVGraphL(gsl_matrix* GraphL, double lam1, double norm)
{
	lambdaVxVGraphL=GraphL;
	lambda_V = lam1; 
	VGraphLNorm = norm; // precomputed Spectral norm of GraphL
	return 0;
}

void* 
NMFLearner::recur_helper(void* arguments) 
{
    struct nmf_arg_struct *args = (nmf_arg_struct *) arguments;
	NMF* bnmf = static_cast<NMF*>(args->bnmf);
	NMF* nmf_UV = static_cast<NMF*>(args->nmf_UV);
	NMF* nmf_WH = static_cast<NMF*>(args->nmf_WH);
	gsl_matrix* GsldataMat = args->datMat;
	gsl_matrix* A = args->A;
	gsl_matrix* P = args->P;
	gsl_matrix* U = args->U;
	gsl_matrix* V = args->V;
	gsl_matrix* W = args->W;
	gsl_matrix* H = args->H;
	bnmf->fit(GsldataMat, A, P);
	nmf_UV->fit(A, V, U);
	nmf_WH->fit(P, W, H);
	
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, V, U, 0, A);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, W, H, 0, P);
	gsl_matrix* R = gsl_matrix_alloc(GsldataMat->size1,GsldataMat->size2);
	gsl_matrix_memcpy(R, GsldataMat);
	double X_frob_2 = utils::get_frobenius_norm(R);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, -1, A, P, 1, R);
	double err = utils::get_frobenius_norm(R);
	gsl_matrix_free(R);
	double relative_e = sqrt(err / X_frob_2);
	(bnmf->reconstruction_err_).push_back(err);
	(bnmf->relative_err_).push_back(relative_e);
	
	printf("final iter error=%f\n",err);
	void* dummyPointer;
	return dummyPointer;
}

int
NMFLearner::inferRecur(gsl_matrix* GsldataMat)
{
	pthread_t threads[numNcaRep];
	allA = new gsl_matrix*[numNcaRep];
	allP = new gsl_matrix*[numNcaRep];
	allU = new gsl_matrix*[numNcaRep];
	allV = new gsl_matrix*[numNcaRep];
	allW = new gsl_matrix*[numNcaRep];
	allH = new gsl_matrix*[numNcaRep];
    nmfArray = new BASENMF*[numNcaRep];
	NMF** nmf_UV = new NMF*[numNcaRep];
	NMF** nmf_WH = new NMF*[numNcaRep];
	int nRegs = n_components;
	int nTgts = GsldataMat->size1;
	int nsamples = GsldataMat->size2;
	for (int i = 0; i < numNcaRep; ++i) {
		// nmf.cpp does X = U^TV
		allA[i] = gsl_matrix_calloc(nRegs, nTgts);
		allP[i] = gsl_matrix_calloc(nRegs, nsamples);
		allU[i] = gsl_matrix_calloc(k1, nTgts);
		allV[i] = gsl_matrix_calloc(k1, nRegs);
		allW[i] = gsl_matrix_calloc(k2, nRegs);
		allH[i] = gsl_matrix_calloc(k2, nsamples);
		nmfArray[i] = new NMF(nRegs, init, max_iter, random_state+i, verbose, tol, permuted_update);
		nmf_UV[i] = new NMF(k1, init, UV_maxIter, random_state+i, false, tol, permuted_update);
		nmf_WH[i] = new NMF(k2, init, WH_maxIter, random_state+i, false, tol, permuted_update);
	}
	struct nmf_arg_struct args[numNcaRep];
	for (int i = 0; i < numNcaRep; ++i) {
		args[i].A = allA[i];
		args[i].P = allP[i];
		args[i].U = allU[i];
		args[i].V = allV[i];
		args[i].W = allW[i];
		args[i].H = allH[i];
		args[i].bnmf = nmfArray[i];
		args[i].datMat = GsldataMat;
		args[i].nmf_UV = nmf_UV[i];
		args[i].nmf_WH = nmf_WH[i];
        if (pthread_create(&threads[i], NULL, recur_helper, (void *) &args[i]) != 0){
			printf("Not able to create thread\n");
			return -1;
		}
    }

	for (int i = 0; i < numNcaRep; ++i) {
        if(pthread_join(threads[i], NULL) != 0){
			printf("Error joining\n");
			return -1;
		}
    }
	ret_mat[0] = allA;
	ret_mat[1] = allP;
	ret_mat[2] = allU;
	ret_mat[3] = allV;
	ret_mat[4] = allW;
	ret_mat[5] = allH;

	for (int i = 0; i < numNcaRep; ++i) {
		delete nmf_UV[i];
		delete nmf_WH[i];
	}

	delete [] nmf_UV;
	delete [] nmf_WH;
	
	return 0;
}

void* 
NMFLearner::hals_prior_helper(void* arguments) 
{
    struct nmf_arg_struct *args = (nmf_arg_struct *) arguments;
	NMFPRIOR* bnmf = static_cast<NMFPRIOR*>(args->bnmf);
	gsl_matrix* GsldataMat = args->datMat;
	gsl_matrix* A = args->A;
	gsl_matrix* P = args->P;
	gsl_matrix* prior = args->prior;
	double lam1 = args->lambda1;
	bnmf->fit(GsldataMat, A, P, prior, lam1);

	void* dummyPointer;
	return dummyPointer;
}

int
NMFLearner::inferHals_prior(gsl_matrix* GsldataMat)
{
	pthread_t threads[numNcaRep];
	allA = new gsl_matrix*[numNcaRep];
	allP = new gsl_matrix*[numNcaRep];
    nmfArray = new BASENMF*[numNcaRep];

	int nRegs = n_components;
	int nTgts = GsldataMat->size1;
	int nsamples = GsldataMat->size2;
	for (int i = 0; i < numNcaRep; ++i) {
		// nmf.cpp does X = U^TV
		allA[i] = gsl_matrix_calloc(nRegs, nTgts);
		allP[i] = gsl_matrix_calloc(nRegs, nsamples);
		
		nmfArray[i] = new NMFPRIOR(nRegs, init, max_iter, random_state+i, verbose, tol, permuted_update);
	}
	struct nmf_arg_struct args[numNcaRep];
	for (int i = 0; i < numNcaRep; ++i) {
		args[i].A = allA[i];
		args[i].P = allP[i];
		args[i].bnmf = nmfArray[i];
		args[i].datMat = GsldataMat;
		args[i].prior = priorMat;
		args[i].lambda1 = lambda_1;

        if (pthread_create(&threads[i], NULL, hals_prior_helper, (void *) &args[i]) != 0){
			printf("Not able to create thread\n");
			return -1;
		}
    }

	for (int i = 0; i < numNcaRep; ++i) {
        if(pthread_join(threads[i], NULL) != 0){
			printf("Error joining\n");
			return -1;
		}
    }
	ret_mat[0] = allA;
	ret_mat[1] = allP;

	
	return 0;
}

void* 
NMFLearner::ibpg2f_helper(void* arguments) 
{
    struct nmf_arg_struct *args = (nmf_arg_struct *) arguments;
	IBPG_2FACTOR* bnmf = static_cast<IBPG_2FACTOR*>(args->bnmf);
	gsl_matrix* GsldataMat = args->datMat;
	gsl_matrix* A = args->A;
	gsl_matrix* P = args->P;
	gsl_matrix* prior = args->prior;
	double lam1 = args->lambda1;
	double lambda_p = args->lambdap;
	bnmf->fit(GsldataMat, A, P, prior, lam1, lambda_p);

	void* dummyPointer;
	return dummyPointer;
}

int
NMFLearner::inferIBPG2F(gsl_matrix* GsldataMat, double lambdap)
{
	pthread_t threads[numNcaRep];
	allA = new gsl_matrix*[numNcaRep];
	allP = new gsl_matrix*[numNcaRep];
    nmfArray = new BASENMF*[numNcaRep];

	int nRegs = n_components;
	int nTgts = GsldataMat->size1;
	int nsamples = GsldataMat->size2;
	for (int i = 0; i < numNcaRep; ++i) {
		// nmf.cpp does X = U^TV
		allA[i] = gsl_matrix_calloc(nTgts, nRegs);
		allP[i] = gsl_matrix_calloc(nRegs, nsamples);
		nmfArray[i] = new IBPG_2FACTOR(nRegs, init, max_iter, random_state+i, verbose, tol);
	}
	struct nmf_arg_struct args[numNcaRep];
	for (int i = 0; i < numNcaRep; ++i) {
		args[i].A = allA[i];
		args[i].P = allP[i];
		args[i].bnmf = nmfArray[i];
		args[i].datMat = GsldataMat;
		args[i].prior = priorMat;
		args[i].lambda1 = lambda_1;
		args[i].lambdap = lambdap;

        if (pthread_create(&threads[i], NULL, ibpg2f_helper, (void *) &args[i]) != 0){
			printf("Not able to create thread\n");
			return -1;
		}
    }

	for (int i = 0; i < numNcaRep; ++i) {
        if(pthread_join(threads[i], NULL) != 0){
			printf("Error joining\n");
			return -1;
		}
    }
	ret_mat[0] = allA;
	ret_mat[1] = allP;

	
	return 0;
}

void* 
NMFLearner::ibpg2fct_helper(void* arguments) 
{
    struct nmf_ct_arg_struct *args = (nmf_ct_arg_struct *) arguments;
	IBPG_2FCT* bnmf = static_cast<IBPG_2FCT*>(args->bnmf);
	map<int, gsl_matrix*>* inputmat = args->datMat;
	map<int, gsl_matrix*>* A = args->A;
	map<int, gsl_matrix*>* P = args->P;
	map<int, pair<gsl_matrix*, double>>* prior = args->prior;
	map<int, map<int, double>>* cell_neighborhood = args->cell_neighborhood;
	double lam1 = args->lambda1;
	double lambda_p = args->lambdap;
	double lambda_E = args->lambda_E;
	double lambda_Cminus = args->lambda_Cminus;
	bnmf->fit(inputmat, A, P, prior, lambda_p, cell_neighborhood, lambda_E, lambda_Cminus, lam1);

	void* dummyPointer;
	return dummyPointer;
}

int
NMFLearner::inferIBPG2FCT(map<int, gsl_matrix*>* inputmat, double lam_p,
				   map<int, map<int, double>>* cell_neighborhood, double lam_E, double lam_Cminus)
{
	pthread_t threads[numNcaRep];
	int nTgts = 0;
	int nRegs = n_components;
	nmfArray = new BASENMF*[numNcaRep];
	{
		map<int, gsl_matrix*>::iterator it;
		nTgts = ((inputmat->begin())->second)->size1;
		for (int i = 0; i < numNcaRep; ++i) {
			for(it = inputmat->begin(); it != inputmat->end(); ++it){
				int nsamples = (it->second)->size2;
				allA_CT[i][it->first] = gsl_matrix_calloc(nTgts, nRegs);
				allP_CT[i][it->first] = gsl_matrix_calloc(nRegs, nsamples);
			}
			nmfArray[i] = new IBPG_2FCT(nRegs, init, max_iter, random_state+i, verbose, tol);
		}
	}

	struct nmf_ct_arg_struct args[numNcaRep];
	for (int i = 0; i < numNcaRep; ++i) {
		args[i].A = &(allA_CT[i]);
		args[i].P = &(allP_CT[i]);
		args[i].bnmf = nmfArray[i];
		args[i].datMat = inputmat;
		args[i].prior = p_net;
		args[i].lambda1 = lambda_1;
		args[i].lambdap = lam_p;
		args[i].lambda_E = lam_E;
		args[i].lambda_Cminus = lam_Cminus;
		args[i].cell_neighborhood = cell_neighborhood;

        if (pthread_create(&threads[i], NULL, ibpg2fct_helper, (void *) &args[i]) != 0){
			printf("Not able to create thread\n");
			return -1;
		}
    }

	for (int i = 0; i < numNcaRep; ++i) {
        if(pthread_join(threads[i], NULL) != 0){
			printf("Error joining\n");
			return -1;
		}
    }
	// ret_mat[0] = allA;
	// ret_mat[1] = allP;

	
	return 0;
}

void* 
NMFLearner::KL_NMF2fct_helper(void* arguments) 
{
    struct nmf_ct_arg_struct *args = (nmf_ct_arg_struct *) arguments;
	KL_NMF_2FCT* bnmf = static_cast<KL_NMF_2FCT*>(args->bnmf);
	map<int, gsl_matrix*>* inputmat = args->datMat;
	map<int, gsl_matrix*>* A = args->A;
	map<int, gsl_matrix*>* P = args->P;
	map<int, pair<gsl_matrix*, double>>* prior = args->prior;
	map<int, map<int, double>>* cell_neighborhood = args->cell_neighborhood;
	double lam1 = args->lambda1;
	double lambda_p = args->lambdap;
	double lambda_E = args->lambda_E;
	double lambda_Cminus = args->lambda_Cminus;
	bnmf->set_total_ncarun(args->npthread);
	bnmf->fit(inputmat, A, P, prior, lambda_p, cell_neighborhood, lambda_E, lambda_Cminus, lam1);
	void* dummyPointer;
	return dummyPointer;
}

int
NMFLearner::inferKL_NMF2fct(map<int, gsl_matrix*>* inputmat, double lam_p,
				   map<int, map<int, double>>* cell_neighborhood, double lam_E, double lam_Cminus)
{
	pthread_t threads[numNcaRep];
	int nTgts = 0;
	int nRegs = n_components;
	nmfArray = new BASENMF*[numNcaRep];
	{
		map<int, gsl_matrix*>::iterator it;
		nTgts = ((inputmat->begin())->second)->size1;
		for (int i = 0; i < numNcaRep; ++i) {
			for(it = inputmat->begin(); it != inputmat->end(); ++it){
				int nsamples = (it->second)->size2;
				allA_CT[i][it->first] = gsl_matrix_calloc(nTgts, nRegs);
				allP_CT[i][it->first] = gsl_matrix_calloc(nRegs, nsamples);
			}
			nmfArray[i] = new KL_NMF_2FCT(nRegs, init, max_iter, random_state+i, verbose, tol);
		}
	}

	struct nmf_ct_arg_struct args[numNcaRep];
	for (int i = 0; i < numNcaRep; ++i) {
		args[i].A = &(allA_CT[i]);
		args[i].P = &(allP_CT[i]);
		args[i].bnmf = nmfArray[i];
		args[i].datMat = inputmat;
		args[i].prior = p_net;
		args[i].lambda1 = lambda_1;
		args[i].lambdap = lam_p;
		args[i].lambda_E = lam_E;
		args[i].lambda_Cminus = lam_Cminus;
		args[i].cell_neighborhood = cell_neighborhood;
		args[i].npthread = numNcaRep;

        if (pthread_create(&threads[i], NULL, KL_NMF2fct_helper, (void *) &args[i]) != 0){
			printf("Not able to create thread\n");
			return -1;
		}
    }

	for (int i = 0; i < numNcaRep; ++i) {
        if(pthread_join(threads[i], NULL) != 0){
			printf("Error joining\n");
			return -1;
		}
    }
	// ret_mat[0] = allA;
	// ret_mat[1] = allP;

	
	return 0;
}

void* 
NMFLearner::nmqf_helper(void* arguments) 
{
    struct nmf_arg_struct *args = (nmf_arg_struct *) arguments;
	NMQF* bnmf = static_cast<NMQF*>(args->bnmf);
	gsl_matrix* GsldataMat = args->datMat;
	gsl_matrix* A = args->A;
	gsl_matrix* P = args->P;
	gsl_matrix* U = args->U;
	gsl_matrix* V = args->V;
	gsl_matrix* W = args->W;
	gsl_matrix* H = args->H;
	gsl_matrix* prior = args->prior;
	double lam1 = args->lambda1;
	bnmf->fit(GsldataMat, U, V, W, H, prior, lam1);

	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, U, V, 0, A);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, W, H, 0, P);
	void* dummyPointer;
	return dummyPointer;
}


int
NMFLearner::inferNMQF(gsl_matrix* GsldataMat)
{
	pthread_t threads[numNcaRep];
	allA = new gsl_matrix*[numNcaRep];
	allP = new gsl_matrix*[numNcaRep];
	allU = new gsl_matrix*[numNcaRep];
	allV = new gsl_matrix*[numNcaRep];
	allW = new gsl_matrix*[numNcaRep];
	allH = new gsl_matrix*[numNcaRep];
    nmfArray = new BASENMF*[numNcaRep];
	int nRegs = n_components;
	int nTgts = GsldataMat->size1;
	int nsamples = GsldataMat->size2;
	for (int i = 0; i < numNcaRep; ++i) {
		allA[i] = gsl_matrix_calloc(nTgts, nRegs);
		allP[i] = gsl_matrix_calloc(nRegs, nsamples);
		allU[i] = gsl_matrix_calloc(nTgts, k1);
		allV[i] = gsl_matrix_calloc(nRegs, k1);
		allW[i] = gsl_matrix_calloc(nRegs, k2);
		allH[i] = gsl_matrix_calloc(nsamples, k2);
		nmfArray[i] = new NMQF(nRegs, k1, k2, init, max_iter, random_state+i, verbose, tol, permuted_update);
	}
	struct nmf_arg_struct args[numNcaRep];
	for (int i = 0; i < numNcaRep; ++i) {
		args[i].A = allA[i];
		args[i].P = allP[i];
		args[i].U = allU[i];
		args[i].V = allV[i];
		args[i].W = allW[i];
		args[i].H = allH[i];
		args[i].bnmf = nmfArray[i];
		args[i].datMat = GsldataMat;
		args[i].prior = priorMat;
		args[i].lambda1 = lambda_1;
        if (pthread_create(&threads[i], NULL, nmqf_helper, (void *) &args[i]) != 0){
			printf("Not able to create thread\n");
			return -1;
		}
    }

	for (int i = 0; i < numNcaRep; ++i) {
        if(pthread_join(threads[i], NULL) != 0){
			printf("Error joining\n");
			return -1;
		}
    }
	ret_mat[0] = allA;
	ret_mat[1] = allP;
	ret_mat[2] = allU;
	ret_mat[3] = allV;
	ret_mat[4] = allW;
	ret_mat[5] = allH;
	
	return 0;
}

void* 
NMFLearner::ibpg_helper(void* arguments) 
{
    struct nmf_arg_struct *args = (nmf_arg_struct *) arguments;
	IBPG* bnmf = static_cast<IBPG*>(args->bnmf);
	gsl_matrix* GsldataMat = args->datMat;
	
	gsl_matrix* A = args->A;
	gsl_matrix* P = args->P;
	gsl_matrix* U = args->U;
	gsl_matrix* V = args->V;
	gsl_matrix* W = args->W;
	gsl_matrix* H = args->H;
	gsl_matrix* prior = args->prior;
	double lam1 = args->lambda1;
	bnmf->setLambdaUxUGraphL((args->UGreg).lambdaxGraphL, (args->UGreg).lambda, (args->UGreg).GraphLNorm);
	bnmf->setLambdaVxVGraphL((args->VGreg).lambdaxGraphL, (args->VGreg).lambda, (args->VGreg).GraphLNorm);
	bnmf->fit(GsldataMat, U, V, W, H, prior, lam1);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, U, V, 0, A);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, W, H, 0, P);
	void* dummyPointer;
	return dummyPointer;
}


int
NMFLearner::inferIBPG(gsl_matrix* GsldataMat)
{
	pthread_t threads[numNcaRep];
	allA = new gsl_matrix*[numNcaRep];
	allP = new gsl_matrix*[numNcaRep];
	allU = new gsl_matrix*[numNcaRep];
	allV = new gsl_matrix*[numNcaRep];
	allW = new gsl_matrix*[numNcaRep];
	allH = new gsl_matrix*[numNcaRep];
    nmfArray = new BASENMF*[numNcaRep];
	int nRegs = n_components;
	int nTgts = GsldataMat->size1;
	int nsamples = GsldataMat->size2;
	for (int i = 0; i < numNcaRep; ++i) {
		allA[i] = gsl_matrix_calloc(nTgts, nRegs);
		allP[i] = gsl_matrix_calloc(nRegs, nsamples);
		allU[i] = gsl_matrix_calloc(nTgts, k1);
		allV[i] = gsl_matrix_calloc(nRegs, k1);
		allW[i] = gsl_matrix_calloc(nRegs, k2);
		allH[i] = gsl_matrix_calloc(nsamples, k2);
		nmfArray[i] = new IBPG(nRegs, k1, k2, init, max_iter, random_state+i, verbose, tol, permuted_update);
	}
	struct nmf_arg_struct args[numNcaRep];
	for (int i = 0; i < numNcaRep; ++i) {
		args[i].A = allA[i];
		args[i].P = allP[i];
		args[i].U = allU[i];
		args[i].V = allV[i];
		args[i].W = allW[i];
		args[i].H = allH[i];
		args[i].bnmf = nmfArray[i];
		args[i].datMat = GsldataMat;
		args[i].prior = priorMat;
		args[i].lambda1 = lambda_1;
		args[i].UGreg.lambdaxGraphL = lambdaUxUGraphL;
		args[i].UGreg.lambda = lambda_U;
		args[i].UGreg.GraphLNorm = UGraphLNorm;
		args[i].VGreg.lambdaxGraphL = lambdaVxVGraphL;
		args[i].VGreg.lambda = lambda_V;
		args[i].VGreg.GraphLNorm = VGraphLNorm;
        if (pthread_create(&threads[i], NULL, ibpg_helper, (void *) &args[i]) != 0){
			printf("Not able to create thread\n");
			return -1;
		}
    }

	for (int i = 0; i < numNcaRep; ++i) {
        if(pthread_join(threads[i], NULL) != 0){
			printf("Error joining\n");
			return -1;
		}
    }
	ret_mat[0] = allA;
	ret_mat[1] = allP;
	ret_mat[2] = allU;
	ret_mat[3] = allV;
	ret_mat[4] = allW;
	ret_mat[5] = allH;
	
	return 0;
}


int
NMFLearner::inferNMF(int numRep, gsl_matrix* GsldataMat)
{
	if(numRep <=0){
		printf("Invalid number of reps: %d\n", numRep);
		return -1;
	}
	numNcaRep = numRep;

	if(nmft == recur){
		inferRecur(GsldataMat);
	} else if (nmft == nmqf){
		inferNMQF(GsldataMat);
	} else if (nmft == ibpg){
		inferIBPG(GsldataMat);
	} else if (nmft == hals_prior) {
		inferHals_prior(GsldataMat);
	} else if (nmft == ibpg_2f) {
		inferIBPG2F(GsldataMat, prior_info_public.lambda_p);
	} else {
		return -1;
	}
	
	return 0;
}

// cell type specific
int
NMFLearner::inferNMF_ct(int numRep, map<int, gsl_matrix*>* inputmat, double lam_p,
				   map<int, map<int, double>>* cell_neighborhood, double lam_E, double lam_Cminus)
{
	if(numRep <=0){
		printf("Invalid number of reps: %d\n", numRep);
		return -1;
	}
	numNcaRep = numRep;

	if (nmft == ibpg_2fct) {
		inferIBPG2FCT(inputmat, lam_p, cell_neighborhood, lam_E, lam_Cminus);
	} else if (nmft == KL_NMF_2fct) {
		inferKL_NMF2fct(inputmat, lam_p, cell_neighborhood, lam_E, lam_Cminus);
	} else {
		return -1;
	}
	
	return 0;
}

// pair<map<int, gsl_matrix**>*, map<int, vector<double>*>*>*
// NMFLearner::get_nmf_out()
// {
// 	return &factor_mat_w_loss;
// }

int
NMFLearner::writeToDir(char* outname)
{
	int status = 0;
	char mkdirCmd[1024];
	FILE* f;
	for (int i = 0; i < numNcaRep; ++i) {
		sprintf(mkdirCmd,"mkdir -p %s/rep_%d", outname, i);
		system(mkdirCmd);

		char aoutname[1024];

		// print metrics
		sprintf(aoutname,"%s/rep_%d/squared_loss.txt", outname, i);
		f = fopen(aoutname,"w");
		list<double>* cur_rl = &(nmfArray[i]->reconstruction_err_);
		for (list<double>::iterator it = cur_rl->begin(); it != cur_rl->end(); ++it) {
        	if((status=fprintf(f, "%f\n", *it))<0)
				return -1;
    	}
		fclose(f);

		sprintf(aoutname,"%s/rep_%d/relative_error.txt", outname, i);
		f = fopen(aoutname,"w");
		cur_rl = &(nmfArray[i]->relative_err_);
		for (list<double>::iterator it = cur_rl->begin(); it != cur_rl->end(); ++it) {
        	if((status=fprintf(f, "%f\n", *it))<0)
				return -1;
    	}
		fclose(f);

		sprintf(aoutname,"%s/rep_%d/cpu_time.txt", outname, i);
		f = fopen(aoutname,"w");
		cur_rl = &(nmfArray[i]->cpu_time_);
		for (list<double>::iterator it = cur_rl->begin(); it != cur_rl->end(); ++it) {
        	if((status=fprintf(f, "%f\n", *it))<0)
				return -1;
    	}
		fclose(f);
		
		
		sprintf(aoutname,"%s/rep_%d/A.txt", outname, i);
		f = fopen(aoutname,"w");
		print_matrix(f, allA[i]);
		fclose(f);
		sprintf(aoutname,"%s/rep_%d/P.txt", outname, i);
		f = fopen(aoutname,"w");
		print_matrix(f, allP[i]);
		fclose(f);

		if ((nmft == hals_prior) || (nmft == ibpg_2f)){ // only two factors
			continue;
		}

		sprintf(aoutname,"%s/rep_%d/U.txt", outname, i);
		f = fopen(aoutname,"w");
		print_matrix(f, allU[i]);
		fclose(f);
		sprintf(aoutname,"%s/rep_%d/V.txt", outname, i);
		f = fopen(aoutname,"w");
		print_matrix(f, allV[i]);
		fclose(f);
		sprintf(aoutname,"%s/rep_%d/W.txt", outname, i);
		f = fopen(aoutname,"w");
		print_matrix(f, allW[i]);
		fclose(f);
		sprintf(aoutname,"%s/rep_%d/H.txt", outname, i);
		f = fopen(aoutname,"w");
		print_matrix(f, allH[i]);
		fclose(f);

	}
	return 0;
}

int
NMFLearner::writeToDir(char* outname, map<int, string>* IDToCell){
	int status = 0;
	char mkdirCmd[1024];
	FILE* f;
	for (int i = 0; i < numNcaRep; ++i) {
		sprintf(mkdirCmd,"mkdir -p %s/rep_%d", outname, i);
		system(mkdirCmd);

		char aoutname[1024];

		// print metrics
		sprintf(aoutname,"%s/rep_%d/squared_loss.txt", outname, i);
		f = fopen(aoutname,"w");
		list<double>* cur_rl = &(nmfArray[i]->reconstruction_err_);
		for (list<double>::iterator it = cur_rl->begin(); it != cur_rl->end(); ++it) {
        	if((status=fprintf(f, "%f\n", *it))<0)
				return -1;
    	}
		fclose(f);

		sprintf(aoutname,"%s/rep_%d/relative_error.txt", outname, i);
		f = fopen(aoutname,"w");
		cur_rl = &(nmfArray[i]->relative_err_);
		for (list<double>::iterator it = cur_rl->begin(); it != cur_rl->end(); ++it) {
        	if((status=fprintf(f, "%f\n", *it))<0)
				return -1;
    	}
		fclose(f);

		sprintf(aoutname,"%s/rep_%d/cpu_time.txt", outname, i);
		f = fopen(aoutname,"w");
		cur_rl = &(nmfArray[i]->cpu_time_);
		for (list<double>::iterator it = cur_rl->begin(); it != cur_rl->end(); ++it) {
        	if((status=fprintf(f, "%f\n", *it))<0)
				return -1;
    	}
		fclose(f);
		
		map<int, string>::iterator it;
		for(it = IDToCell->begin(); it != IDToCell->end(); ++it)
		{

			sprintf(mkdirCmd,"mkdir -p %s/rep_%d/%s", outname, i, (it->second).c_str());
			system(mkdirCmd);
			sprintf(aoutname,"%s/rep_%d/%s/A.txt", outname, i, (it->second).c_str());
			f = fopen(aoutname,"w");
			print_matrix(f, allA_CT[i][it->first]);
			fclose(f);
			sprintf(aoutname,"%s/rep_%d/%s/P.txt", outname, i, (it->second).c_str());
			f = fopen(aoutname,"w");
			print_matrix(f, allP_CT[i][it->first]);
			fclose(f);
		}
		
	}
	return 0;
}

int 
NMFLearner::print_matrix(FILE *f, const gsl_matrix *m)
{
	int status, n = 0;
	for (size_t i = 0; i < m->size1; i++) {
		for (size_t j = 0; j < m->size2; j++) {
			if ((status = fprintf(f, "%g ", gsl_matrix_get(m, i, j))) < 0)
					return -1;
			n += status;
		}
		if ((status = fprintf(f, "\n")) < 0)
			return -1;
		n += status;
	}
	return n;
}