#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_eigen.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <ctime>

#include "initialization.h"
#include "ibpg_2factor.h"
#include "utils.h"

IBPG_2FACTOR::IBPG_2FACTOR(int k, init_method initMethod, int maxIter, int seed, bool verb, double termTol){
	n_components = k;
	init = initMethod;
	max_iter = maxIter;
	random_state = seed;
	verbose = verb;
	tol = termTol;
	list<double> reconstruction_err_;
	list<double> relative_err_;
	list<double> cpu_time_;
	beta = 1.01; // control the influence of the proximal point term
	deltaw = 0.99; // control the gradient term. Also control the proximal point term. set to 0 and you will get gradient descent. Set to 0.99 for IBPG-A
}

IBPG_2FACTOR::~IBPG_2FACTOR() {};


int IBPG_2FACTOR::record_progress(double rc_l, double relative_e, double time){
	reconstruction_err_.push_back(rc_l);
	relative_err_.push_back(relative_e);
	cpu_time_.push_back(time);
	return 0;
}

int IBPG_2FACTOR::update_A(){
	gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1, P, 0, PPT);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, X, P, 0, gradAy); 

	gsl_vector* eval = gsl_vector_alloc(n_components); 
	gsl_matrix* tmp_PPT = gsl_matrix_alloc(n_components, n_components); 
	gsl_matrix_memcpy(tmp_PPT, PPT);
    gsl_eigen_symm_workspace* w = gsl_eigen_symm_alloc(n_components); 
    gsl_eigen_symm(tmp_PPT, eval, w); //lower triangular + diagonal destroyed
    gsl_eigen_symm_free(w);
	gsl_matrix_free(tmp_PPT);

	double L_A = max(gsl_vector_max(eval), 0.01); 
	gsl_vector_free(eval);

	double gamma;
	if(A_first_it){
		gamma = what; 
	} else {
		gamma = min( what, deltaw*sqrt( L_A_prev/L_A ) ); 
	}
	int j = 0;

	while(j < inneriterA){
		gsl_matrix* tmp_An = A;
		A = Aold;
		gsl_matrix_sub(A, tmp_An); // -1 * change in A
		gsl_matrix_scale(A, -1*gamma); 
		gsl_matrix_memcpy(Ay, tmp_An);
		gsl_matrix_add(Ay, A);
        gsl_blas_dsymm(CblasRight, CblasLower, 1, PPT, Ay, -1, gradAy); 
		gsl_matrix_add(gradAy, Prior);

		gsl_matrix_scale(gradAy, 1/L_A);
		gsl_matrix_scale(A, beta);
		gsl_matrix_sub(A, gradAy);
		gsl_matrix_add(A, tmp_An);

        size_t tda = A->tda;
        double* data = A->data;
		for(int r = 0; r < A->size1; ++r){
			for(int c = 0; c < n_components; ++c){
				double* val =  &data[r * tda + c];
				if (*val < 0) {
					*val = 0;
				}
			}
		}
		Aold = tmp_An;
		++j;
	}
	// for(int k = 0; k < k1; ++k){
	// 	gsl_vector_view u_k = gsl_matrix_column(U, k);
	// 	if (gsl_vector_isnull(&u_k.vector)) {
	// 		gsl_vector_set_all(&u_k.vector, 0.01);
	// 	}
	// }
	L_A_prev = L_A; // update gradient Lipschitz constant

    return 0;
}
		
int IBPG_2FACTOR::update_P(){
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1, A, 0, ATA);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, A, X, 0, gradPy); 

	gsl_vector* eval = gsl_vector_alloc(n_components); 
	gsl_matrix* tmp_ATA = gsl_matrix_alloc(n_components, n_components); 
	gsl_matrix_memcpy(tmp_ATA, ATA);
    gsl_eigen_symm_workspace* w = gsl_eigen_symm_alloc(n_components); 
    gsl_eigen_symm(tmp_ATA, eval, w); //lower triangular + diagonal destroyed
    gsl_eigen_symm_free(w);
	gsl_matrix_free(tmp_ATA);

	double L_P = max(gsl_vector_max(eval), 0.01); 
	gsl_vector_free(eval);

	double gamma;
	if(A_first_it){
		gamma = what; 
	} else {
		gamma = min( what, deltaw*sqrt( L_P_prev/L_P ) ); 
	}
	int j = 0;

	while(j < inneriterP){
		gsl_matrix* tmp_Pn = P;
		P = Pold;
		gsl_matrix_sub(P, tmp_Pn); // -1 * change in A
		gsl_matrix_scale(P, -1*gamma); 
		gsl_matrix_memcpy(Py, tmp_Pn);
		gsl_matrix_add(Py, P);
        gsl_blas_dsymm(CblasLeft, CblasLower, 1, ATA, Py, -1, gradPy); 
		gsl_matrix_add_constant(gradPy, lambda_p);

		gsl_matrix_scale(gradPy, 1/L_P);
		gsl_matrix_scale(P, beta);
		gsl_matrix_sub(P, gradPy);
		gsl_matrix_add(P, tmp_Pn);

        size_t tda = P->tda;
        double* data = P->data;
		for(int r = 0; r < n_components; ++r){
			for(int c = 0; c < P->size2; ++c){
				double* val =  &data[r * tda + c];
				if (*val < 0) {
					*val = 0;
				}
			}
		}
		Pold = tmp_Pn;
		++j;
	}
	// for(int k = 0; k < k1; ++k){
	// 	gsl_vector_view u_k = gsl_matrix_column(U, k);
	// 	if (gsl_vector_isnull(&u_k.vector)) {
	// 		gsl_vector_set_all(&u_k.vector, 0.01);
	// 	}
	// }
	L_P_prev = L_P; // update gradient Lipschitz constant
    return 0;
}

int IBPG_2FACTOR::update(){
    update_A();
    update_P();
    return 0;
}

double IBPG_2FACTOR::calculate_objective(){
    gsl_matrix_memcpy(R, X);
	// A and P are calculated at the beginning of fit() and in the middle of update()
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1, A, P, 1, R);
	double sum = utils::get_frobenius_norm(R);
    cur_relative_err = sqrt(sum / X_frob);
    sum *= 0.5;
    size_t rows = A->size1;
    double* data = A->data; // Pointer to the raw data array
    double* pdata = Prior->data;
    size_t tda = A->tda;
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < n_components; ++j) {
            sum += data[i * tda + j] * pdata[i * tda + j];
        }
    }
    data = P-> data;
    tda = P->tda;
    size_t cols = P->size2;
    for (size_t i = 0; i < n_components; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            sum += data[i * tda + j] * lambda_p;
        }
    }
    return sum;
}

int IBPG_2FACTOR::fit(gsl_matrix* inputmat, gsl_matrix* A_in, gsl_matrix* P_in, gsl_matrix* p_net, const double lam_1, const double lam_p){
    X = inputmat;
	X_frob = utils::get_frobenius_norm(X);
	int n = X->size1;
	int m = X->size2;
	A = A_in;
	P = P_in;
    if ((A->size2 != n_components) || (P->size1 != n_components)) {
		cout << "The inner dimension of the two factor matrices should agree with the supplied number of factors (regulators)." << endl;
		return 1;
	}

    if (verbose) {
		cout << "Initializing..." << endl;
	}	
	if (init == random_init) {
		init::random(A, P, random_state);
	} else {
        gsl_matrix* H = gsl_matrix_calloc(m, n_components);
		init::nndsvd2(X, A, H, random_state);
        gsl_matrix_transpose_memcpy(P, H);
        gsl_matrix_free(H);
	}

    R = gsl_matrix_alloc(n,m);

    //For updating A
    Aold = gsl_matrix_alloc(n, n_components);
    gsl_matrix_memcpy(Aold, A); 
	PPT = gsl_matrix_calloc(n_components, n_components);
	gradAy = gsl_matrix_calloc(n, n_components);
	Ay = gsl_matrix_calloc(n, n_components);

    //For updating P
    Pold = gsl_matrix_alloc(n_components, m);
    gsl_matrix_memcpy(Pold, P);
	ATA = gsl_matrix_calloc(n_components, n_components);
	gradPy = gsl_matrix_calloc(n_components, m);
	Py = gsl_matrix_calloc(n_components, m);

    lambda_p = lam_p;
    Prior = gsl_matrix_alloc(n, n_components);
    gsl_matrix_memcpy(Prior, p_net); // note that p_net has entries of edge confidence ranged from 0 to 1
	gsl_matrix_add_constant(Prior, -1); 
	gsl_matrix_scale(Prior, -1*lam_1); // make edge confidence to penalty terms. Edge not appearing or with low confidence will be penalized

    // TODO change this
    inneriterA = 1;
    inneriterP = 1;
    tau = 1;
	prev_tau = 1;

    // Below ensures that the right matrices are pointed to by the input pointers, not discarded.
    gsl_matrix* Aold_tmp = Aold;
    gsl_matrix* Pold_tmp = Pold; 

    const gsl_rng_type* T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	ri = gsl_rng_alloc(T);
	gsl_rng_set(ri, random_state);
	clock_t c_start = clock();


    int below_tol_it = 30; // change this
	int below_tol_it_counter = 1;
    A_first_it = true;
    double old_error = calculate_objective();
    double old_slope;
    record_progress(old_error, cur_relative_err, 0);
    for (int n_iter =0; n_iter < max_iter; ++n_iter){
        tau = 0.5 * (1+ sqrt( 1 + 4*pow(prev_tau, 2) ));
		what = (prev_tau-1)/tau;
		prev_tau = tau;
		update();
		clock_t c_cur = clock();
        A_first_it = false;
		double error = calculate_objective();
		double slope = old_error - error;
		double elapse = (c_cur-c_start) / CLOCKS_PER_SEC;
		record_progress(error, cur_relative_err, elapse);
		if (verbose) {
			cout << "Itr " << n_iter+1 << " err2 = " << error << ", slope = " << slope 
				<< ", relative = " << cur_relative_err << ", cpu(s) = " << elapse << endl;
		}
		if (slope < tol) {
			if (below_tol_it_counter >= below_tol_it){
				if (verbose) {
					cout << "Converged at iteration " << n_iter+1 << endl;	
				}
				break;
			}
			++below_tol_it_counter;
		} else {
			old_error = error;
			old_slope = slope;
		}
	}
    // Ensure that input pointers point to the correct objects
	if(Aold != Aold_tmp){
		gsl_matrix_memcpy(Aold, A);
		Aold = Aold_tmp;
	}
	if(Pold != Pold_tmp){
		gsl_matrix_memcpy(Pold, P);
		Pold = Pold_tmp;
	}


    gsl_matrix_free(Prior);
    // Free A
    gsl_matrix_free(Aold);
	gsl_matrix_free(PPT);
	gsl_matrix_free(gradAy);
	gsl_matrix_free(Ay);
    // Free P
    gsl_matrix_free(Pold);
	gsl_matrix_free(ATA);
	gsl_matrix_free(gradPy);
	gsl_matrix_free(Py);

    gsl_rng_free(ri);
    return 0;
}