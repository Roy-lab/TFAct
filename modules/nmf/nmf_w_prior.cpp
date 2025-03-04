#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <ctime>

#include "initialization.h"
#include "nmf_w_prior.h"
#include "utils.h"

NMFPRIOR::NMFPRIOR(int k, init_method initMethod, int maxIter, int seed, bool verb, double termTol) {
	n_components = k;
	init = initMethod;
	max_iter = maxIter;
	random_state = seed;
	verbose = verb;
	tol = termTol;
	list<double> reconstruction_err_;
	list<double> relative_err_;
	list<double> cpu_time_;
	update_permuted = false;
}

NMFPRIOR::NMFPRIOR(int k, init_method initMethod, int maxIter, int seed, bool verb, double termTol, bool perm_up){
	n_components = k;
	init = initMethod;
	max_iter = maxIter;
	random_state = seed;
	verbose = verb;
	tol = termTol;
	list<double> reconstruction_err_;
	update_permuted = perm_up;
}

NMFPRIOR::~NMFPRIOR() {};

int NMFPRIOR::update_kth_block_of_U(int k){
	gsl_vector_view u_k = gsl_matrix_row(U, k);
	gsl_vector_view v_k = gsl_matrix_row(V, k);
	gsl_vector_view p_k = gsl_matrix_row(Prior, k);
	double v_norm = pow(gsl_blas_dnrm2(&v_k.vector), 2);
	gsl_blas_dgemv(CblasNoTrans, 1, R, &v_k.vector, 0, &u_k.vector);
	gsl_vector_sub(&u_k.vector, &p_k.vector);
	gsl_vector_scale(&u_k.vector, 1/v_norm);

	int n = R->size1;
	for (int i = 0; i < n; ++i) {
		double *val = &((&u_k.vector)->data[i]);
		if (*val < 0) {
			*val = 0;
		}
	}
	if (gsl_vector_isnull(&u_k.vector)) {
		gsl_vector_set_all(&u_k.vector, 0.01);
	}
	return 0;
}

int NMFPRIOR::update_kth_block_of_V(int k){
	gsl_vector_view u_k = gsl_matrix_row(U, k);
	gsl_vector_view v_k = gsl_matrix_row(V, k);
	double u_norm = pow(gsl_blas_dnrm2(&u_k.vector), 2);
	gsl_blas_dgemv(CblasTrans, 1/u_norm, R, &u_k.vector, 0, &v_k.vector);
	int m = R->size2;
	for (int i = 0; i < m; ++i) {
		double *val = &((&v_k.vector)->data[i]);
		if (*val < 0) {
			*val = 0;
		}
	}
	if (gsl_vector_isnull(&v_k.vector)) {
		gsl_vector_set_all(&v_k.vector, 0.01);
	}
return 0;
}

int NMFPRIOR::update() {
	int update_order[n_components];
	for (int k = 0; k < n_components; ++k) {
		update_order[k] = k;
	}
	if(update_permuted){
		gsl_ran_shuffle(ri, update_order, n_components, sizeof(int));
	}
	int k;
	for (int i = 0; i < n_components; ++i) {
		k = update_order[i];
		gsl_vector_view u_k = gsl_matrix_row(U, k);
		gsl_vector_view v_k = gsl_matrix_row(V, k);
		gsl_blas_dger(1, &u_k.vector, &v_k.vector, R);
		update_kth_block_of_U(k);
		update_kth_block_of_V(k);
		gsl_blas_dger(-1, &u_k.vector, &v_k.vector, R);
	}
    normalize_and_scale();
	return 0;
}

int NMFPRIOR::normalize_and_scale() {
	for (int k = 0; k < n_components; ++k) {
		gsl_vector_view u_k = gsl_matrix_row(U, k);
		gsl_vector_view v_k = gsl_matrix_row(V, k);
		double norm = gsl_blas_dnrm2(&v_k.vector);
		gsl_vector_scale(&v_k.vector, 1/norm);
		gsl_vector_scale(&u_k.vector, norm);
	}
	return 0;
}

double NMFPRIOR::calculate_objective() {
	gsl_matrix_memcpy(R, X);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, -1, U, V, 1, R);
	double error = utils::get_frobenius_norm(R);
	// reconstruction_err_.push_back(error);
	return error;
}

int NMFPRIOR::record_progress(double rc_l, double relative_e, double time){
	reconstruction_err_.push_back(rc_l);
	relative_err_.push_back(relative_e);
	cpu_time_.push_back(time);
	return 0;
}

int NMFPRIOR::fit(gsl_matrix* inputmat, gsl_matrix* W, gsl_matrix* H,  gsl_matrix* p_net, const double lam_1) {
	X = inputmat;
	double X_frob = utils::get_frobenius_norm(X);
	int n = X->size1;
	int m = X->size2;
	U = W;
	V = H;
	R = gsl_matrix_alloc(n,m);

	if ((U->size1 != n_components) || (V->size1 != n_components)) {
		cout << "The first dimension of U and V (i.e. their number of rows) should equal the number of components specified when instantiating NMF." << endl;
		return 1;
	} 

	if (verbose) {
		cout << "Initializing..." << endl;
	}	
	if (init == random_init) {
		init::random(U, V, random_state);
	} else {
		init::nndsvd(X, U, V, random_state);
	}
	
	Prior = gsl_matrix_alloc(n_components, n);
	gsl_matrix_transpose_memcpy(Prior, p_net);
	gsl_matrix_add_constant(Prior, -1);
	gsl_matrix_scale(Prior, -1*lam_1);


	double old_error = calculate_objective();
	if(verbose){
		printf("initial error: %f\n", old_error);
	}
	double old_slope;
	{
		double relativeErr = sqrt(old_error / X_frob);
		record_progress(old_error, relativeErr, 0);
	}
	const gsl_rng_type* T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	ri = gsl_rng_alloc(T);
	gsl_rng_set(ri, random_state);
	clock_t c_start = clock();

	int below_tol_it = 10; // change this
	int below_tol_it_counter = 1;
	
	for (int n_iter =0; n_iter < max_iter; n_iter++){
		update();
		clock_t c_cur = clock();
		double error = calculate_objective();
		double slope = old_error - error;
		double relativeErr = sqrt(error / X_frob);
		double elapse = (c_cur-c_start) / CLOCKS_PER_SEC;
		record_progress(error, relativeErr, elapse);
		if (verbose) {
			cout << "Itr " << n_iter+1 << " err2 = " << error << ", slope = " << slope 
				<< ", relative = " << relativeErr << ", cpu(s) = " << elapse << endl;
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
	gsl_matrix_free(R);
	gsl_rng_free(ri);

	/*
	for (int i = 0; i < n_components; ++i) {
		for (int j = 0; j < 10; j++) {
			cout << U->data[i * U->tda + j] << ", ";
		}
		cout << endl;
	}
	*/

	return 0;
}
