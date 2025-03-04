#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <ctime>

#include "initialization.h"
#include "nmqf.h"
#include "utils.h"

NMQF::NMQF(int k, init_method initMethod, int maxIter, int seed, bool verb, double termTol) {
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

NMQF::NMQF(int k, int k_one, int k_two, init_method initMethod, int maxIter, int seed, bool verb, double termTol, bool perm_up) {
	n_components = k;
	k1 = k_one;
	k2 = k_two;
	init = initMethod;
	max_iter = maxIter;
	random_state = seed;
	verbose = verb;
	tol = termTol;
	list<double> reconstruction_err_;
	list<double> relative_err_;
	update_permuted = perm_up;
	maxInnerDim = max(k, k1);
}

NMQF::~NMQF() {};

int NMQF::update_kth_block_of_U(int k){
	gsl_vector_view u_k = gsl_matrix_column(U, k);
	double q_norm = pow(gsl_blas_dnrm2(q), 2);
	if(regA){
		gsl_vector_view v_k = gsl_matrix_column(V, k);
		gsl_blas_dgemv(CblasNoTrans, -1/max(q_norm,0.01), Prior, &v_k.vector, 0, &u_k.vector);
		gsl_blas_dgemv(CblasNoTrans, 1/max(q_norm,0.01), R, q, 1, &u_k.vector);
	} else {
		gsl_blas_dgemv(CblasNoTrans, 1/max(q_norm,0.01), R, q, 0, &u_k.vector);
	}
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

// int NMQF::update_kth_block_of_V(int k){
// 	gsl_vector_view u_k = gsl_matrix_column(U, k);
// 	gsl_vector_view v_k = gsl_matrix_column(V, k);
// 	double u_norm = pow(gsl_blas_dnrm2(&u_k.vector), 2);
// 	gsl_blas_dgemv(CblasTrans, 1/u_norm, R, &u_k.vector, 0, &v_k.vector);
// 	int m = R->size2;
// 	for (int i = 0; i < m; i++) {
// 		double *val = &((&v_k.vector)->data[i]);
// 		if (*val < 0) {
// 			*val = 0;
// 		}
// 	}
// 	if (gsl_vector_isnull(&v_k.vector)) {
// 		gsl_vector_set_all(&v_k.vector, 0.01);
// 	}
// 	return 0;
// }

int NMQF::update_kth_block_of_H(int k){
	gsl_vector_view h_k = gsl_matrix_column(H, k);
	double p_norm = pow(gsl_blas_dnrm2(p), 2);
	gsl_blas_dgemv(CblasTrans, 1/max(p_norm,0.01), R, p, 0, &h_k.vector);
	int m = R->size2;
	for (int i = 0; i < m; i++) {
		double *val = &((&h_k.vector)->data[i]);
		if (*val < 0) {
			*val = 0;
		}
	}
	if (gsl_vector_isnull(&h_k.vector)) {
		gsl_vector_set_all(&h_k.vector, 0.01);
	}
	return 0;
}

int NMQF::update_ij_entry_of_V(int i, int j){
	gsl_vector_view u_l = gsl_matrix_column(U, j);
	gsl_vector_view P_k_row = gsl_matrix_row(P, i);
	double* Vij = gsl_matrix_ptr(V, i, j);
	double u_norm = pow(gsl_blas_dnrm2(&u_l.vector), 2);
	double p_norm = pow(gsl_blas_dnrm2(&P_k_row.vector), 2);
	gsl_blas_dgemv(CblasNoTrans, 1, R, &P_k_row.vector, 0, p);
	if(regA){
		gsl_vector_view Prior_k_col = gsl_matrix_column(Prior, i);
		gsl_vector_sub(p, &Prior_k_col.vector);
	}
	gsl_blas_ddot(&u_l.vector, p, Vij);
	*Vij /= max(u_norm*p_norm, 0.01);
	if (*Vij < 0) {
		*Vij = 0;
	}
	return 0;
}

int NMQF::update_ij_entry_of_W(int i, int j){
	gsl_vector_view a_k = gsl_matrix_column(A, i);
	gsl_vector_view h_l = gsl_matrix_column(H, j);
	double* Wij = gsl_matrix_ptr(W, i, j);
	double a_norm = pow(gsl_blas_dnrm2(&a_k.vector), 2);
	double h_norm = pow(gsl_blas_dnrm2(&h_l.vector), 2);
	gsl_blas_dgemv(CblasNoTrans, 1, R, &h_l.vector, 0, p);
	gsl_blas_ddot(&a_k.vector, p, Wij);
	*Wij /= max(h_norm*a_norm, 0.01);
	if (*Wij < 0) {
		*Wij = 0;
	}
	return 0;
}

int NMQF::update() {
	int update_order[maxInnerDim];
	int update_order_component[n_components]; // for choosing Vkl, Wkl
	for (int k = 0; k < k1; ++k) {
		update_order[k] = k;
	}
	for (int k = 0; k < n_components; ++k) {
		update_order_component[k] = k;
	}
	if(update_permuted){
		gsl_ran_shuffle(ri, update_order, k1, sizeof(int));
	}

	for (int i = 0; i < k1; ++i) { //update U
		int k = update_order[i];
		gsl_vector_view u_k = gsl_matrix_column(U, k);
		gsl_vector_view v_k = gsl_matrix_column(V, k);
		gsl_blas_dgemv(CblasTrans, 1, P, &v_k.vector, 0, q);
		gsl_blas_dger(1, &u_k.vector, q, R);
		update_kth_block_of_U(k);
		// update_kth_block_of_V(k);
		gsl_blas_dger(-1, &u_k.vector, q, R);
	}
	if(update_permuted){
		gsl_ran_shuffle(ri, update_order, k1, sizeof(int));
		gsl_ran_shuffle(ri, update_order_component, n_components, sizeof(int));
	}
	for (int j = 0; j < k1; ++j) {//update Vkl
		int l = update_order[j];
		for (int i = 0; i < n_components; ++i) { 
			int k = update_order_component[i];
			gsl_vector_view u_l = gsl_matrix_column(U, l);
			gsl_vector_view P_k_row = gsl_matrix_row(P, k);
			double* Vkl = gsl_matrix_ptr(V, k, l);
			gsl_blas_dger(*Vkl, &u_l.vector, &P_k_row.vector, R);
			update_ij_entry_of_V(k,l);
			gsl_blas_dger(-1*(*Vkl), &u_l.vector, &P_k_row.vector, R);
		}
		gsl_vector_view v_k = gsl_matrix_column(V, l);
		if (gsl_vector_isnull(&v_k.vector)) {// ensure that dividing by V's column norm is defined
			gsl_vector_set_all(&v_k.vector, 0.001);
		}
	}
	
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, U, V, 0, A); // UV are done

	for (int k = 0; k < k2; ++k) {
		update_order[k] = k;
	}
	if(update_permuted){
		gsl_ran_shuffle(ri, update_order, k2, sizeof(int));
		gsl_ran_shuffle(ri, update_order_component, n_components, sizeof(int));
	}
	for (int i = 0; i < k2; ++i) { //update H
		int k = update_order[i];
		gsl_vector_view w_k = gsl_matrix_column(W, k);
		gsl_vector_view h_k = gsl_matrix_column(H, k);
		gsl_blas_dgemv(CblasNoTrans, 1, A, &w_k.vector, 0, p);
		gsl_blas_dger(1, p, &h_k.vector, R);
		// update_kth_block_of_U(k);
		update_kth_block_of_H(k);
		gsl_blas_dger(-1, p, &h_k.vector, R);
	}
	for (int j = 0; j < k2; ++j) {//update Wkl
		int l = update_order[j];
		for (int i = 0; i < n_components; ++i) { 
			int k = update_order_component[i];
			gsl_vector_view a_k = gsl_matrix_column(A, k);
			gsl_vector_view h_l = gsl_matrix_column(H, l);
			double* Wkl = gsl_matrix_ptr(W, k, l);
			gsl_blas_dger(*Wkl, &a_k.vector, &h_l.vector, R);
			update_ij_entry_of_W(k,l);
			gsl_blas_dger(-1*(*Wkl), &a_k.vector, &h_l.vector, R);
		}
		gsl_vector_view w_l = gsl_matrix_column(W, l);
		if (gsl_vector_isnull(&w_l.vector)) {// ensure that dividing by W's column norm is defined
			gsl_vector_set_all(&w_l.vector, 0.001);
		}
	}
	
	//// The below line is performed in normalize_and_scale. If don't want to normalize and scale, remove comment for below
	// gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, W, H, 0, P); // WH are done
    normalize_and_scale();
	return 0;
}

double NMQF::get_absolute_sum(const gsl_matrix* matrix){
	double sum = 0.0;
    size_t rows = matrix->size1;
    size_t cols = matrix->size2;
    double* data = matrix->data;  
    size_t tda = matrix->tda;

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            double value = data[i * tda + j];  
            sum += abs(value); // actually abs is not needed as all values are nonnegative
        }
    }
	return sum;
}

int NMQF::normalize_and_scale() { 
	// redistribute weights for WH^T 
	for (int k = 0; k < k2; ++k) {
		gsl_vector_view w_k = gsl_matrix_column(W, k);
		gsl_vector_view h_k = gsl_matrix_column(H, k);
		double norm = gsl_blas_dnrm2(&h_k.vector);
		gsl_vector_scale(&h_k.vector, 1/norm);
		gsl_vector_scale(&w_k.vector, norm);
	}

	// For UV^T
	for (int k = 0; k < k1; ++k) {
		gsl_vector_view u_k = gsl_matrix_column(U, k);
		gsl_vector_view v_k = gsl_matrix_column(V, k);
		double norm = gsl_blas_dnrm2(&u_k.vector);
		gsl_vector_scale(&u_k.vector, 1/norm);
		gsl_vector_scale(&v_k.vector, norm);
	}

	// scale the full matrices 
	double W_abs_sum = get_absolute_sum(W);
	double V_abs_sum = get_absolute_sum(W);
	double scaling_factor = sqrt(W_abs_sum/V_abs_sum);
	gsl_matrix_scale(W, 1/scaling_factor);
	gsl_matrix_scale(V, scaling_factor);
	
	// Update A and P
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, U, V, 0, A);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, W, H, 0, P);
	return 0;
}

double NMQF::calculate_objective() {
	gsl_matrix_memcpy(R, X);
	// A and P are calculated at the beginning of fit() and in the middle of update()
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1, A, P, 1, R);
	double error = utils::get_frobenius_norm(R);
	// reconstruction_err_.push_back(error);
	return error;
}

int NMQF::record_progress(double rc_l, double relative_e, double time){
	reconstruction_err_.push_back(rc_l);
	relative_err_.push_back(relative_e);
	cpu_time_.push_back(time);
	return 0;
}

int NMQF::fit(gsl_matrix* inputmat, gsl_matrix* U_in, gsl_matrix* V_in, gsl_matrix* W_in, gsl_matrix* H_in, gsl_matrix* p_net, const double lam_1) 
{
	X = inputmat;
	double X_frob = utils::get_frobenius_norm(X);
	int n = X->size1;
	int m = X->size2;
	U = U_in;
	V = V_in;
	W = W_in;
	H = H_in;
	R = gsl_matrix_alloc(n,m);
	A = gsl_matrix_alloc(n,n_components);
	P = gsl_matrix_alloc(n_components,m);
	q = gsl_vector_calloc(m);
	p = gsl_vector_calloc(n);

	if (lam_1>0){
		regA = true;
		Prior = gsl_matrix_alloc(n,n_components);
		gsl_matrix_set_all(Prior, 1);
		gsl_matrix_sub(Prior, p_net);
		gsl_matrix_scale(Prior, lam_1);
	} else {
		regA = false;
	}

	if ((W->size1 != n_components) || (V->size1 != n_components)) {
		cout << "The first dimension of V and W (i.e. their number of rows) should equal the number of components (transcription factors) specified when instantiating NMF." << endl;
		return 1;
	} 
	if ((U->size2 != k1) || (V->size2 != k1)) {
		cout << "The second dimension of U and V (i.e. their number of columns) should equal the number of left components specified when instantiating NMF." << endl;
		return 1;
	} 
	if ((W->size2 != k2) || (H->size2 != k2)) {
		cout << "The second dimension of U and V (i.e. their number of columns) should equal the number of right components specified when instantiating NMF." << endl;
		return 1;
	} 

	if (verbose) {
		cout << "Initializing..." << endl;
	}	
	if (init == random_init) {
		vector<gsl_matrix*> mat_init_vec;
		mat_init_vec.push_back(U);
		mat_init_vec.push_back(V);
		mat_init_vec.push_back(W);
		mat_init_vec.push_back(H);
		init::random(mat_init_vec, random_state);
	} else {
		init::nndsvd(X, U, V, W, H, random_state);
	}
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, U, V, 0, A); // initialize A and P
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, W, H, 0, P); 

	double old_error = calculate_objective();
	printf("initial error: %f\n", old_error);
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
	int tol_it = 10; // TODO decrease expected to be more stable than ibpg
	clock_t c_start = clock();
	for (int n_iter =0; n_iter < max_iter; n_iter++){
		update();
		clock_t c_cur = clock();
		double error = calculate_objective();
		double slope = old_error - error;
		double relativeErr = sqrt(error / X_frob); //utils::get_frobenius_norm get squared norm
		double elapse = (c_cur-c_start) / CLOCKS_PER_SEC;
		record_progress(error, relativeErr, elapse);
		if (verbose) {
			cout << "Itr " << n_iter+1 << " err2 = " << error << ", slope = " << slope 
				<< ", relative = " << relativeErr << ", cpu(s) = " << elapse << endl;
		}
		if ((slope < tol) && (tol_it < 0)) {
			if (verbose) {
				cout << "Converged at iteration " << n_iter+1 << endl;	
			}
			break;
		} else if (slope < tol) {
			--tol_it;
		}
		old_error = error;
		old_slope = slope;
	}
	gsl_matrix_free(R);
	gsl_matrix_free(A);
	gsl_matrix_free(P);
	gsl_vector_free(p);
	gsl_vector_free(q);
	gsl_rng_free(ri);

	/*
	for (int i = 0; i < n_components; i++) {
		for (int j = 0; j < 10; j++) {
			cout << U->data[i * U->tda + j] << ", ";
		}
		cout << endl;
	}
	*/

	return 0;
}
