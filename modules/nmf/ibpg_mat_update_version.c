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
#include "ibpg.h"
#include "utils.h"

IBPG::IBPG(int k, int k_one, int k_two, init_method initMethod, int maxIter, int seed, bool verb, double termTol, bool perm_up) {
	n_components = k;
	k1 = k_one;
	k2 = k_two;
	init = initMethod;
	max_iter = maxIter;
	random_state = seed;
	verbose = verb;
	tol = termTol;
	list<double> reconstruction_err_;
	update_permuted = perm_up;
	maxInnerDim = max(k, k1);
	delta = 0.01;
	alpha = 0.5; // control # inner updates
	beta = 1.01; // control the influence of the proximal point term
	deltaw = 0.99; // control the gradient term. Also control the proximal point term.
}

IBPG::~IBPG() {};

int IBPG::update_U(){
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, V, P, 0, Q);
	gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1, Q, 0, QQT);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, X, Q, 0, XQT);

	gsl_vector* eval = gsl_vector_alloc(k1); 
	gsl_matrix* tmp_QQT = gsl_matrix_alloc(k1, k1); 
	gsl_matrix_memcpy(tmp_QQT, QQT);
    gsl_eigen_symm_workspace* w = gsl_eigen_symm_alloc(k1); 
    gsl_eigen_symm(tmp_QQT, eval, w); //lower triangular + diagonal destroyed
    gsl_eigen_symm_free(w);
	gsl_matrix_free(tmp_QQT);

	double L_U = max(gsl_vector_max(eval), 0.01);
	gsl_vector_free(eval);

	double gamma;
	if(U_first_it){
		// U_first_it = false; // update after the first update()
		gamma = what; 
	} else {
		gamma = min( what, deltaw*sqrt( L_U_prev/L_U ) ); 
	}
	int j = 0;

	while(j < inneriterU){
		gsl_matrix* tmp_Un = U;
		U = Uold;
		gsl_matrix_sub(U, tmp_Un); //U = -Unold
		gsl_matrix_scale(U, -1*gamma); // U = ww(i)*Unold
		gsl_matrix_memcpy(Uy, tmp_Un);
		gsl_matrix_add(Uy, U);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Uy, QQT, 0, gradUy);
		gsl_matrix_sub(gradUy, XQT); //Note that we have done ||QQT||_2 instead of ||2*QQT||_2
		gsl_matrix_scale(gradUy, 1/L_U);
		gsl_matrix_scale(U, beta);
		gsl_matrix_sub(U, gradUy);
		gsl_matrix_add(U, tmp_Un);
		for(int r = 0; r < U->size1; ++r){
			for(int c = 0; c < k1; ++c){
				double* val =  gsl_matrix_ptr(U, r, c);
				if (*val < 0) {
					*val = 0;
				}
			}
		}
		Uold = tmp_Un;
		++j;
	}
	for(int k = 0; k < k1; ++k){
		gsl_vector_view u_k = gsl_matrix_column(U, k);
		if (gsl_vector_isnull(&u_k.vector)) {
			gsl_vector_set_all(&u_k.vector, 0.01);
		}
	}
	L_U_prev = L_U; // update gradient Lipschitz constant
	return 0;
}

int IBPG::update_H(){
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, W, 0, S);
	gsl_blas_dsyrk(CblasLower, CblasTrans, 1, S, 0, STS);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, X, S, 0, XTS);

	gsl_vector* eval = gsl_vector_alloc(k2); 
	gsl_matrix* tmpSTS = gsl_matrix_alloc(k2, k2); 
	gsl_matrix_memcpy(tmpSTS, STS);
    gsl_eigen_symm_workspace* w = gsl_eigen_symm_alloc(k2); 
    gsl_eigen_symm(tmpSTS, eval, w); //lower triangular + diagonal destroyed
    gsl_eigen_symm_free(w);
	gsl_matrix_free(tmpSTS);

	double L_H = max(gsl_vector_max(eval), 0.01);
	gsl_vector_free(eval);

	double gamma;
	if(U_first_it){
		gamma = what; 
	} else {
		gamma = min( what, deltaw*sqrt( L_H_prev/L_H ) ); 
	}
	int j = 0;

	while(j < inneriterH){
		gsl_matrix* tmp_Hn = H;
		H = Hold;
		gsl_matrix_sub(H, tmp_Hn); //H = -Hnold
		gsl_matrix_scale(H, -1*gamma); // H = wh(i)*Hnold
		gsl_matrix_memcpy(Hy, tmp_Hn);
		gsl_matrix_add(Hy, H);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Hy, STS, 0, gradHy);
		gsl_matrix_sub(gradHy, XTS); //Note that we have done ||STS||_2 instead of ||2*STS||_2
		gsl_matrix_scale(gradHy, 1/L_H);
		gsl_matrix_scale(H, beta);
		gsl_matrix_sub(H, gradHy);
		gsl_matrix_add(H, tmp_Hn);
		for(int r = 0; r < H->size1; ++r){
			for(int c = 0; c < k2; ++c){
				double* val =  gsl_matrix_ptr(H, r, c);
				if (*val < 0) {
					*val = 0;
				}
			}
		}
		Hold = tmp_Hn;
		++j;
	}
	for(int k = 0; k < k2; ++k){
		gsl_vector_view h_k = gsl_matrix_column(H, k);
		if (gsl_vector_isnull(&h_k.vector)) {
			gsl_vector_set_all(&h_k.vector, 0.01);
		}
	}
	L_H_prev = L_H; // update gradient Lipschitz constant
	return 0;
}

int IBPG::update_kth_block_of_V(int k){
	gsl_vector_view u_k = gsl_matrix_column(U, k);
	gsl_vector_view v_k = gsl_matrix_column(V, k);
	gsl_vector_view v_k_old = gsl_matrix_column(Vold, k);
	double u_norm = pow(gsl_blas_dnrm2(&u_k.vector), 2);
	gsl_blas_dgemv(CblasTrans, 1, R, &u_k.vector, 0, RTu);
	gsl_blas_dgemv(CblasNoTrans, 1, P, RTu, 0, PRTu);

	double L_Vk = PPT_norm; // actually 0.5 * L / u_norm^T
	double cur_L_Vk = max((L_Vk * u_norm),0.01);
	double gamma;
	if(U_first_it){
		gamma = what; 
	} else {
		gamma = min( what, deltaw*sqrt( L_Vk_prev.at(k)/cur_L_Vk) ); 
	}
	int j = 0;

	while(j < inneriterV){
		gsl_vector_sub(&v_k_old.vector, &v_k.vector); 
		gsl_vector_scale(&v_k_old.vector, -1*gamma); 
		gsl_vector_memcpy(vy, &v_k.vector);
		gsl_vector_add(vy, &v_k_old.vector);
		gsl_blas_dgemv(CblasNoTrans, 1, PPT, vy, 0, gradvy); 
		gsl_vector_sub(gradvy, PRTu); 
		gsl_vector_scale(gradvy, 1/L_Vk);
		gsl_vector_scale(&v_k_old.vector, beta);
		gsl_vector_sub(&v_k_old.vector, gradvy);
		gsl_vector_add(&v_k_old.vector, &v_k.vector);
		
		gsl_vector_swap(&v_k.vector, &v_k_old.vector);
		for(int r = 0; r < n_components; ++r){
			if (gsl_vector_get(&v_k.vector, r) < 0) {
				gsl_vector_set(&v_k.vector, r, 0);
			}
		}
		++j;
	}

	if (gsl_vector_isnull(&v_k.vector)) {
		gsl_vector_set_all(&v_k.vector, 0.01);
	}
	
	L_Vk_prev.at(k) = cur_L_Vk; // update gradient Lipschitz constant
	return 0;
}

int IBPG::update_kth_block_of_W(int k){
	gsl_vector_view h_k = gsl_matrix_column(H, k);
	gsl_vector_view w_k = gsl_matrix_column(W, k);
	gsl_vector_view w_k_old = gsl_matrix_column(Wold, k);
	double h_norm = pow(gsl_blas_dnrm2(&h_k.vector), 2);
	gsl_blas_dgemv(CblasNoTrans, 1, R, &h_k.vector, 0, Rh);
	gsl_blas_dgemv(CblasTrans, 1, A, Rh, 0, ATRh);

	double L_Wk = ATA_norm; // actually 0.5 * L / u_norm^T
	double cur_L_Wk = max((L_Wk * h_norm), 0.01);
	double gamma;
	if(U_first_it){
		gamma = what; 
	} else {
		gamma = min( what, deltaw*sqrt( L_Wk_prev.at(k)/ cur_L_Wk) ); 
	}
	int j = 0;

	while(j < inneriterW){
		gsl_vector_sub(&w_k_old.vector, &w_k.vector); 
		gsl_vector_scale(&w_k_old.vector, -1*gamma); 
		gsl_vector_memcpy(wy, &w_k.vector);
		gsl_vector_add(wy, &w_k_old.vector);
		gsl_blas_dgemv(CblasNoTrans, 1, ATA, wy, 0, gradwy); 
		gsl_vector_sub(gradwy, ATRh); 
		gsl_vector_scale(gradwy, 1/L_Wk);
		gsl_vector_scale(&w_k_old.vector, beta);
		gsl_vector_sub(&w_k_old.vector, gradwy);
		gsl_vector_add(&w_k_old.vector, &w_k.vector);
		
		gsl_vector_swap(&w_k.vector, &w_k_old.vector);
		for(int r = 0; r < n_components; ++r){
			if (gsl_vector_get(&w_k.vector, r) < 0) {
				gsl_vector_set(&w_k.vector, r, 0);
			}
		}
		++j;
	}

	if (gsl_vector_isnull(&w_k.vector)) {
		gsl_vector_set_all(&w_k.vector, 0.01);
	}
	
	L_Wk_prev.at(k) = cur_L_Wk; // update gradient Lipschitz constant
	return 0;
}

int IBPG::update_kth_block_of_U(int k){
	gsl_vector_view u_k = gsl_matrix_column(U, k);
	// gsl_vector_view v_k = gsl_matrix_column(V, k);
	double q_norm = pow(gsl_blas_dnrm2(q), 2);
	gsl_blas_dgemv(CblasNoTrans, 1/max(q_norm,0.01), R, q, 0, &u_k.vector);
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

int IBPG::update_kth_block_of_H(int k){
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

int IBPG::update_ij_entry_of_V(int i, int j){
	gsl_vector_view u_l = gsl_matrix_column(U, j);
	gsl_vector_view P_k_row = gsl_matrix_row(P, i);
	double* Vij = gsl_matrix_ptr(V, i, j);
	double u_norm = pow(gsl_blas_dnrm2(&u_l.vector), 2);
	double p_norm = pow(gsl_blas_dnrm2(&P_k_row.vector), 2);
	gsl_blas_dgemv(CblasNoTrans, 1/max(p_norm, 0.01), R, &P_k_row.vector, 0, p);
	gsl_blas_ddot(&u_l.vector, p, Vij);
	*Vij /= u_norm;
	if (*Vij < 0) {
		*Vij = 0;
	}
	return 0;
}

int IBPG::update_ij_entry_of_W(int i, int j){
	gsl_vector_view a_k = gsl_matrix_column(A, i);
	gsl_vector_view h_l = gsl_matrix_column(H, j);
	double* Wij = gsl_matrix_ptr(W, i, j);
	double a_norm = pow(gsl_blas_dnrm2(&a_k.vector), 2);
	double h_norm = pow(gsl_blas_dnrm2(&h_l.vector), 2);
	gsl_blas_dgemv(CblasNoTrans, 1/max(a_norm, 0.01), R, &h_l.vector, 0, p);
	gsl_blas_ddot(&a_k.vector, p, Wij);
	*Wij /= h_norm;
	if (*Wij < 0) {
		*Wij = 0;
	}
	return 0;
}

int IBPG::update_w_col() {
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

	// // Update U columns
	// for (int i = 0; i < k1; ++i) {
	// int k = update_order[i];
	// gsl_vector_view u_k = gsl_matrix_column(U, k);
	// gsl_vector_view v_k = gsl_matrix_column(V, k);
	// gsl_blas_dgemv(CblasTrans, 1, P, &v_k.vector, 0, q);
	// gsl_blas_dger(1, &u_k.vector, q, R);
	// update_kth_block_of_U(k);
	// // update_kth_block_of_V(k);
	// gsl_blas_dger(-1, &u_k.vector, q, R);
	// }

	// update whole U
	update_U();
	gsl_matrix_memcpy(R, X);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1, U, Q, 1, R);// Ensure that R is correct
	
	if(update_permuted){
		gsl_ran_shuffle(ri, update_order, k1, sizeof(int));
		gsl_ran_shuffle(ri, update_order_component, n_components, sizeof(int));
	}
	// update columns of V
	gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1, P, 0, PPT);
	{
		gsl_vector* eval = gsl_vector_alloc(n_components); 
		gsl_matrix* tmpPPT = gsl_matrix_alloc(n_components, n_components); 
		gsl_matrix_memcpy(tmpPPT, PPT);
		gsl_eigen_symm_workspace* w = gsl_eigen_symm_alloc(n_components); 
		gsl_eigen_symm(tmpPPT, eval, w); //lower triangular + diagonal destroyed
		gsl_eigen_symm_free(w);
		gsl_matrix_free(tmpPPT);
		PPT_norm = max(gsl_vector_max(eval), 0.01);
		gsl_vector_free(eval);
	}
	for (int j = 0; j < k1; ++j) {
		int l = update_order[j];
		// for (int i = 0; i < n_components; ++i) { //update Vkl
		// 	int k = update_order_component[i];
		// 	gsl_vector_view u_l = gsl_matrix_column(U, l);
		// 	gsl_vector_view P_k_row = gsl_matrix_row(P, k);
		// 	double* Vkl = gsl_matrix_ptr(V, k, l);
		// 	gsl_blas_dger(*Vkl, &u_l.vector, &P_k_row.vector, R);
		// 	update_ij_entry_of_V(k,l);
		// 	gsl_blas_dger(-1*(*Vkl), &u_l.vector, &P_k_row.vector, R);
		// }
		gsl_vector_view v_l = gsl_matrix_column(V, l);
		gsl_vector_view u_l = gsl_matrix_column(U, l);
		gsl_blas_dgemv(CblasTrans, 1, P, &v_l.vector, 0, q);
		gsl_blas_dger(1, &u_l.vector, q, R);
		update_kth_block_of_V(l);
		gsl_blas_dger(-1, &u_l.vector, q, R);
		
		if (gsl_vector_isnull(&v_l.vector)) {// ensure that dividing by V's column norm is defined
			gsl_vector_set_all(&v_l.vector, 0.01);
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
	
	// for (int i = 0; i < k2; ++i) { //update H
	// 	int k = update_order[i];
	// 	gsl_vector_view w_k = gsl_matrix_column(W, k);
	// 	gsl_vector_view h_k = gsl_matrix_column(H, k);
	// 	gsl_blas_dgemv(CblasNoTrans, 1, A, &w_k.vector, 0, p);
	// 	gsl_blas_dger(1, p, &h_k.vector, R);
	// 	// update_kth_block_of_U(k);
	// 	update_kth_block_of_H(k);
	// 	gsl_blas_dger(-1, p, &h_k.vector, R);
	// }

	// Update whole H
	update_H();
	gsl_matrix_memcpy(R, X);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1, S, H, 1, R);

	// update columns of W
	gsl_blas_dsyrk(CblasLower, CblasTrans, 1, A, 0, ATA);
	{
		gsl_vector* eval = gsl_vector_alloc(n_components); 
		gsl_matrix* tmpATA = gsl_matrix_alloc(n_components, n_components); 
		gsl_matrix_memcpy(tmpATA, ATA);
		gsl_eigen_symm_workspace* w = gsl_eigen_symm_alloc(n_components); 
		gsl_eigen_symm(tmpATA, eval, w); //lower triangular + diagonal destroyed
		gsl_eigen_symm_free(w);
		gsl_matrix_free(tmpATA);
		ATA_norm = max(gsl_vector_max(eval), 0.01);
		gsl_vector_free(eval);
	}
	for (int j = 0; j < k2; ++j) {
		int l = update_order[j];
		// for (int i = 0; i < n_components; ++i) { //update Wkl
		// 	int k = update_order_component[i];
		// 	gsl_vector_view a_k = gsl_matrix_column(A, k);
		// 	gsl_vector_view h_l = gsl_matrix_column(H, l);
		// 	double* Wkl = gsl_matrix_ptr(W, k, l);
		// 	gsl_blas_dger(*Wkl, &a_k.vector, &h_l.vector, R);
		// 	update_ij_entry_of_W(k,l);
		// 	gsl_blas_dger(-1*(*Wkl), &a_k.vector, &h_l.vector, R);
		// }
		
		gsl_vector_view w_l = gsl_matrix_column(W, l);
		gsl_vector_view h_l = gsl_matrix_column(H, l);
		gsl_blas_dgemv(CblasNoTrans, 1, A, &w_l.vector, 0, p);
		gsl_blas_dger(1, p, &h_l.vector, R);
		update_kth_block_of_W(l);
		gsl_blas_dger(-1, p, &h_l.vector, R);
		if (gsl_vector_isnull(&w_l.vector)) {// ensure that dividing by W's column norm is defined
			gsl_vector_set_all(&w_l.vector, 0.01);
		}
	}
	
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, W, H, 0, P); // WH are done
    // normalize_and_scale();
	return 0;
}

int IBPG::update() {
	// int update_order[maxInnerDim];
	// int update_order_component[n_components]; // for choosing Vkl, Wkl
	// for (int k = 0; k < k1; ++k) {
	// 	update_order[k] = k;
	// }
	// for (int k = 0; k < n_components; ++k) {
	// 	update_order_component[k] = k;
	// }
	// if(update_permuted){
	// 	gsl_ran_shuffle(ri, update_order, k1, sizeof(int));
	// }

	// update whole U
	update_U();
	gsl_matrix_memcpy(R, X);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1, U, Q, 1, R);// Ensure that R is correct
	
	// if(update_permuted){
	// 	gsl_ran_shuffle(ri, update_order, k1, sizeof(int));
	// 	gsl_ran_shuffle(ri, update_order_component, n_components, sizeof(int));
	// }
	// // update columns of V
	// gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1, P, 0, PPT);
	// {
	// 	gsl_vector* eval = gsl_vector_alloc(n_components); 
	// 	gsl_matrix* tmpPPT = gsl_matrix_alloc(n_components, n_components); 
	// 	gsl_matrix_memcpy(tmpPPT, PPT);
	// 	gsl_eigen_symm_workspace* w = gsl_eigen_symm_alloc(n_components); 
	// 	gsl_eigen_symm(tmpPPT, eval, w); //lower triangular + diagonal destroyed
	// 	gsl_eigen_symm_free(w);
	// 	gsl_matrix_free(tmpPPT);
	// 	PPT_norm = max(gsl_vector_max(eval), 0.01);
	// 	gsl_vector_free(eval);
	// }
	// for (int j = 0; j < k1; ++j) {
	// 	int l = update_order[j];
	// 	gsl_vector_view v_l = gsl_matrix_column(V, l);
	// 	gsl_vector_view u_l = gsl_matrix_column(U, l);
	// 	gsl_blas_dgemv(CblasTrans, 1, P, &v_l.vector, 0, q);
	// 	gsl_blas_dger(1, &u_l.vector, q, R);
	// 	update_kth_block_of_V(l);
	// 	gsl_blas_dger(-1, &u_l.vector, q, R);
		
	// 	if (gsl_vector_isnull(&v_l.vector)) {// ensure that dividing by V's column norm is defined
	// 		gsl_vector_set_all(&v_l.vector, 0.01);
	// 	}
	// }
	
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, U, V, 0, A); // UV are done

	// for (int k = 0; k < k2; ++k) {
	// 	update_order[k] = k;
	// }
	// if(update_permuted){
	// 	gsl_ran_shuffle(ri, update_order, k2, sizeof(int));
	// 	gsl_ran_shuffle(ri, update_order_component, n_components, sizeof(int));
	// }
	
	// Update whole H
	update_H();
	gsl_matrix_memcpy(R, X);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1, S, H, 1, R);

	// // update columns of W
	// gsl_blas_dsyrk(CblasLower, CblasTrans, 1, A, 0, ATA);
	// {
	// 	gsl_vector* eval = gsl_vector_alloc(n_components); 
	// 	gsl_matrix* tmpATA = gsl_matrix_alloc(n_components, n_components); 
	// 	gsl_matrix_memcpy(tmpATA, ATA);
	// 	gsl_eigen_symm_workspace* w = gsl_eigen_symm_alloc(n_components); 
	// 	gsl_eigen_symm(tmpATA, eval, w); //lower triangular + diagonal destroyed
	// 	gsl_eigen_symm_free(w);
	// 	gsl_matrix_free(tmpATA);
	// 	ATA_norm = max(gsl_vector_max(eval), 0.01);
	// 	gsl_vector_free(eval);
	// }
	// for (int j = 0; j < k2; ++j) {
	// 	int l = update_order[j];
	// 	// for (int i = 0; i < n_components; ++i) { //update Wkl
	// 	// 	int k = update_order_component[i];
	// 	// 	gsl_vector_view a_k = gsl_matrix_column(A, k);
	// 	// 	gsl_vector_view h_l = gsl_matrix_column(H, l);
	// 	// 	double* Wkl = gsl_matrix_ptr(W, k, l);
	// 	// 	gsl_blas_dger(*Wkl, &a_k.vector, &h_l.vector, R);
	// 	// 	update_ij_entry_of_W(k,l);
	// 	// 	gsl_blas_dger(-1*(*Wkl), &a_k.vector, &h_l.vector, R);
	// 	// }
		
	// 	gsl_vector_view w_l = gsl_matrix_column(W, l);
	// 	gsl_vector_view h_l = gsl_matrix_column(H, l);
	// 	gsl_blas_dgemv(CblasNoTrans, 1, A, &w_l.vector, 0, p);
	// 	gsl_blas_dger(1, p, &h_l.vector, R);
	// 	update_kth_block_of_W(l);
	// 	gsl_blas_dger(-1, p, &h_l.vector, R);
	// 	if (gsl_vector_isnull(&w_l.vector)) {// ensure that dividing by W's column norm is defined
	// 		gsl_vector_set_all(&w_l.vector, 0.01);
	// 	}
	// }
	
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, W, H, 0, P); // WH are done
    // normalize_and_scale();
	return 0;
}


int IBPG::normalize_and_scale() {
	for (int k = 0; k < k1; ++k) {
		gsl_vector_view u_k = gsl_matrix_column(U, k);
		gsl_vector_view v_k = gsl_matrix_column(V, k);
		double norm = gsl_blas_dnrm2(&u_k.vector);
		gsl_vector_scale(&u_k.vector, 1/norm);
		gsl_vector_scale(&v_k.vector, norm);
	}
	for (int k = 0; k < k2; ++k) {
		gsl_vector_view w_k = gsl_matrix_column(W, k);
		gsl_vector_view h_k = gsl_matrix_column(H, k);
		double norm = gsl_blas_dnrm2(&h_k.vector);
		gsl_vector_scale(&h_k.vector, 1/norm);
		gsl_vector_scale(&w_k.vector, norm);
	}
	return 0;
}

double IBPG::calculate_objective() {
	gsl_matrix_memcpy(R, X);
	// A and P are calculated at the beginning of fit() and in the middle of update()
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1, A, P, 1, R);
	double error = utils::get_frobenius_norm(R);
	reconstruction_err_.push_back(error);
	return error;
}

int IBPG::fit(gsl_matrix* inputmat, gsl_matrix* U_in, gsl_matrix* V_in, gsl_matrix* W_in, gsl_matrix* H_in) {
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
	//For updating U
	Uold = gsl_matrix_alloc(n, k1);
	Q = gsl_matrix_calloc(k1, m);
	QQT = gsl_matrix_calloc(k1, k1);
	XQT = gsl_matrix_calloc(n, k1);
	gradUy = gsl_matrix_calloc(n, k1);
	Uy = gsl_matrix_calloc(n, k1);
	// For V col
	Vold = gsl_matrix_alloc(n_components, k1);
	PPT = gsl_matrix_calloc(n_components, n_components);
	RTu = gsl_vector_calloc(m);
	PRTu = gsl_vector_calloc(n_components);
	gradvy = gsl_vector_calloc(n_components);
	vy = gsl_vector_calloc(n_components);
	for(int i = 0; i < n_components; ++i){
		L_Vk_prev.push_back(-1);
	}
	// For W col
	Wold = gsl_matrix_alloc(n_components, k2);
	ATA = gsl_matrix_calloc(n_components, n_components);
	Rh = gsl_vector_calloc(n);
	ATRh = gsl_vector_calloc(n_components);
	gradwy = gsl_vector_calloc(n_components);
	wy = gsl_vector_calloc(n_components);
	for(int i = 0; i < n_components; ++i){
		L_Wk_prev.push_back(-1);
	}
	//For updating H
	Hold = gsl_matrix_alloc(m, k2);
	S = gsl_matrix_calloc(n, k2);
	STS = gsl_matrix_calloc(k2, k2);
	XTS = gsl_matrix_calloc(m, k2);
	gradHy = gsl_matrix_calloc(m, k2);
	Hy = gsl_matrix_calloc(m, k2);

	// Ensure that input pointers point to the correct objects
	gsl_matrix* Uold_tmp = Uold;
	gsl_matrix* Vold_tmp = Vold;
	gsl_matrix* Wold_tmp = Wold;
	gsl_matrix* Hold_tmp = Hold;

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
	gsl_matrix_memcpy(Uold, U); 
	gsl_matrix_memcpy(Hold, H);
	
	// set up other IBPG parameters
	int count_nzero = 0;
	for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (gsl_matrix_get(X, i, j) != 0) {
                ++count_nzero;
            }
        }
    }
	int totol_entries = n*k1 + k1*n_components + n_components*k2 + k2*m;
	inneriterU = 1 + alpha * (count_nzero - n*k1 + totol_entries) / (n*k1 + n);
	// inneriterV = 1 + alpha * (count_nzero - n_components + totol_entries) / (min(m,n)*k1); //column update
	// inneriterW = 1 + alpha * (count_nzero - n_components + totol_entries) / (min(m,n)*k2); //column 
	inneriterV = 1 + alpha * (count_nzero - n_components + totol_entries) / (n_components*k1 + n_components); //matrix update
	inneriterW = 1 + alpha * (count_nzero - n_components + totol_entries) / (n_components*k2 + n_components); //matrix 
	inneriterH = 1 + alpha * (count_nzero - k2*m + totol_entries) / (k2*m + m);
	tau = 1;
	prev_tau = 1;

	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, U, V, 0, A);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, W, H, 0, P);
	double old_error = calculate_objective();
	printf("initial error: %f\n", old_error);
	double old_slope;
	const gsl_rng_type* T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	ri = gsl_rng_alloc(T);
	gsl_rng_set(ri, random_state);
	int tol_it = 50;
	// Note that the behavior of the first iteration is different. See IBPG_A.m implementation.
	U_first_it = true; 
	clock_t c_start = clock();
	for (int n_iter =0; n_iter < max_iter; n_iter++){
		tau = 0.5 * (1+ sqrt( 1 + 4*pow(prev_tau, 2) ));
		what = (prev_tau-1)/tau;
		prev_tau = tau;
		update();
		clock_t c_cur = clock();
		U_first_it = false;
		double error = calculate_objective();
		double slope = old_error - error;
		double relativeErr = sqrt(error / X_frob); //utils::get_frobenius_norm get squared norm
		double elapse = (c_cur-c_start) / CLOCKS_PER_SEC;
		//
		// if(relativeErr < 2){
		// 	deltaw = relativeErr / 2;
		// } else {
		// 	deltaw = 0.99;
		// }
		//
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
	// Ensure that input pointers point to the correct objects
	if(Uold != Uold_tmp){
		Uold = Uold_tmp;
		gsl_matrix_memcpy(Uold, U);
	}
	if(Vold != Vold_tmp){
		Vold = Vold_tmp;
		gsl_matrix_memcpy(Vold, V);
	}
	if(Wold != Wold_tmp){
		Wold = Wold_tmp;
		gsl_matrix_memcpy(Wold, W);
	}
	if(Hold != Hold_tmp){
		Hold = Hold_tmp;
		gsl_matrix_memcpy(Hold, H);
	}
	gsl_matrix_free(R);
	gsl_matrix_free(A);
	gsl_matrix_free(P);
	gsl_vector_free(p);
	gsl_vector_free(q);
	//U
	gsl_matrix_free(Uold);
	gsl_matrix_free(Q);
	gsl_matrix_free(QQT);
	gsl_matrix_free(XQT);
	gsl_matrix_free(gradUy);
	gsl_matrix_free(Uy);
	//V
	gsl_matrix_free(Vold);
	gsl_matrix_free(PPT);
	gsl_vector_free(RTu);
	gsl_vector_free(PRTu);
	gsl_vector_free(vy);
	gsl_vector_free(gradvy);
	//W
	gsl_matrix_free(Wold);
	gsl_matrix_free(ATA);
	gsl_vector_free(Rh);
	gsl_vector_free(ATRh);
	gsl_vector_free(wy);
	gsl_vector_free(gradwy);
	//H
	gsl_matrix_free(Hold);
	gsl_matrix_free(S);
	gsl_matrix_free(STS);
	gsl_matrix_free(XTS);
	gsl_matrix_free(gradHy);
	gsl_matrix_free(Hy);

	gsl_rng_free(ri);
	return 0;
}
