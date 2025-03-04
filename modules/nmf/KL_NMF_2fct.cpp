#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_eigen.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <ctime>
#include <omp.h>

#include "initialization.h"
#include "KL_NMF_2fct.h"
#include "utils.h"

KL_NMF_2FCT::KL_NMF_2FCT(int k, init_method initMethod, int maxIter, int seed, bool verb, double termTol){
	n_components = k;
	init = initMethod;
	max_iter = maxIter;
	random_state = seed;
	verbose = verb;
	tol = termTol;
	list<double> reconstruction_err_;
	list<double> relative_err_;
	list<double> cpu_time_;
	epsilon = 1e-15;
}

KL_NMF_2FCT::~KL_NMF_2FCT() {};


int KL_NMF_2FCT::record_progress(double rc_l, double relative_e, double time){
	reconstruction_err_.push_back(rc_l);
	relative_err_.push_back(relative_e);
	cpu_time_.push_back(time);
	return 0;
}

int KL_NMF_2FCT::update_A(int ct){
	gsl_matrix* P = (*Ps)[ct];
	gsl_matrix* X = (*Xs)[ct];
	gsl_matrix* R = Rs[ct];
	
	int inneriterA = inneriterAs[ct];

	bool hasPrior = false;
	map<int, gsl_matrix*>::iterator it=Priors.find(ct);
	if(it != Priors.end()){
		hasPrior = true;
	}

	int j = 0;
	map<int, double> P_rowsum;
	double* data = P->data;
	size_t tda = P->tda;
	for (size_t i = 0; i < n_components; ++i) {
		P_rowsum[i] = 0;
		for (size_t j = 0; j < P->size2; ++j) {
			P_rowsum[i] += data[i * tda + j];
		}
	}

	while(j < inneriterA){
		// gsl_matrix_sub(gradAy, A_neighbor_sum[ct]);
		// gsl_matrix_scale(gradAy, 1/(L_A + As_normalizing_coeff[ct]));
		
        size_t tda = (*As)[ct]->tda;
        double* data = (*As)[ct]->data;

		#pragma omp parallel for collapse(2)
		for(int r = 0; r < ((*As)[ct])->size1; ++r){
			for(int c = 0; c < n_components; ++c){
				double cur_Xrowsum = gsl_vector_get(Xs_rowsums[ct], r);
				double* cur_val =  &data[r * tda + c];
				double b = cur_Xrowsum / *cur_val + P_rowsum[c] - gsl_matrix_get(A_neighbor_sum[ct], r, c);
				for (size_t pc = 0; pc < P->size2; ++pc) { // remember that R = AP
					b -= gsl_matrix_get(X, r, pc) * gsl_matrix_get(P, c, pc) / gsl_matrix_get(R, r, pc);
				}
				if(hasPrior){ // branch prediction is not bad here
					b += gsl_matrix_get(Priors[ct], r, c);
				} else {
					b += lambda_Cminus;
				}
				if (lambda_E > 0){
					*cur_val = (-b + sqrt(b*b + 4 * As_normalizing_coeff[ct] * cur_Xrowsum))/(2 * As_normalizing_coeff[ct]);
				} else {
					*cur_val = cur_Xrowsum / b;
				}
				
				if (*cur_val < epsilon) {
					*cur_val = epsilon;
				}
			}
		}
		++j;
		// update R
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, (*As)[ct], P, 0, R);
	}

    return 0;
}
		
int KL_NMF_2FCT::update_P(int ct){
	gsl_matrix* A = (*As)[ct];
	gsl_matrix* X = (*Xs)[ct];
	gsl_matrix* R = Rs[ct];
	int inneriterP = inneriterPs[ct];

	map<int, double> A_colsum;
	double* data = A->data;
	size_t tda = A->tda;
	for (size_t j = 0; j < n_components; ++j) {
		A_colsum[j] = 0;
		for (size_t i = 0; i < A->size1; ++i) {
			A_colsum[j] += data[i * tda + j];
		}
	}

	int j = 0;

	while(j < inneriterP){
	    size_t tda = ((*Ps)[ct])->tda;
	    double* data = ((*Ps)[ct])->data;
		
		#pragma omp parallel for collapse(2)
		for(int c = 0; c < ((*Ps)[ct])->size2; ++c){
			for(int r = 0; r < n_components; ++r){
				double cur_Xcolsum = gsl_vector_get(Xs_colsums[ct], c);
				double* cur_val =  &data[r * tda + c];
				double b = A_colsum[r];
				for (size_t ar = 0; ar < A->size1; ++ar) { // remember that R = AP
					b -= gsl_matrix_get(X, ar, c) * gsl_matrix_get(A, ar, r) / gsl_matrix_get(R, ar, c);
				}
				b = (b+lambda_p) * (*cur_val);
				*cur_val /= 1 + b/cur_Xcolsum;
				if (*cur_val < epsilon) {
					*cur_val = epsilon;
				}
			}
		}
		++j;
		// update R; don't need at the last iteration as calculate_objective will do its job.
		if(j < inneriterP){
			gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, (*Ps)[ct], 0, R);
		}
	}

    return 0;
}

int KL_NMF_2FCT::subtract_A_from_neighbor_sum(int ct){
	map<int, double>::iterator it;
	for(it = (cell_neighbor[ct]).begin(); it != (cell_neighbor[ct]).end(); ++it){
		int cur_neighbor = it->first;
		double cur_weight = it->second;
		gsl_matrix_memcpy(tmpA, (*As)[ct]);
		gsl_matrix_scale(tmpA, cur_weight);
		gsl_matrix_sub(A_neighbor_sum[cur_neighbor], tmpA);
	}
	return 0;
}

int KL_NMF_2FCT::add_A_to_neighbor_sum(int ct){
	map<int, double>::iterator it;
	for(it = (cell_neighbor[ct]).begin(); it != (cell_neighbor[ct]).end(); ++it){
		int cur_neighbor = it->first;
		double cur_weight = it->second;
		gsl_matrix_memcpy(tmpA, (*As)[ct]);
		gsl_matrix_scale(tmpA, cur_weight);
		gsl_matrix_add(A_neighbor_sum[cur_neighbor], tmpA);
	}
	return 0;
}


int KL_NMF_2FCT::update(){
	map<int, gsl_matrix*>::iterator it;
	for(it = Xs->begin(); it != Xs->end(); ++it){
		subtract_A_from_neighbor_sum(it->first);
		update_A(it->first);
		add_A_to_neighbor_sum(it->first);
    	update_P(it->first);
	}
    return 0;
}

double KL_NMF_2FCT::calculate_objective(){
	cur_relative_err = 0;
	double sum = 0;
	map<int, gsl_matrix*>::iterator it;
	for(it = Xs->begin(); it != Xs->end(); ++it){
		sum += calculate_objective(it->first);
	}
	map<int, map<int, double>>::iterator it2;
	map<int, double>::iterator neighbors_it;
	for(it2 = cell_neighbor.begin(); it2 != cell_neighbor.end(); ++it2){
		map<int, double>* neighbor_lst = &(it2->second);
		int ct = it2->first;
		for(neighbors_it = neighbor_lst->begin(); neighbors_it != neighbor_lst->end(); ++neighbors_it){
			int ct_cur_neighbor = neighbors_it->first;
			if(ct < ct_cur_neighbor){
				sum += calculate_among_cell_loss((*As)[ct], (*As)[ct_cur_neighbor], neighbors_it->second);
			}
		}
	}
	return sum;
}

// cell 1 A, cell 2 A, weights proportional to similarity (note that lambda_E is already incorporated)
double KL_NMF_2FCT::calculate_among_cell_loss(gsl_matrix* A1, gsl_matrix* A2, double weight){
	gsl_matrix_memcpy(tmpA, A1);
	gsl_matrix_sub(tmpA, A2);
	return utils::get_frobenius_norm(tmpA) * 0.5 * weight;
}

double KL_NMF_2FCT::calculate_objective(int ct){
	gsl_matrix* X = (*Xs)[ct]; 
	gsl_matrix* R = Rs[ct]; // R is AP, not X - AP this time
	gsl_matrix* A = (*As)[ct]; 
	gsl_matrix* P = (*Ps)[ct];
    // gsl_matrix_memcpy(R, X);
	// A and P are calculated at the beginning of fit() and in the middle of update()
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, P, 0, R);
	double sum = utils::get_KL(X, R);
    cur_relative_err += sum;
    // sum *= 0.5;
    size_t rows = A->size1;
    double* data = A->data; // Pointer to the raw data array
	size_t tda = A->tda;
	map<int,gsl_matrix*>::iterator it = Priors.find(ct);
	if(it != Priors.end()){
		gsl_matrix* Prior = Priors[ct];
		double* pdata = Prior->data;
		for (size_t i = 0; i < rows; ++i) {
			for (size_t j = 0; j < n_components; ++j) {
				sum += data[i * tda + j] * pdata[i * tda + j];
			}
		}
	} else {
		for (size_t i = 0; i < rows; ++i) {
			for (size_t j = 0; j < n_components; ++j) {
				sum += data[i * tda + j] * lambda_Cminus;
			}
		}
	}
    data = P->data;
    tda = P->tda;
    size_t cols = P->size2;
    for (size_t i = 0; i < n_components; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            sum += data[i * tda + j] * lambda_p;
        }
    }
    return sum;
}

// reminder: lambda in the Framework.C is lambda_C
int KL_NMF_2FCT::fit(map<int, gsl_matrix*>* inputmat, map<int, gsl_matrix*>* A_in, map<int, gsl_matrix*>* P_in, 
				   map<int, pair<gsl_matrix*, double>>* p_net, double lam_p, 
				   map<int, map<int, double>>* cell_neighborhood, double lam_E, double lam_Cminus, double lambda1){
    Xs = inputmat;
	// X_KL = 0;
	int n = ((Xs->begin())->second)->size1;
	map<int, int> ms;
	// int m = X->size2;
	map<int, gsl_matrix*> A_in_copy = *A_in;
	map<int, gsl_matrix*> P_in_copy = *P_in;
	As = &A_in_copy;
	Ps = &P_in_copy;
	tmpA = gsl_matrix_alloc(n, n_components);
	if (verbose) {
		cout << "Initializing..." << endl;
	}

	{
		map<int, gsl_matrix*>::iterator it;
		for(it = Xs->begin(); it != Xs->end(); ++it){
			int ct = it->first;
			ms[ct] = (it->second)->size2;
			// X_KL += utils::get_KL(it->second);
			gsl_matrix* A = (*As)[ct];
			gsl_matrix* P = (*Ps)[ct];
			if ((A->size2 != n_components) || (P->size1 != n_components)) {
				cout << "The inner dimension of the two factor matrices should agree with the supplied number of factors (regulators)." << endl;
				return 1;
			}
			// weights initialization
			if (init == random_init) {
				init::random(A, P, random_state);
			} else {
				gsl_matrix* H = gsl_matrix_calloc(ms[ct], n_components);
				init::nndsvd2((*Xs)[ct], A, H, random_state);
				gsl_matrix_transpose_memcpy(P, H);
				gsl_matrix_free(H);
			}
			double* data = A->data;
			size_t tda = A->tda; 
			for (size_t i = 0; i < n; ++i) {
				for (size_t j = 0; j < n_components; ++j) {
					double* value = &data[i * tda + j];
					if (*value < epsilon) {
						*value = epsilon;
					}
				}
			}
			data = P->data;
			tda = P->tda; 
			for (size_t i = 0; i < n_components; ++i) {
				for (size_t j = 0; j < ms[ct]; ++j) {
					double* value = &data[i * tda + j];
					if (*value < epsilon) {
						*value = epsilon;
					}
				}
			}

			Rs[ct] = gsl_matrix_alloc(n, ms[ct]);

			//For updating A
			Xs_rowsums[ct] = gsl_vector_calloc(n);

			//For updating P
			Xs_colsums[ct] = gsl_vector_calloc(ms[ct]);

			data = (it->second)->data;
			tda = (it->second)->tda;
			for (size_t i = 0; i < n; ++i) {
				for (size_t j = 0; j < ms[ct]; ++j) {
					double value = data[i * tda + j];
					gsl_vector_set(Xs_rowsums[ct], i, gsl_vector_get(Xs_rowsums[ct], i) + value);
					gsl_vector_set(Xs_colsums[ct], j, gsl_vector_get(Xs_colsums[ct], j) + value);
				}
			}
			
			// For updating among cell type A's terms
			A_neighbor_sum[ct] = gsl_matrix_calloc(n, n_components);

			// set up prior penalties
			map<int, pair<gsl_matrix*, double>>::iterator it_neighbor = (*p_net).find(ct);

    		if (it_neighbor != (*p_net).end()) {
				Priors[ct] = gsl_matrix_alloc(n, n_components);
				gsl_matrix_memcpy(Priors[ct], ((*p_net)[ct]).first); // note that p_net has entries of edge confidence ranged from 0 to 1
				gsl_matrix_add_constant(Priors[ct], -1); 
				gsl_matrix_scale(Priors[ct], -1 * lambda1 * ((*p_net)[ct]).second); // make edge confidence to penalty terms. Edge not appearing or with low confidence will be penalized
    		}

		}
	}

	/*
	* Setting up hyperparameters
	*/
    lambda_p = lam_p;
	lambda_E = lam_E;
	lambda_Cminus = lam_Cminus;
	// TODO change this
	{
		map<int, gsl_matrix*>::iterator it;
		for(it = Xs->begin(); it != Xs->end(); ++it) {
			inneriterAs[it->first]=1;
			inneriterPs[it->first]=1;
		}
	}

	// set up edge similarity (assume that the user supplies distances). Note that lambda_E is included
	cell_neighbor = *cell_neighborhood;
	{
		map<int, map<int, double>>::iterator it;
		map<int, double>::iterator it2;
		for(it = cell_neighbor.begin(); it != cell_neighbor.end(); ++it){
			for(it2 = (it->second).begin(); it2 != (it->second).end(); ++it2){
				it2->second = lambda_E / (1 + it2->second);
			}
		}
	}

	// now initializing among-cell type terms
	{
		map<int, gsl_matrix*>::iterator it;
		map<int, double>::iterator it2;
		for(it = Xs->begin(); it != Xs->end(); ++it){
			int ct = it->first;
			size_t tda = (A_neighbor_sum[ct])->tda;
        	double* data = (A_neighbor_sum[ct])->data;
			int rows = (A_neighbor_sum[ct])->size1;
			int cols = (A_neighbor_sum[ct])->size2;
			for(it2 = (cell_neighbor[ct]).begin(); it2 != (cell_neighbor[ct]).end(); ++it2){
				int cur_neighbor = it2->first;
				double cur_weight = it2->second;
				double* neighbor_dat = ((*As)[cur_neighbor])->data;
				
				for(int i = 0; i < rows; ++i){
					for(int j = 0; j < cols; ++j){
						data[i * tda + j] += neighbor_dat[i * tda + j] * cur_weight;
					}
				}
				map<int, double>::iterator As_normalize_it = As_normalizing_coeff.find(ct);
				if(As_normalize_it != As_normalizing_coeff.end()){
					As_normalizing_coeff[ct] += cur_weight;
				} else {
					As_normalizing_coeff[ct] = cur_weight;
				}
				
			}
		}
	}
	
	clock_t c_start = clock();

    int below_tol_it = 10; // change this
	int below_tol_it_counter = 1;

    double old_error = calculate_objective(); // this also populates R with AP
    double old_slope;
    record_progress(old_error, cur_relative_err, 0);
    for (int n_iter =0; n_iter < max_iter; ++n_iter){
		update();
		clock_t c_cur = clock();
		double error = calculate_objective();
		double slope = old_error - error;
		double elapse = (c_cur-c_start) / CLOCKS_PER_SEC;
		record_progress(error, cur_relative_err, elapse);
		if (verbose) {
			cout << "Itr " << n_iter+1 << " err2 = " << error << ", slope = " << slope 
				<< ", divergence = " << cur_relative_err << ", cpu(s) = " << elapse << endl;
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

	// freeing stuffs and postprocessing
	{
		map<int, gsl_matrix*>::iterator it;
		for(it = Xs->begin(); it != Xs->end(); ++it){
			int ct = it->first;
			gsl_matrix_free(Rs[ct]);
			gsl_matrix_free(Priors[ct]);
			// Free A
			gsl_vector_free(Xs_rowsums[ct]);

			// // Free P
			gsl_vector_free(Xs_colsums[ct]);


			// freeing among cell type terms
			gsl_matrix_free(A_neighbor_sum[ct]);
		}
	}
	gsl_matrix_free(tmpA);
    return 0;
}