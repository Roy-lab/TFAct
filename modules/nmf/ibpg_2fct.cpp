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
#include "ibpg_2fct.h"
#include "utils.h"

IBPG_2FCT::IBPG_2FCT(int k, init_method initMethod, int maxIter, int seed, bool verb, double termTol){
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
	epsilon = 1e-3; // pre-conditioning 
}

IBPG_2FCT::~IBPG_2FCT() {};


int IBPG_2FCT::record_progress(double rc_l, double relative_e, double time){
	reconstruction_err_.push_back(rc_l);
	relative_err_.push_back(relative_e);
	cpu_time_.push_back(time);
	return 0;
}

int IBPG_2FCT::update_A(int ct){
	gsl_matrix* P = (*Ps)[ct];
	gsl_matrix* X = (*Xs)[ct];
	// gsl_matrix* Aold = Aolds[ct];
	gsl_matrix* PPT = PPTs[ct];
	gsl_matrix* gradAy = gradAys[ct];
	gsl_matrix* Ay = Ays[ct];
	double& L_A_prev = L_A_prevs[ct];
	int inneriterA = inneriterAs[ct];

	bool hasPrior = false;
	map<int, gsl_matrix*>::iterator it=Priors.find(ct);
	if(it != Priors.end()){
		hasPrior = true;
	}

	gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1, P, 0, PPT);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, X, P, 0, gradAy); 

	gsl_vector* eval = gsl_vector_alloc(n_components); 
	gsl_matrix* tmp_PPT = gsl_matrix_alloc(n_components, n_components); 
	gsl_matrix_memcpy(tmp_PPT, PPT);
    gsl_eigen_symm_workspace* w = gsl_eigen_symm_alloc(n_components); 
	gsl_vector_view diag = gsl_matrix_diagonal(tmp_PPT); // preconditioning. The QR step just halts.
    gsl_vector_add_constant(&diag.vector, epsilon);
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
		gsl_matrix* tmp_An = (*As)[ct];
		(*As)[ct] = Aolds[ct];
		gsl_matrix_sub((*As)[ct], tmp_An); // -1 * change in A
		gsl_matrix_scale((*As)[ct], -1*gamma); 
		gsl_matrix_memcpy(Ay, tmp_An);
		gsl_matrix_add(Ay, (*As)[ct]);
        gsl_blas_dsymm(CblasRight, CblasLower, 1, PPT, Ay, -1, gradAy); 
		if(hasPrior){
			gsl_matrix_add(gradAy, Priors[ct]);
		} else {
			gsl_matrix_add_constant(gradAy, lambda_Cminus);
		}
		gsl_matrix_sub(gradAy, A_neighbor_sum[ct]);

		gsl_matrix_scale(gradAy, 1/(L_A + As_normalizing_coeff[ct]));

		gsl_matrix_scale((*As)[ct], beta);
		gsl_matrix_add((*As)[ct], tmp_An);
		gsl_matrix_scale((*As)[ct], L_A/(L_A + As_normalizing_coeff[ct]));
		gsl_matrix_sub((*As)[ct], gradAy);
		
        size_t tda = (*As)[ct]->tda;
        double* data = (*As)[ct]->data;
		for(int r = 0; r < ((*As)[ct])->size1; ++r){
			for(int c = 0; c < n_components; ++c){
				double* val =  &data[r * tda + c];
				if (*val < 0) {
					*val = 0;
				}
			}
		}
		Aolds[ct] = tmp_An;
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
		
int IBPG_2FCT::update_P(int ct){
	gsl_matrix* A = (*As)[ct];
	// gsl_matrix* P = (*Ps)[ct];
	gsl_matrix* X = (*Xs)[ct];
	// gsl_matrix* Pold = Polds[ct];
	gsl_matrix* ATA = ATAs[ct];
	gsl_matrix* gradPy = gradPys[ct];
	gsl_matrix* Py = Pys[ct];
	double& L_P_prev = L_P_prevs[ct];
	int inneriterP = inneriterPs[ct];

	gsl_blas_dsyrk(CblasLower, CblasTrans, 1, A, 0, ATA);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, A, X, 0, gradPy); 

	gsl_vector* eval = gsl_vector_alloc(n_components); 
	gsl_matrix* tmp_ATA = gsl_matrix_alloc(n_components, n_components); 
	gsl_matrix_memcpy(tmp_ATA, ATA);
	gsl_eigen_symm_workspace* w = gsl_eigen_symm_alloc(n_components); 
	gsl_vector_view diag = gsl_matrix_diagonal(tmp_ATA); // preconditioning.
    gsl_vector_add_constant(&diag.vector, epsilon);
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
	    gsl_matrix* tmp_Pn = (*Ps)[ct];
	    (*Ps)[ct] = Polds[ct];
	    gsl_matrix_sub((*Ps)[ct], tmp_Pn); // -1 * change in A
	    gsl_matrix_scale((*Ps)[ct], -1*gamma); 
	    gsl_matrix_memcpy(Py, tmp_Pn);
	    gsl_matrix_add(Py, (*Ps)[ct]);
	    gsl_blas_dsymm(CblasLeft, CblasLower, 1, ATA, Py, -1, gradPy); 
	    gsl_matrix_add_constant(gradPy, lambda_p);

	    gsl_matrix_scale(gradPy, 1/L_P);
	    gsl_matrix_scale((*Ps)[ct], beta);
	    gsl_matrix_sub((*Ps)[ct], gradPy);
	    gsl_matrix_add((*Ps)[ct], tmp_Pn);

	    size_t tda = ((*Ps)[ct])->tda;
	    double* data = ((*Ps)[ct])->data;
	    for(int r = 0; r < n_components; ++r){
		    for(int c = 0; c < ((*Ps)[ct])->size2; ++c){
			    double* val =  &data[r * tda + c];
			    if (*val < 0) {
				    *val = 0;
			    }
		    }
	    }
	    Polds[ct] = tmp_Pn;
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

int IBPG_2FCT::subtract_A_from_neighbor_sum(int ct){
	gsl_matrix* tmpA = (gradAys.begin())->second; // as this matrix is only used in the update, we can treat it as a tmp
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

int IBPG_2FCT::add_A_to_neighbor_sum(int ct){
	gsl_matrix* tmpA = (gradAys.begin())->second; // as this matrix is only used in the update, we can treat it as a tmp
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


int IBPG_2FCT::update(){
	map<int, gsl_matrix*>::iterator it;
	for(it = Xs->begin(); it != Xs->end(); ++it){
		subtract_A_from_neighbor_sum(it->first);
		update_A(it->first);
		add_A_to_neighbor_sum(it->first);
    	update_P(it->first);
	}
    return 0;
}

double IBPG_2FCT::calculate_objective(){
	cur_relative_err = 0;
	double sum = 0;
	map<int, gsl_matrix*>::iterator it;
	for(it = Xs->begin(); it != Xs->end(); ++it){
		sum += calculate_objective(it->first);
	}
	cur_relative_err = sqrt(cur_relative_err/X_frob);
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
double IBPG_2FCT::calculate_among_cell_loss(gsl_matrix* A1, gsl_matrix* A2, double weight){
	gsl_matrix* tmpA = (gradAys.begin())->second; // as this matrix is only used in the update, we can treat it as a tmp
	gsl_matrix_memcpy(tmpA, A1);
	gsl_matrix_sub(tmpA, A2);
	return utils::get_frobenius_norm(tmpA) * 0.5 * weight;
}

double IBPG_2FCT::calculate_objective(int ct){
	gsl_matrix* X = (*Xs)[ct]; 
	gsl_matrix* R = Rs[ct];
	gsl_matrix* A = (*As)[ct]; 
	gsl_matrix* P = (*Ps)[ct];
    gsl_matrix_memcpy(R, X);
	// A and P are calculated at the beginning of fit() and in the middle of update()
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1, A, P, 1, R);
	double sum = utils::get_frobenius_norm(R);
    cur_relative_err += sum;
    sum *= 0.5;
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
int IBPG_2FCT::fit(map<int, gsl_matrix*>* inputmat, map<int, gsl_matrix*>* A_in, map<int, gsl_matrix*>* P_in, 
				   map<int, pair<gsl_matrix*, double>>* p_net, double lam_p, 
				   map<int, map<int, double>>* cell_neighborhood, double lam_E, double lam_Cminus, double lambda1){
    Xs = inputmat;
	X_frob = 0;
	int n = ((Xs->begin())->second)->size1;
	map<int, int> ms;
	// int m = X->size2;
	map<int, gsl_matrix*> A_in_copy = *A_in;
	map<int, gsl_matrix*> P_in_copy = *P_in;
	As = &A_in_copy;
	Ps = &P_in_copy;
	if (verbose) {
		cout << "Initializing..." << endl;
	}
	map<int, gsl_matrix*> Aold_tmps;
    map<int, gsl_matrix*> Pold_tmps;
	{
		map<int, gsl_matrix*>::iterator it;
		for(it = Xs->begin(); it != Xs->end(); ++it){
			int ct = it->first;
			ms[ct] = (it->second)->size2;
			X_frob += utils::get_frobenius_norm(it->second);
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
			Rs[ct] = gsl_matrix_alloc(n, ms[ct]);

			//For updating A
			Aolds[ct] = gsl_matrix_alloc(n, n_components);
			gsl_matrix_memcpy(Aolds[ct], A); 
			PPTs[ct] = gsl_matrix_calloc(n_components, n_components);
			gradAys[ct] = gsl_matrix_calloc(n, n_components);
			Ays[ct] = gsl_matrix_calloc(n, n_components);

			//For updating P
			Polds[ct] = gsl_matrix_alloc(n_components, ms[ct]);
			gsl_matrix_memcpy(Polds[ct], P);
			ATAs[ct] = gsl_matrix_calloc(n_components, n_components);
			gradPys[ct] = gsl_matrix_calloc(n_components, ms[ct]);
			Pys[ct] = gsl_matrix_calloc(n_components, ms[ct]);
			
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
			
			// Below ensures that the right matrices are pointed to by the input pointers, not discarded.
    		Aold_tmps[ct] = Aolds[ct];
    		Pold_tmps[ct] = Polds[ct]; 
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
    // inneriterA = 1;
    // inneriterP = 1;
    tau = 1;
	prev_tau = 1;

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
				const double* neighbor_dat = ((*As)[cur_neighbor])->data;
				
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

	// freeing stuffs and postprocessing
	{
		map<int, gsl_matrix*>::iterator it;
		for(it = Xs->begin(); it != Xs->end(); ++it){
			int ct = it->first;
			// Ensure that input pointers point to the correct objects
			if(Aolds[ct] != Aold_tmps[ct]){
				gsl_matrix_memcpy(Aolds[ct], (*As)[ct]);
				Aolds[ct] = Aold_tmps[ct];
			}
			if(Polds[ct] != Pold_tmps[ct]){
				gsl_matrix_memcpy(Polds[ct], (*Ps)[ct]);
				Polds[ct] = Pold_tmps[ct];
			}
			gsl_matrix_free(Rs[ct]);
			gsl_matrix_free(Priors[ct]);
			// Free A
			gsl_matrix_free(Aolds[ct]);
			gsl_matrix_free(PPTs[ct]);
			gsl_matrix_free(gradAys[ct]);
			gsl_matrix_free(Ays[ct]);
			// Free P
			gsl_matrix_free(Polds[ct]);
			gsl_matrix_free(ATAs[ct]);
			gsl_matrix_free(gradPys[ct]);
			gsl_matrix_free(Pys[ct]);

			// freeing among cell type terms
			gsl_matrix_free(A_neighbor_sum[ct]);
		}
	}

    gsl_rng_free(ri);
    return 0;
}
