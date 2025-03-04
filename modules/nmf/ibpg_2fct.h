#ifndef _ibpg_2fct_
#define _ibpg_2fct_
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <list>
#include <vector>
#include <map>

#include "basenmf.h"
using namespace std;

class IBPG_2FCT : public BASENMF {
	public:
		IBPG_2FCT(int, init_method, int, int, bool, double);
		~IBPG_2FCT();
		int fit(map<int, gsl_matrix*>*, map<int, gsl_matrix*>*, map<int, gsl_matrix*>*, map<int, pair<gsl_matrix*, double>>*, double, map<int, map<int, double>>*, double, double, double);
		bool update_permuted;

	private:
		double X_frob; // actually the sum of squared frobenius norm of expression matrix
		double cur_relative_err; // treated as concatenating expression matrices to form a large one, and concatenate each X_i - A_iP_i
		map<int, gsl_matrix*>* Xs;
		map<int, gsl_matrix*> Rs;
		map<int, gsl_matrix*>* As;
		map<int, gsl_matrix*>* Ps;
		map<int, gsl_matrix*> Priors; // edge confidence * (1 - edge confidence). Note that this maps may have different size than Rs, As, and Ps
		double lambda_p;
		double lambda_E;
		double lambda_Cminus;
		double epsilon;
		map<int, map<int, double>> cell_neighbor;

		//  for readability, if a matrix is followed by T, then it means transpose. 
		// 	if y, then it has backslash on its top (extrapolation point 1))

		// update A
		map<int, gsl_matrix*> Aolds;
		map<int, gsl_matrix*> PPTs;
		map<int, gsl_matrix*> Ays;
		map<int, gsl_matrix*> gradAys; // To store the proximal gradient
		map<int, double> L_A_prevs;
		bool A_first_it; // govern all

		// update P
		map<int, gsl_matrix*> Polds;
		map<int, gsl_matrix*> ATAs;
		map<int, gsl_matrix*> Pys;
		map<int, gsl_matrix*> gradPys; // To store the proximal gradient
		map<int, double> L_P_prevs;
 
		// For updating among cell type A's terms
		map<int, gsl_matrix*> A_neighbor_sum;
		// remember that when updating A, we divide it by the function's gradient's lipschitz constant plus the sum of its incident edge weights
		map<int, double> As_normalizing_coeff; 
		
		double beta;
		double deltaw;
		map<int, int> inneriterAs;
		map<int, int> inneriterPs; 

		double tau;
		double prev_tau;
		double what; //see https://github.com/LeThiKhanhHien/IBPG/blob/master/NMF/IBPG_A.m


		int update_A(int); // int is cell type
		int update_P(int);

		int update();
		double calculate_objective();
		double calculate_objective(int);
		int subtract_A_from_neighbor_sum(int);
		int add_A_to_neighbor_sum(int);

		int record_progress(double, double, double);
		double calculate_among_cell_loss(gsl_matrix*, gsl_matrix*, double);
		gsl_rng* ri;
};
#endif
