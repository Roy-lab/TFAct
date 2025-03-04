#ifndef _KL_NMF_2FCT_
#define _KL_NMF_2FCT_
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <list>
#include <vector>
#include <map>

#include "basenmf.h"
using namespace std;

class KL_NMF_2FCT : public BASENMF {
	public:
		KL_NMF_2FCT(int, init_method, int, int, bool, double);
		~KL_NMF_2FCT();
		int fit(map<int, gsl_matrix*>*, map<int, gsl_matrix*>*, map<int, gsl_matrix*>*, map<int, pair<gsl_matrix*, double>>*, double, map<int, map<int, double>>*, double, double, double);

	private:
		double X_KL; // actually the sum of squared frobenius norm of expression matrix
		double cur_relative_err; // treated as concatenating expression matrices to form a large one, and concatenate each X_i - A_iP_i
		map<int, gsl_matrix*>* Xs;
		map<int, gsl_vector*> Xs_rowsums;
		map<int, gsl_vector*> Xs_colsums;

		map<int, gsl_matrix*> Rs; // AP, not X - AP
		map<int, gsl_matrix*>* As;
		map<int, gsl_matrix*>* Ps;
		map<int, gsl_matrix*> Priors; // edge confidence * (1 - edge confidence). Note that this maps may have different size than Rs, As, and Ps
		gsl_matrix* tmpA;
		double lambda_p;
		double lambda_E;
		double lambda_Cminus;
		double epsilon; // KL NMF min val
		map<int, map<int, double>> cell_neighbor;

 
		// For updating among cell type A's terms
		map<int, gsl_matrix*> A_neighbor_sum;
		// remember that when updating A, we have 2\lambda_E\sum... in the denominator. We have a denominator value for each cell type.
		map<int, double> As_normalizing_coeff; 

		map<int, int> inneriterAs;
		map<int, int> inneriterPs; 


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
