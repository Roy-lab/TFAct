#ifndef _ibpg_2factor_
#define _ibpg_2factor_
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <list>
#include <vector>

#include "basenmf.h"
using namespace std;

class IBPG_2FACTOR : public BASENMF {
	public:
		IBPG_2FACTOR(int, init_method, int, int, bool, double);
		~IBPG_2FACTOR();
		// int setLambdaUxUGraphL(gsl_matrix*, double, double);
		// int setLambdaVxVGraphL(gsl_matrix*, double, double);
		int fit(gsl_matrix*, gsl_matrix*, gsl_matrix*, gsl_matrix*, double, double);
		bool update_permuted;

	private:
		double X_frob; // actually squared frobenius norm
		double cur_relative_err;
		gsl_matrix* R;
		gsl_matrix* A; 
		gsl_matrix* P;
		gsl_matrix* Prior; // lam_1 * (1 - edge confidence)
		double lambda_p;
		
		//  for readability, if a matrix is followed by T, then it means transpose. 
		// 	if y, then it has backslash on its top (extrapolation point 1))

		// update A
		gsl_matrix* Aold;
		gsl_matrix* PPT;
		gsl_matrix* Ay;
		gsl_matrix* gradAy; // To store the proximal gradient
		double L_A_prev;
		bool A_first_it; // govern all

		gsl_matrix* Pold;
		gsl_matrix* ATA;
		gsl_matrix* Py;
		gsl_matrix* gradPy; // To store the proximal gradient
		double L_P_prev;

		// // Update U
		// gsl_matrix* Uold;
		// gsl_matrix* Q; // V^TWH^T
		// gsl_matrix* QQT;
		// gsl_matrix* XQT;
		// gsl_matrix* gradUy;
		// gsl_matrix* Uy; // extrapolate point 1 \.
		// gsl_matrix* PriorV; // soft threshold term
		// gsl_matrix* lambda_LxUy;
		// double L_U_prev;
		// bool U_first_it; // govern all

		// // Update V 
		// gsl_matrix* Vold;
		// gsl_matrix* PPT;
		// gsl_matrix* UTU;
		// gsl_matrix* XTU;
		// gsl_matrix* PXTU;
		// gsl_matrix* gradVy;
		// gsl_matrix* Vy;
		// gsl_matrix* PriorTU; 
		// gsl_matrix* lambda_LxVy;
		// double L_V_prev;
		// // // Update V col
		// // gsl_vector* RTu;
		// // gsl_vector* PRTu; // compute RTu first
		// // gsl_vector* vy;
		// // gsl_vector* gradvy;
		// // double PPT_norm; // ||WH^THW^T||. we note that all L_V has the form 2||ui||^2||WH^THW^T||
		// // vector<double> L_Vk_prev;

		// // Update W
		// gsl_matrix* Wold;
		// gsl_matrix* ATA;
		// gsl_matrix* HTH;
		// gsl_matrix* XH;
		// gsl_matrix* ATXH;
		// gsl_matrix* gradWy;
		// gsl_matrix* Wy;
		// double L_W_prev;
		// // // Update W
		// // gsl_vector* Rh;
		// // gsl_vector* ATRh; // compute Rh first
		// // gsl_vector* wy;
		// // gsl_vector* gradwy;
		// // double ATA_norm; 
		// // vector<double> L_Wk_prev;

		// // Update H
		// gsl_matrix* Hold;
		// gsl_matrix* S; // UV^TW
		// gsl_matrix* STS;
		// gsl_matrix* XTS;
		// gsl_matrix* gradHy;
		// gsl_matrix* Hy; // extrapolate point 1.
		// double L_H_prev;

		double beta;
		double deltaw;
		int inneriterA;
		int inneriterP; 

		double tau;
		double prev_tau;
		double what; //see https://github.com/LeThiKhanhHien/IBPG/blob/master/NMF/IBPG_A.m

		// // regularization
		// bool regA; //l1 on (1-P) \circ (UV^T)
		// bool gregU;
		// bool gregV;
		// gsl_matrix* lambdaUxUGraphL;
		// double lambda_UxUGraphLNorm;
		// gsl_matrix* lambdaVxVGraphL;
		// double lambda_VxVGraphLNorm;

		int update_A();
		int update_P();

		int update();
		// int normalize_and_scale();
		double calculate_objective();
		int record_progress(double, double, double);
		// int graph_init(gsl_matrix*, gsl_matrix*, gsl_matrix*, gsl_matrix*, gsl_matrix*, int);
		gsl_rng* ri;
};
#endif
