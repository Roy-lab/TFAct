#ifndef _ibpg_
#define _ibpg_
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <list>
#include <vector>

#include "basenmf.h"
using namespace std;

class IBPG : public BASENMF {
	public:
		IBPG(int, int, int, init_method, int, int, bool, double, bool);
		~IBPG();
		int fit(gsl_matrix*, gsl_matrix*, gsl_matrix*, gsl_matrix*, gsl_matrix*);
		bool update_permuted;

	private:
		gsl_matrix* A; //UVT
		gsl_matrix* P; //WHT
		// gsl_matrix* VTWHT;
		// gsl_matrix* UVTW;
		gsl_matrix* U;
		gsl_matrix* V;
		gsl_matrix* W;
		gsl_matrix* H;
		gsl_matrix* R;
		
		gsl_vector* q; //rows of matrix of dimension VTWHT
		gsl_vector* p; //columns of matrix of dimension UVTW. Not P.

		// Update U
		gsl_matrix* Uold;
		gsl_matrix* Q; // V^TWH^T
		gsl_matrix* QQT;
		gsl_matrix* XQT;
		gsl_matrix* gradUy;
		gsl_matrix* Uy; // extrapolate point 1.
		double L_U_prev;
		bool U_first_it; // govern all

		// Update V 
		gsl_matrix* Vold;
		gsl_matrix* PPT;
		gsl_matrix* UTU;
		double L_V_prev;
		// // Update V col
		// gsl_vector* RTu;
		// gsl_vector* PRTu; // compute RTu first
		// gsl_vector* vy;
		// gsl_vector* gradvy;
		// double PPT_norm; // ||WH^THW^T||. we note that all L_V has the form 2||ui||^2||WH^THW^T||
		// vector<double> L_Vk_prev;

		// Update W
		gsl_matrix* Wold;
		gsl_matrix* ATA;
		gsl_matrix* HTH;
		double L_W_prev;
		// // Update W
		// gsl_vector* Rh;
		// gsl_vector* ATRh; // compute Rh first
		// gsl_vector* wy;
		// gsl_vector* gradwy;
		// double ATA_norm; 
		// vector<double> L_Wk_prev;

		// Update H
		gsl_matrix* Hold;
		gsl_matrix* S; // UV^TW
		gsl_matrix* STS;
		gsl_matrix* XTS;
		gsl_matrix* gradHy;
		gsl_matrix* Hy; // extrapolate point 1.
		double L_H_prev;

		int k1;
		int k2;
		int maxInnerDim;

		double delta;
		double alpha;
		double beta;
		double deltaw;
		int inneriterU;
		int inneriterV; //column update
		int inneriterW; //column 
		int inneriterH;
		double tau;
		double prev_tau;
		double what; //see https://github.com/LeThiKhanhHien/IBPG/blob/master/NMF/IBPG_A.m

		int update_kth_block_of_U(int);
		int update_U();
		int update_H();
		int update_kth_block_of_V(int);
		int update_kth_block_of_W(int);
		int update_kth_block_of_H(int);
		int update_ij_entry_of_V(int, int);
		int update_ij_entry_of_W(int, int);
		int update();
		int update_w_col();
		int normalize_and_scale();
		double calculate_objective();
		gsl_rng* ri;
};
#endif
