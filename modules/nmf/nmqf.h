#ifndef _nmqf_
#define _nmqf_
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <list>

#include "basenmf.h"
using namespace std;

class NMQF : public BASENMF {
	public:
		NMQF(int, init_method, int, int, bool, double);
		NMQF(int, int, int, init_method, int, int, bool, double, bool);
		~NMQF();
		int fit(gsl_matrix*, gsl_matrix*, gsl_matrix*, gsl_matrix*, gsl_matrix*, gsl_matrix*, double);
		bool update_permuted;
		double get_absolute_sum(const gsl_matrix*);

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
		int k1;
		int k2;
		int maxInnerDim;
		// regularization
		bool regA; //l1 on (1-P) \circ (UV^T)
		gsl_matrix* Prior; // 1 - actual prior, where edges in the actual prior range from 0 to 1
		int update_kth_block_of_U(int);
		// int update_kth_block_of_V(int);
		int update_kth_block_of_H(int);
		int update_ij_entry_of_V(int, int);
		int update_ij_entry_of_W(int, int);
		int update();
		int normalize_and_scale();
		double calculate_objective();
		int record_progress(double, double, double);
		gsl_rng* ri;
};
#endif
