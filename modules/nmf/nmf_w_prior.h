#ifndef _nmf_w_prior_
#define _nmf_w_prior_
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <list>

#include "basenmf.h"
using namespace std;

class NMFPRIOR : public BASENMF {
	public:
		NMFPRIOR(int, init_method, int, int, bool, double);
		NMFPRIOR(int, init_method, int, int, bool, double, bool);
		~NMFPRIOR();
		int fit(gsl_matrix*, gsl_matrix*, gsl_matrix*, gsl_matrix*, double);
		bool update_permuted;

	private:
		gsl_matrix* U;
		gsl_matrix* V;
		gsl_matrix* R;
		gsl_matrix* Prior; // 1 - actual prior
		int update_kth_block_of_U(int);
		int update_kth_block_of_V(int);
		int update();
		int normalize_and_scale();
		double calculate_objective();
		int record_progress(double, double, double);
		gsl_rng* ri;
};
#endif
