#ifndef _nmf_wo_init_
#define _nmf_wo_init_
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <list>

#include "basenmf.h"
using namespace std;

class NMFWI : public BASENMF {
	public:
		NMFWI(int, int, int, bool, double);
		~NMFWI();
		int fit(gsl_matrix*, gsl_matrix*, gsl_matrix*);

	private:
		gsl_matrix* U;
		gsl_matrix* V;
		gsl_matrix* R;
		int update_kth_block_of_U(int);
		int update_kth_block_of_V(int);
		int update();
		int normalize_and_scale();
		double calculate_objective();
		int record_progress(double, double, double);
};
#endif
