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
		NMQF(int, init_method, int, int, bool, double, bool);
		~NMQF();
		int fit(gsl_matrix*, gsl_matrix*, gsl_matrix*);
		bool update_permuted;

	private:
		gsl_matrix* U;
		gsl_matrix* V;
		gsl_matrix* R;
		int update_kth_block_of_U(int);
		int update_kth_block_of_V(int);
		int update();
		int normalize_and_scale();
		double calculate_objective();
		gsl_rng* ri;
};
#endif
