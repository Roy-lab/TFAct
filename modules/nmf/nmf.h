#ifndef _nmf_
#define _nmf_
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <list>

#include "basenmf.h"
using namespace std;

class NMF : public BASENMF {
	public:
		NMF(int, init_method, int, int, bool, double);
		NMF(int, init_method, int, int, bool, double, bool);
		~NMF();
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
		int record_progress(double, double, double);
		gsl_rng* ri;
};
#endif
