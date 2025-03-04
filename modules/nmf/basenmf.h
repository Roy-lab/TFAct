#ifndef _basenmf_
#define _basenmf_
#include <gsl/gsl_matrix.h>
#include <list>
using namespace std;

enum init_method {nndsvd_init, random_init, graph_init};

class BASENMF {
	protected:
		gsl_matrix* X; //To be factorized
	public:
		virtual double calculate_objective() = 0;
		virtual int record_progress(double, double, double) = 0;
		int n_components;
		init_method init;
		int max_iter;
		int random_state;
		bool verbose;
		double tol;
		list<double> reconstruction_err_;
		list<double> relative_err_;
		list<double> cpu_time_;
		virtual ~BASENMF() {}
};
#endif
