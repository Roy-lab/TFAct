#include <gsl/gsl_matrix.h>
#ifndef _init_
#define _init_
#include <vector>
using namespace std;

namespace init{
	int nndsvd(gsl_matrix*, gsl_matrix*, gsl_matrix*, int);
	int nndsvd_general(gsl_matrix*, gsl_matrix*, gsl_matrix*, int);
	int nndsvd2(gsl_matrix*, gsl_matrix*, gsl_matrix*, int);
	int nndsvd(gsl_matrix*, gsl_matrix*, gsl_matrix*, gsl_matrix*, gsl_matrix*, int);
	int random(gsl_matrix*, gsl_matrix*, int);
	int random(vector<gsl_matrix*>&, int);
};
#endif
