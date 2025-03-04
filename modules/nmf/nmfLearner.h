#ifndef _nmfLearner
#define _nmfLearner
#include <gsl/gsl_matrix.h>
#include <list>
#include <map>
#include <vector>
#include <iostream>
#include <stddef.h>
#include <stdio.h>

#include "nmf.h"
#include "nmqf.h"
#include "ibpg.h"
#include "nmf_w_prior.h"
#include "ibpg_2factor.h"
#include "ibpg_2fct.h"
#include "KL_NMF_2fct.h"

using namespace std;
enum nmfType {recur, nmqf, ibpg, hals_prior, ibpg_2f, ibpg_2fct, KL_NMF_2fct};

class NMFLearner{
    public:
        NMFLearner(int, int, int, init_method, int, int, int, int, bool, double, bool, nmfType);
		NMFLearner(int, init_method, int, int, bool, double, nmfType);
        ~NMFLearner();
        int setPrior(gsl_matrix*, double);
        int setPrior(map<int, pair<gsl_matrix*, double>>*, double);
        int setLambdaUxUGraphL(gsl_matrix*, double, double);
        int setLambdaVxVGraphL(gsl_matrix*, double, double);
        int inferNMF(int, gsl_matrix*);
        int inferNMF_ct(int, map<int, gsl_matrix*>*, double, map<int, map<int, double>>*, double, double);
        int writeToDir(char*);
        int writeToDir(char*, map<int, string>*);
        struct Prior_Info_Public {
            double lambda_p;
        };
        Prior_Info_Public prior_info_public;
        // pair<map<int, gsl_matrix**>*, map<int, vector<double>*>*>* 
        // get_nmf_out();

    private:
        static void* recur_helper(void*);
        static void* nmqf_helper(void*);
        static void* ibpg_helper(void*);
        static void* hals_prior_helper(void*);
        static void* ibpg2f_helper(void*);
        static void* ibpg2fct_helper(void*);
        static void* KL_NMF2fct_helper(void*);
        int print_matrix(FILE*, const gsl_matrix*);
        int inferRecur(gsl_matrix*);
        int inferNMQF(gsl_matrix*);
        int inferIBPG(gsl_matrix*);
        int inferIBPG2F(gsl_matrix*, double);
        int inferHals_prior(gsl_matrix*);
        int inferIBPG2FCT(map<int, gsl_matrix*>*, double, map<int, map<int, double>>*, double, double);
        int inferKL_NMF2fct(map<int, gsl_matrix*>*, double, map<int, map<int, double>>*, double, double);
        struct graph_reg {
            gsl_matrix* GraphL;
            gsl_matrix* lambdaxGraphL; // we only need lambda * GraphL for IBPG
            double lambda;
            double GraphLNorm;
        };
        struct prior_info {
            double lambda_p;
        };
        struct nmf_arg_struct {
            gsl_matrix* A;
            gsl_matrix* P;
            gsl_matrix* U;
            gsl_matrix* V;
            gsl_matrix* W;
            gsl_matrix* H;
            gsl_matrix* datMat;
            gsl_matrix* prior;
            double lambda1;
            double lambdap;
            BASENMF* bnmf;
            BASENMF* nmf_UV;
            BASENMF* nmf_WH;
            graph_reg UGreg;
            graph_reg VGreg;
            prior_info p_info;
        };
        gsl_matrix* priorMat;
        double lambda_1;
        // gsl_matrix* UGraphL;
        gsl_matrix* lambdaUxUGraphL;
        double lambda_U;
        double UGraphLNorm;
        // gsl_matrix* VGraphL;
        gsl_matrix* lambdaVxVGraphL;
        double lambda_V;
        double VGraphLNorm;
        nmfType nmft;
        gsl_matrix** allA;
        gsl_matrix** allP;
        gsl_matrix** allU;
        gsl_matrix** allV;
        gsl_matrix** allW;
        gsl_matrix** allH;

        
        BASENMF** nmfArray;
        init_method init;
        int numNcaRep;
        int n_components;
        int max_iter;
        int UV_maxIter;
        int WH_maxIter;
        int random_state;
        int k1; // inner dimension of UV
        int k2; // inner dimension of WH
        bool verbose;
        bool permuted_update;
        double tol;
        map<int, gsl_matrix**> ret_mat;

        //first map{0:A,1:P,2:Us,3:Vs,4:Ws,5:Hs}. Second map: {rep:reconstruction loss}
        // pair<map<int, gsl_matrix**>*, map<int, vector<double>*>*> factor_mat_w_loss; 

        /*
         * Here we start inferring cell type specific networks
        */
        map<int, map<int, gsl_matrix*>> allA_CT; // cell type specific
        map<int, map<int, gsl_matrix*>> allP_CT;
        map<int, pair<gsl_matrix*, double>>* p_net;

        struct nmf_ct_arg_struct {
            map<int, gsl_matrix*>* A;
            map<int, gsl_matrix*>* P;
            map<int, gsl_matrix*>* datMat;
            map<int, pair<gsl_matrix*, double>>* prior;
            map<int, map<int, double>>* cell_neighborhood;
            double lambda1;
            double lambdap;
            double lambda_E;
            double lambda_Cminus;
            BASENMF* bnmf;
        };

        

};
#endif