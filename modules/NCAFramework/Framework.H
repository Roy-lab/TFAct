#ifndef _FRAMEWORK
#define _FRAMEWORK

#include <vector>
#include <map>
#include <string>

#include "Matrix.H"
#include "NCALearner.H"
#include "VariableManager.H"
#include "EvidenceManager.H"
#include "Graph.H"

#include "../nmf/nmf.h"
#include "../nmf/similarity_mat.h"
#include "../nmf/nmfLearner.h"

using namespace std;

class Framework
{
public:
	Framework();
	~Framework();
	int init(int argc, char** argv);
	int start();
private:
	VariableManager* regMngr;
	VariableManager* tgtMngr;
	VariableManager* allGenesMngr;
	EvidenceManager* evidMngr;
	Graph* priorNet;
	map<string, pair<Graph*, double>> priorNets;
	EdgeList* pnet;
	double lambda; //LASSO penalty
	double l_reg_graph; //lambda for graph regularization on regulators (V)
	double l_tgt_graph; //lambda for graph regularization on target genes (U)
	double lambdap; //for regularizing p
	double lambda_Cminus; // for regularizing cells that don't have priors.
	double lambda_E;
	int CV;//CV
	int numNcaRep;
	int k1;
	int k2;
	int maxIter;
	int seed;
	bool verbosity;
	init_method initMethod;
	// similarity_type sim_type;
	bool randomizedBCD;
	char outdir[1025];

	int printHelp(char* name);

	int inferNCA(EdgeList* initNet);
	int inferNMF();
	int inferNMF_ct();
	nmfType nmft;
	map<string, int> cellToID;
	map<int, string> IDToCell;
	map<int, map<int, double>> speciesTree;
};

#endif //_FRAMEWORK
