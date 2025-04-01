#include "Framework.H"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <unistd.h>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <getopt.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <iomanip>

#include <pthread.h>
#include <semaphore.h>

#include "ThreadManager.H"
#include "MemoryCheck.H"

Framework::Framework()
{
	regMngr = NULL;
	tgtMngr = NULL;
	evidMngr = NULL;
	priorNet = NULL;
	lambda = 0;
	lambdap = 0;
	lambda_Cminus = 0;
	lambda_E = 0.01; // control cell lineage contribution
	CV=0;
	randomizedBCD=false;
}

Framework::~Framework()
{
	if (regMngr!=NULL)
	{
		delete regMngr;
	}
	if (tgtMngr!=NULL)
	{
		delete tgtMngr;
	}
	if (evidMngr!=NULL)
	{
		delete evidMngr;
	}
	if (priorNet!=NULL)
	{
		delete priorNet;
	}
	{
		map<string, pair<Graph*,double>>::iterator it;
		for(it = priorNets.begin(); it != priorNets.end(); ++it){
			delete (it->second).first;
		}
		priorNets.clear();
	}
	
}

int
Framework::inferNCA(EdgeList* initNet)
{
	//In case of duplicate (both expression and TFA regulator),
	//I am taking Max. Maybe I should take mean?
	EdgeList* initNet2 = initNet->removeSuffix("_nca");
	char ioutdir[1025];
	sprintf(ioutdir,"%s/",outdir);
	Graph* updatedPrior = new Graph;
	updatedPrior->setVariableManagers(regMngr,tgtMngr);
	updatedPrior->setNet(initNet2);
	//NCALearner* ncalrnr = new NCALearner(regMngr, tgtMngr, evidMngr, priorNet, lambda, 10, numThread, ioutdir);
	NCALearner* ncalrnr = new NCALearner(regMngr, tgtMngr, evidMngr, updatedPrior, lambda, CV, numNcaRep, ioutdir);
	ncalrnr->start();
	VariableManager* nregMngr;
	EvidenceManager* ntfaMngr;
	ncalrnr->getData(nregMngr,ntfaMngr,"_nca");
	delete ncalrnr;

	regMngr->mergeVarSets(nregMngr);

	delete ntfaMngr;
	delete initNet2;

	return 0;
}

// map<int, gsl_matrix*>* inputmat, map<int, gsl_matrix*>* A_in, map<int, gsl_matrix*>* P_in, 
// 				   map<int, pair<gsl_matrix*, double>>* p_net, double lam_p, 
// 				   map<int, map<int, double>>* cell_neighborhood, double lam_E, double lam_Cminus

int
Framework::inferNMF_ct(){
	char ioutdir[1025];
	sprintf(ioutdir,"%s/",outdir);
	int inferStatus;

	int nRegs = regMngr->getVarCnt();
	int nTgts = tgtMngr->getVarCnt();

	map<int, gsl_matrix*> inputmat;
	map<int, pair<gsl_matrix*, double>> p_net;
	{
		map<int, string>::iterator it;
		map<string, pair<Graph*, double>>::iterator priorit;
		for(it = IDToCell.begin(); it != IDToCell.end(); ++it){
			Matrix* dataMat = evidMngr->getDataMat(it->second); // it allocates a new Matrix. 
			int colCnt = dataMat->getColCnt();
			gsl_matrix* GsldataMat = gsl_matrix_calloc(nTgts, colCnt);
			dataMat->getGslMatrix(GsldataMat);
			inputmat[it->first] = GsldataMat;

			priorit = priorNets.find(it->second);
			if(priorit != priorNets.end()){
				gsl_matrix* priorMat = gsl_matrix_calloc(nTgts, nRegs);
				Matrix* p_matrix = ((priorNets[it->second]).first)->getAdj(); // it allocates a new Matrix. 
				p_matrix->getGslMatrix(priorMat);
				p_net[it->first] = make_pair(priorMat, (priorNets[it->second]).second);
				delete p_matrix;
			}
			delete dataMat;
		}
	}
	delete evidMngr;
	evidMngr = NULL;
	
	verbosity=true; //change this
	NMFLearner* nmfLearner = NULL;
	
	nmfLearner = new NMFLearner(nRegs, initMethod, maxIter, seed, verbosity, 1, nmft);

	// // pair<map<int, gsl_matrix**>*, map<int, vector<double>*>*>* nmfOut;
	nmfLearner->setPrior(&p_net, lambda);

	(nmfLearner->prior_info_public).lambda_p = lambdap;

	inferStatus = nmfLearner->inferNMF_ct(numNcaRep, &inputmat, lambdap, &speciesTree, lambda_E, lambda_Cminus);
	if(inferStatus != 0){
		printf("\n");
		return -1;
	}
	// nmfOut = nmfLearner->get_nmf_out();
	nmfLearner->writeToDir(ioutdir, &IDToCell);
	// gsl_matrix_free(GsldataMat);
	{
		map<int, gsl_matrix*>::iterator it;
		for(it=inputmat.begin(); it != inputmat.end(); ++it){
			gsl_matrix_free(it->second);
		}
		inputmat.clear();
		map<int, pair<gsl_matrix*, double>>::iterator p_net_it;
		for(p_net_it=p_net.begin(); p_net_it != p_net.end(); ++p_net_it){
			gsl_matrix_free((p_net_it->second).first);
		}
		p_net.clear();
	}
	// gsl_matrix_free(priorMat);
	// if (l_tgt_graph > 0){
	// 	gsl_matrix_free(UGraphLaplacian);
	// }
	// if (l_reg_graph > 0){
	// 	gsl_matrix_free(VGraphLaplacian);
	// }
	delete nmfLearner;
	return 0;
}

int
Framework::inferNMF()
{
	if((nmft == ibpg_2fct)||(nmft == KL_NMF_2fct))
	{
		int ret = inferNMF_ct();
		return ret;
	}
	char ioutdir[1025];
	sprintf(ioutdir,"%s/",outdir);
	int inferStatus;

	int nRegs = regMngr->getVarCnt();
	int nTgts = tgtMngr->getVarCnt();

	Matrix* dataMat = evidMngr->getDataMat();
	int colCnt = dataMat->getColCnt();
	gsl_matrix* GsldataMat = gsl_matrix_calloc(nTgts, colCnt);
	dataMat->getGslMatrix(GsldataMat);
	delete evidMngr;
	evidMngr = NULL;
	
	gsl_matrix* priorMat = gsl_matrix_calloc(nTgts, nRegs);
	Matrix* p_matrix = priorNet->getAdj();
	p_matrix->getGslMatrix(priorMat);
	delete priorNet;
	priorNet = NULL;
	
	// seed = static_cast<int>(time(0));
	seed = 1; // change this
	verbosity=true; //change this
	NMFLearner* nmfLearner = NULL;
	if((nmft == recur) || (nmft == nmqf)|| (nmft == ibpg)){
		nmfLearner = new NMFLearner(nRegs, k1, k2, initMethod, maxIter, 50, 50, seed, verbosity, 1, randomizedBCD, nmft);
	} else { // 2-factor
		nmfLearner = new NMFLearner(nRegs, initMethod, maxIter, seed, verbosity, 1, nmft);
	}
	// pair<map<int, gsl_matrix**>*, map<int, vector<double>*>*>* nmfOut;
	nmfLearner->setPrior(priorMat, lambda);

	(nmfLearner->prior_info_public).lambda_p = lambdap;


	gsl_matrix* UGraphLaplacian;
	gsl_matrix* VGraphLaplacian;
	SIMMAT sim;
	double UGraphLnorm = 0;
	double VGraphLnorm = 0;
	// TODO probably different similarity matrix
	
	if (l_tgt_graph > 0){
		UGraphLaplacian = gsl_matrix_alloc(priorMat->size1, priorMat->size1);
		sim.getGLFromPrior(priorMat, CblasLower, UGraphLaplacian, 0);
		UGraphLnorm = sim.simMatSpectralNorm(UGraphLaplacian);
		gsl_matrix_scale(UGraphLaplacian, l_tgt_graph); // scaling by lambda as required by IBPG
	} 
	if (l_reg_graph > 0){
		VGraphLaplacian = gsl_matrix_alloc(priorMat->size2, priorMat->size2);
		sim.getGLFromPrior(priorMat, CblasLower, VGraphLaplacian, 1);
		VGraphLnorm = sim.simMatSpectralNorm(VGraphLaplacian);
		gsl_matrix_scale(VGraphLaplacian, l_reg_graph);
	}
	nmfLearner->setLambdaUxUGraphL(UGraphLaplacian, l_tgt_graph, UGraphLnorm); // note that UGraphLaplacian is scaled by l_tgt_graph
	nmfLearner->setLambdaVxVGraphL(VGraphLaplacian, l_reg_graph, VGraphLnorm);

	inferStatus = nmfLearner->inferNMF(numNcaRep, GsldataMat);
	if(inferStatus != 0){
		printf("\n");
		return -1;
	}
	// nmfOut = nmfLearner->get_nmf_out();
	nmfLearner->writeToDir(ioutdir);
	gsl_matrix_free(GsldataMat);
	gsl_matrix_free(priorMat);
	if (l_tgt_graph > 0){
		gsl_matrix_free(UGraphLaplacian);
	}
	if (l_reg_graph > 0){
		gsl_matrix_free(VGraphLaplacian);
	}
	delete nmfLearner;
	return 0;
}

int
Framework::start()
{
	MemoryCheck mc;
	//Start with MotifTFA
	mc.begin();
	cout << "Starting NCA:" << endl;
	// inferNCA(pnet);
	inferNMF();
	mc.end();
	mc.print("inferNCA");

	return 0;
}

int
Framework::init(int argc, char** argv)
{
	bool inputG=false;
	bool inputR=false;
	bool inputE=false;
	bool inputN=false;
	bool inputL=false;
	bool inputC=false;
	bool inputO=false;
	bool inputT=false;
	bool inputk1=false;
	bool inputk2=false;
	bool inputNCA=false;
	bool inputIter=false;
	bool inputInit=false;
	bool inputRegGraphR=false;
	bool inputTgtGraphR=false;
	bool inputTree=false;
	bool seedset=false;

	char gname[1025];
	char rname[1025];
	char pname[1025];
	char dname[1025];
	char tname[1025];

	int optret='-';
	opterr=1;
	int oldoptind=optind;
	int condCnt=1;
	int long_option_index;

	static struct option long_options[] =
    {
		{"k_1", required_argument, 0, 1000},
		{"k_2", required_argument, 0, 1001},
		{"num_rep", required_argument, 0, 1002},
		{"max_iter", required_argument, 0, 1003},
		{"randomized", no_argument, 0, 1004},
		{"init_m", required_argument, 0, 1005},
		{"l_reg_graph", required_argument, 0, 1006},
		{"l_tgt_graph", required_argument, 0, 1007},
		{"lambda_p", required_argument, 0, 1008},
		{"tree", required_argument, 0, 1009},
		{"lambda_Cminus", required_argument, 0, 1010},
		{"lambda_E", required_argument, 0, 1011},
		{"seed", required_argument, 0, 1012},
		{0, 0, 0, 0}
    }; 

	while((optret=getopt_long(argc,argv,"r:g:p:o:d:l:c:t:h::", long_options, &long_option_index))!=-1)
	{
		if(optret=='?')
		{
			cout << "Option error " << optopt << endl;
			return -1;
		} 
		int c = optret;
		char* my_optarg=NULL;
		// c=*(argv[oldoptind]+1);
		if(optind-oldoptind == 2)        //for -v 1
		{
			my_optarg=argv[oldoptind+1];	
		}
		else                             //for -v1
		{
			my_optarg=argv[oldoptind]+2;
		}
		switch(c)
		{
			case 'h': //Help
			{
				printHelp(argv[0]);
				return -1;
			}
			case 'o': //Output directory
			{
				strcpy(outdir,my_optarg);
				inputO=true;
				break;
			}
			case 'l'://Lambda
			{
				lambda = atof(my_optarg);
				inputL=true;
				break;
			}
			case 'c'://CV
			{
				CV = atoi(my_optarg);
				inputC=true;
				break;
			}
			case 'd'://Expression
			{
				strcpy(dname,my_optarg);
				inputE=true;
				break;
			}
			case 'p'://Prior net
			{
				strcpy(pname,my_optarg);
				inputN=true;
				break;
			}
			case 'g'://Target genes
			{
				strcpy(gname,my_optarg);
				inputG=true;
				break;
			}
			case 'r'://Regulators
			{
				strcpy(rname,my_optarg);
				inputR=true;
				break;
			}
			case 't'://nmf types
			{	
				if (strcmp(my_optarg, "recur") == 0) {
                    nmft = recur;
                } else if (strcmp(my_optarg, "nmqf") == 0) {
                    nmft = nmqf;
                } else if (strcmp(my_optarg, "ibpg") == 0) {
                    nmft = ibpg;
				} else if (strcmp(my_optarg, "hals_prior") == 0) {
					nmft = hals_prior;
                }else if (strcmp(my_optarg, "ibpg_2f") == 0) {
					nmft = ibpg_2f;
                }else if (strcmp(my_optarg, "ibpg_2fct") == 0) { // cell type specific
					nmft = ibpg_2fct;
                }else if (strcmp(my_optarg, "KL_NMF_2fct") == 0) { // cell type specific
					nmft = KL_NMF_2fct;
                }else {
                    cerr << "Please supply a valid NMF method name." << endl;
                    return -1;
                }
				inputT=true;
				break;
			}
			case 1000://k1
			{
				k1 = atoi(optarg);
				if(k1 < 1){
					cerr << "Please provide a valid k_1" << endl;
					return -1;
				}
				inputk1=true;
				break;
			}
			case 1001://k2
			{
				k2 = atoi(optarg);
				if(k2 < 1){
					cerr << "Please provide a valid k_2" << endl;
					return -1;
				}
				inputk2=true;
				break;
			}
			case 1002://num_rep
			{
				numNcaRep = atoi(optarg);
				if(numNcaRep < 1){
					cerr << "Please provide a valid num_rep" << endl;
					return -1;
				}
				inputNCA=true;
				break;
			}
			case 1003://max_iter
			{
				maxIter = atoi(optarg);
				if(maxIter < 1){
					cerr << "Please provide a valid max iteration" << endl;
					return -1;
				}
				inputIter=true;
				break;
			}
			case 1004://randomized
			{
				randomizedBCD = true;
				break;
			}
			case 1005://init_m
			{
				if (strcmp(my_optarg, "nndsvd") == 0) {
                    initMethod = nndsvd_init;
                } else if (strcmp(my_optarg, "random") == 0) {
                    initMethod = random_init;
                }else {
                    cerr << "Please supply a valid initialization method name." << endl;
                    return -1;
                }
				inputInit = true;
				break;
			}
			case 1006://lambda on reg (matrix U graph regularizer)
			{
				l_reg_graph=atof(my_optarg);
				if (l_reg_graph < 0) {
                    cerr << "Please supply a nonnegative graph regularization lambda for regulators." << endl;
                    return -1;
                } else if (l_reg_graph > 0){
					printf("Graph regularize on regulators.\n");
				}
				inputRegGraphR = true;
				break;
			}
			case 1007://lambda on tgt (matrix V graph regularizer)
			{
				l_tgt_graph=atof(my_optarg);
				if (l_tgt_graph < 0) {
                    cerr << "Please supply a nonnegative graph regularization lambda for genes." << endl;
                    return -1;
                } else if (l_tgt_graph > 0){
					printf("Graph regularize on genes.\n");
				}
				inputTgtGraphR = true;
				break;
			}
			case 1008:
			{
				lambdap=atof(my_optarg);
				if (lambdap < 0) {
                    cerr << "Please supply a nonnegative lambda for regularizing TFA." << endl;
                    return -1;
                } else if (lambdap > 0){
					printf("l1 regularizes TFA.\n");
				}
				break;
			}
			case 1009: // species/cell type tree. 
			{
				strcpy(tname,my_optarg);
				inputTree=true;
				break;
			}
			case 1010:
			{
				lambda_Cminus=atof(my_optarg);
				if (lambda_Cminus < 0) {
                    cerr << "Please supply a nonnegative lambda for regularizing cell clusters that don't have prior information." << endl;
                    return -1;
                } else if (lambda_Cminus > 0){
					printf("l1 regularizes cell clusters that don't have prior information.\n");
				}
				break;
			}
			case 1011: // control contribution of the cell lineage (distances between cell types)
			{
				lambda_E=atof(my_optarg);
				if (lambda_E < 0) {
                    cerr << "Please supply a nonnegative lambda for controling the influence of the cell lineage." << endl;
                    return -1;
                } else if (lambda_E > 0){
					printf("l2 regularizes differences between cell clusters that are neighbors.\n");
				}
				break;
			}
			case 1012: // control contribution of the cell lineage (distances between cell types)
			{
				seed=atoi(my_optarg);
				seedset=true;
				break;
			}
			default:
			{
				cerr << "Unhandled option " << c  << endl;
				printHelp(argv[0]);
				return -1;
			}
		}
		oldoptind=optind;
	}
	if (!inputE)
	{
		cerr << "Expression matrix was not provided!" << endl;
		printHelp(argv[0]);
		return -1;
	}
	if (!inputN)
	{
		cerr << "Prior network was not provided!" << endl;
		printHelp(argv[0]);
		return -1;
	}
	if (!inputG)
	{
		cerr << "List of target genes was not provided!" << endl;
		printHelp(argv[0]);
		return -1;
	}
	if (!inputR)
	{
		cerr << "List of regulators was not provided!" << endl;
		printHelp(argv[0]);
		return -1;
	}
	if (!inputO)
	{
		cerr << "Output directory was not provided!" << endl;
		printHelp(argv[0]);
		return -1;
	}
	if (!inputT)
	{
		cerr << "nmf type not provided. Default to nmqf." << endl;
		nmft = nmqf;
	}
	if ((!inputTree) && ((nmft == ibpg_2fct) || (nmft == KL_NMF_2fct))){
		cerr << "Please provide a tab separated cell lineage tree." << endl;
		return -1;
	}
	if (!inputL && !inputC)
	{
		cerr << "Lambda for soft thresholding on matrix A was not provided" << endl;
		cerr << "Setting lambda to 0" << endl;
		lambda = 0;
	}
	if (!inputRegGraphR)
	{
		l_reg_graph = 0;
	}
	if (!inputTgtGraphR)
	{
		l_tgt_graph = 0;
	}
	if (inputL && inputC)
	{
		cerr << "Both Lambda and CV were provided" << endl;
		cerr << "Ignoring CV" << endl;
		CV=0;
	}
	if ((!inputk1) && (nmft == ibpg))
	{
		cerr << "k_1 not provided. Set to 8." << endl;
		k1 = 8;
	}
	if ((!inputk2) && (nmft == ibpg))
	{
		cerr << "k_2 not provided. Set to 8." << endl;
		k2 = 8;
	}
	if (!inputNCA)
	{
		numNcaRep = 1;
	}
	if (!inputIter)
	{
		cerr << "number of iterations not provided. Set to 300." << endl;
		maxIter = 300;
	}
	if (!inputInit){
		initMethod = nndsvd_init;
	}
	if (!seedset){
		cerr << "Seed not provided. Set to 1." << endl;
		seed = 1;
	}

	regMngr = new VariableManager;
	regMngr->readVariables(rname);

	tgtMngr = new VariableManager;
	tgtMngr->readVariables(gname);

	evidMngr = new EvidenceManager;
	evidMngr->setVariableManager(tgtMngr);
	if((nmft == ibpg_2fct) || (nmft == KL_NMF_2fct)){
		evidMngr->loadEvidenceFromMultipleFiles_Continuous(dname);
		{
			vector<string> allcell = evidMngr->getSpeciesNames();
			for(int i = 0; i < allcell.size(); ++i){
				cellToID[allcell[i]] = i;
				IDToCell[i] = allcell[i];
			}
		}
		ifstream inFile(pname);
		char buffer[2049];
		double cf;
		while(inFile.good()){
			inFile.getline(buffer,2048);
			if(strlen(buffer)<=0)
			{
				continue;
			}
			string cellName;
			char* priorName;
			char* tok=strtok(buffer,"\t");
			int tokCnt=0;
			while(tok!=NULL)
			{
				if(tokCnt==0)
				{
					cellName.append(tok);
				}
				else if(tokCnt==1)
				{
					priorName = tok;

				} 
				else if(tokCnt==2)
				{
					cf = atof(tok);
					Graph* prior_net = new Graph;
					prior_net->setVariableManagers(regMngr,tgtMngr);
					prior_net->readNet(priorName);
					pair<Graph*, double> p_to_insert(prior_net, cf);
					priorNets[cellName] = p_to_insert;
					break;
				}
				tok=strtok(NULL,"\t");
				tokCnt++;
			}
		}
		// load (undirected) tree
		inFile.close();  
		inFile.clear();
		inFile.open(tname);
		while(inFile.good()){
			inFile.getline(buffer,2048);
			if(strlen(buffer)<=0)
			{
				continue;
			}
			string cell1Name;
			string cell2Name;
			double dist;
			char* tok=strtok(buffer,"\t");
			int tokCnt=0;
			while(tok!=NULL)
			{
				if(tokCnt==0)
				{
					cell1Name.append(tok);
				}
				else if (tokCnt==1)
				{
					cell2Name.append(tok);
				}
				else if (tokCnt==2)
				{
					dist = atof(tok);
					speciesTree[cellToID[cell1Name]][cellToID[cell2Name]] = dist;
					speciesTree[cellToID[cell2Name]][cellToID[cell1Name]] = dist;
					break;
				}
				tok=strtok(NULL,"\t");
				tokCnt++;
			}
		}


	} else {
		evidMngr->loadEvidenceFromFile_Continuous(dname);
		priorNet = new Graph;
		priorNet->setVariableManagers(regMngr,tgtMngr);
		priorNet->readNet(pname);

		pnet = priorNet->getNetMap();
		pnet->addSuffix("_nca");
	}
	return 0;
}

int
Framework::printHelp(char* name)
{
	cerr << "Usage for " << name << endl;
	cerr << left << setw(20) << "-h" << "print this message" << endl;
	cerr << left << setw(20) << "-d" << "Input expression, tab delimited, [GeneName]\\t[val]\\t[val]... Pass a 2-column file of [celltype]\\t[path] if doing cell-type specific inference (TFAct)." << endl;
	cerr << left << setw(20) << "-r" << "List of regulators." << endl;
	cerr << left << setw(20) << "-g" << "List of target genes." << endl;
	cerr << left << setw(20) << "-p" << "Prior network, [TF]\\t[TG]\\t[Confidence]. Pass a 3-column file of [celltype]\\t[path]\\t[prior quality (0~1)] if doing cell-type specific inference." << endl;
	cerr << left << setw(20) << "-l" << "Penalty for l1 regularization on A. Higher means imposing a heavier penalty on non prior edges (default 0)" << endl;
	cerr << left << setw(20) << "-c" << "Number of cross validation folds (instead of fixed Lambda). Not implemented." << endl;
	cerr << left << setw(20) << "-o" << "Output directory. This will overwrite existing files." << endl;
	cerr << left << setw(20) << "-t" << "type of nmf. Support recur, nmqf, ibpg, hals_prior, ibpg_2f, ibpg_2fct (TFAct-F), KL_NMF_2fct (TFAct-KL)" << endl;
	cerr << left << setw(20) << "--k_1" << "(applicable to recur, nmqf, ibpg) dimension of rows of U in UV^TWH^T (default 8)" << endl;
	cerr << left << setw(20) << "--k_2" << "(applicable to recur, nmqf, ibpg) dimension of rows of W in UV^TWH^T (default 8)" << endl;
	cerr << left << setw(20) << "--num_rep" << "the # of NMF with different initializations (default 1)" << endl;
	cerr << left << setw(20) << "--max_iter" << "the # max iteration cycles to update all entries (default 300)" << endl;
	cerr << left << setw(20) << "--randomized" << "randomly choosing blocks to update (flag. default false. Applicable to nmqf, ibpg only)" << endl;
	cerr << left << setw(20) << "--init_m" << "initialization method. (default nndsvd, and should use random for TFAct)" << endl;
	cerr << left << setw(20) << "--l_reg_graph" << "(applicable to recur, nmqf, ibpg) graph regularization on matrix V (tf/regulators). (default 0)" << endl;
	cerr << left << setw(20) << "--l_tgt_graph" << "(applicable to recur, nmqf, ibpg) graph regularization on matrix U (genes). (default 0)" << endl;
	cerr << left << setw(20) << "--lambda_p" << "l1 regularization on the TFA matrix. (default 0)" << endl;
	cerr << left << setw(20) << "--tree" << "3-column, tab delim cell type lineage tree (Note that edges are distances, not similarity). Format is child\\tparent\\tdistance" << endl;
	cerr << left << setw(20) << "--seed" << "Seed. An integer." << endl;
	return 0;
}

/*
int
Framework::makeRandomSubs(int maxCnt)
{
	Matrix* d = allEvidMngr->getDataMat();
	int maxDim = d->getColCnt();
	for (int r=0;r<maxCnt;r++)
	{
		vector<int> allids;
		for (int i=0;i<maxDim;i++)
		{
			allids.push_back(i);
		}
		random_shuffle(allids.begin(),allids.end());
		vector<int> remids;
		for (int i=0;i<maxDim/2;i++)
		{
			remids.push_back(allids[i]);
		}
		allRMIDs.push_back(remids);
	}
	return 0;
}
*/
