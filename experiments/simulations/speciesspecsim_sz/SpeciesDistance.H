#ifndef _SPECIES_DATA_
#define _SPECIES_DATA_

#include <map>
#include <string>
#include <vector>
#include "gsl/gsl_randist.h"
using namespace std;
class SpeciesDistance
{
	public:
	SpeciesDistance();
	~SpeciesDistance();
	int readSpeciesTree(const char*);
	double getEdgeStatusProb(map<string,int>& edgeStatus);
	double getProbGainGain(string&);
	double getProbGainLoss(string&);
	double getProbLossGain(string&);
	double getProbLossLoss(string&);
	int genEdgeStatus(gsl_rng *r,int status, map<string,int>& edgeStatus);
	// add by Shilu to set prior probability of an edge in the root p(root=1)
    int setProot(double p);
	struct Species;

	struct Species 
	{
		string name;
		Species* parent;
        //update by Shilu
        //Species* leftchild;
        //Species* rightchild;
        vector<Species*> children;
		double p_maintain_edge;
		double p_gain;	
	};
	
	Species* getRoot();
	typedef enum
	{
		ON_ON=0,
		ON_OFF,
		OFF_ON,
		OFF_OFF
		
	}EdgeStatusType;
	double getSubTree(bool,Species*,map<string,int>&);
	int sampleEdgeStatus_Subtree(Species* node, int rootstatus, int parentstatus, map<string,int>& edgeStatus);
	map<EdgeStatusType,int>& getEdgeStatus();
	private:
	map<string,Species*> speciesSet;
	Species* root;
	gsl_rng* rndmizer;
	map<SpeciesDistance::EdgeStatusType,int> edgeTypeAnalysis;
	//Shilu  store prior
    map<int,double> binEdgeKeyProbMap;
    double proot;
};
#endif
