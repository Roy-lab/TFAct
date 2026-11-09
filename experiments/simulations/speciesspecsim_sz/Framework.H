/*Here is how the simulator will work: We will read an input network, which will
 * give us the probability of edge presence and absence. We will also have
 * the conditional distributions in the SpeciesDistance class. We will use the ancient
 * input network to compute the probability of an edge presence/absence in the first
 * ancestor. We will also accept as input the number of edge moves we will attempt.
 * After generating the networks, we will output for each species
 * a) network: pdint.txt, genes.txt and proteins.txt in  a format that can be read by renco*/

#ifndef _FRAMEWORK_
#define _FRAMEWORK_
#include <map>
#include <string>
#include "gsl/gsl_randist.h"
using namespace std;

class Graph;
class SpeciesDistance;
class Framework
{
	public:
	Framework();
	~Framework();
	int readAncientNetwork(const char*);
	int readSpeciesDistance(const char*, double);
	int setMoves(int);
	int generateNetworks();
	int dumpNetworks(const char*);
	int readTFList(const char*);
	int estimatePhyloTopoDist();
	private:
	int estimateEdgePresenceProb();
	int dumpNetwork(Graph*, const char*);
	double getNwSim(Graph*,Graph*);
	Graph ancientNetwork;
	map<string,Graph*> extantNetworks;
	double edgePresenceProb;
	SpeciesDistance spDist;
	int moveCnt;
	map<string,int> tflist;
	gsl_rng* r;
	map<string,double> phyloDist;
	map<string,double> nwDist;
	Graph speciesGraph;
	
};
#endif
