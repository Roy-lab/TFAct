#include <iostream>
#include <fstream>
#include <queue>
#include <string>
#include <string.h>
#include "Vertex.H"
#include "Graph.H"
#include "SpeciesDistance.H"
#include "Framework.H"


Framework::Framework()
{
}

Framework::~Framework()
{
}

int 
Framework::readAncientNetwork(const char* aFName)
{
    //read edges in aFName and make graph:
	ancientNetwork.makeGraph(aFName);
	return 0;
}

int
Framework::readTFList(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);		
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string tfname;
		tfname.append(buffer);
		tflist[tfname]=0;   //map<string,int> tflist;
	}
	inFile.close();
	cout << "Read " << tflist.size() << " TFs" << endl;
	return 0;
}

int 
Framework::readSpeciesDistance(const char* aFName, double proot)
{
    spDist.setProot(proot);
	spDist.readSpeciesTree(aFName);
	return 0;
}

int 
Framework::setMoves(int aCnt)
{
	moveCnt=aCnt;
	return 0;
}

/*For each regulator-target state in the ancestral network, propagate its
* state to the leave nodes.*/

int 
Framework::generateNetworks()
{
	r=gsl_rng_alloc(gsl_rng_default);
	map<string,Vertex*>& vertexSet=ancientNetwork.getVertexSet();
	for(map<string,int>::iterator vIter=tflist.begin();vIter!=tflist.end();vIter++)
	{
		Vertex* tf=NULL;
		map<string,Vertex*>* neighbors=NULL;
		if(vertexSet.find(vIter->first)!=vertexSet.end())
		{
			tf=vertexSet[vIter->first];
			neighbors=&(tf->getImmediateNeighbours());
		}
		//For now assume that the tf must be connected to soembody in the original network
		if(tf==NULL)
		{
			continue;
		}
		for(map<string,Vertex*>::iterator nIter=vertexSet.begin();nIter!=vertexSet.end();nIter++)
		{
			if(strcmp(vIter->first.c_str(),nIter->first.c_str())==0)
			{
				continue;
			}
			// if tf --> nIter->first in original network edgeStatus=1
			int edgeStatus=0;
			if(neighbors->find(nIter->first)==neighbors->end())
			{
				edgeStatus=0;
			}
			else
			{
				edgeStatus=1;
			}
			map<string,int> edgeStatusInTree;
                        //cout << "Edge status of: " << vIter->first.c_str() << "->" << nIter->first.c_str() << endl;
			spDist.genEdgeStatus(r,edgeStatus,edgeStatusInTree);
			for(map<string,int>::iterator eIter=edgeStatusInTree.begin();eIter!=edgeStatusInTree.end();eIter++)
			{
				if(eIter->second==0)
				{
					continue;
				}
				Graph* currentNw=NULL;
				if(extantNetworks.find(eIter->first)==extantNetworks.end())
				{
					currentNw=new Graph;
					extantNetworks[eIter->first]=currentNw;
				}
				else
				{
					currentNw=extantNetworks[eIter->first];
				}
				currentNw->addEdge(vIter->first.c_str(),nIter->first.c_str());
				//cout << eIter->first << " add edge: " << vIter->first.c_str() << "->" << nIter->first.c_str() << endl;
 			}
		}
	}
	map<SpeciesDistance::EdgeStatusType,int>& edgeTypes=spDist.getEdgeStatus();
	for(map<SpeciesDistance::EdgeStatusType,int>::iterator eIter=edgeTypes.begin();eIter!=edgeTypes.end();eIter++)
	{
		switch(eIter->first)
		{
			case SpeciesDistance::ON_ON:
			{
				cout <<"Edge ON_ON " << eIter->second << endl;
				break;
			}
			case SpeciesDistance::ON_OFF:
			{
				cout <<"Edge ON_OFF " << eIter->second << endl;
				break;
			}
			case SpeciesDistance::OFF_ON:	
			{
				cout <<"Edge OFF_ON " << eIter->second << endl;
				break;
			}
			case SpeciesDistance::OFF_OFF:
			{
				cout <<"Edge OFF_OFF " << eIter->second << endl;
				break;
			}
		}
	}
	return 0;
}

int
Framework::dumpNetworks(const char* outDir)
{
    cout << "dumpNetworks: " << extantNetworks.size() << endl;
	for(map<string,Graph*>::iterator gIter=extantNetworks.begin();gIter!=extantNetworks.end();gIter++)
	{
		char aDirName[1024];
		sprintf(aDirName,"%s/%s",outDir,gIter->first.c_str());
		char commandName[1024];
		sprintf(commandName,"mkdir -p %s",aDirName);
		system(commandName);
		dumpNetwork(gIter->second,aDirName);
	}	
	return 0;
}


int 
Framework::dumpNetwork(Graph* graph, const char* outputDir)
{
	map<string,Vertex*>& vertexSet=graph->getVertexSet();
	char aFName[1024];
	sprintf(aFName,"%s/nw.tsv",outputDir);
	ofstream nwFile(aFName);
	cout << "dumpNetworks: " << aFName;
	sprintf(aFName,"%s/nw.sif",outputDir);
	ofstream regFile(aFName);
	sprintf(aFName,"%s/genes.txt",outputDir);
	ofstream gFile(aFName);
	sprintf(aFName,"%s/proteins.txt",outputDir);
	ofstream pFile(aFName);
	int nedge=0;
	for(map<string,Vertex*>::iterator vIter=vertexSet.begin();vIter!=vertexSet.end();vIter++)
	{
		Vertex* v=vIter->second;
		gFile << vIter->first << endl;
		pFile << vIter->first << endl;
		map<string,Vertex*>& targets=v->getImmediateNeighbours();
		for(map<string,Vertex*>::iterator nIter=targets.begin();nIter!=targets.end();nIter++)
		{
			double pval=gsl_ran_flat(r,0,1);
			if(pval<=0.5)
			{
				regFile << vIter->first << "\tactivates\t" << nIter->first << endl;
			}
			else
			{
				regFile << vIter->first <<"\trepresesses\t" << nIter->first << endl;
			}
			nwFile << vIter->first<<"\t" << nIter->first << endl;
			nedge++;
		}
	}
	nwFile.close();
	gFile.close();
	pFile.close();	
	cout <<" Number of edges: " << nedge << endl;
	return 0;
}

//update by Shilu
int 
Framework::estimatePhyloTopoDist()
{
	SpeciesDistance::Species* root=spDist.getRoot();
	queue<SpeciesDistance::Species*> specQueue;
	specQueue.push(root);
	while(!specQueue.empty())
	{
		SpeciesDistance::Species* node=specQueue.front();
		specQueue.pop();
        for(int i=0;i<node->children.size();i++)
        {
            SpeciesDistance::Species* child=node->children[i];
            speciesGraph.addEdge(node->name.c_str(),child->name.c_str());
            speciesGraph.addEdge(child->name.c_str(),node->name.c_str());
            if(child->children.size()>0)
            {
                specQueue.push(child);
            }
        }
	}
	speciesGraph.obtainConnectivity();
	for(map<string,Graph*>::iterator aIter=extantNetworks.begin();aIter!=extantNetworks.end();aIter++)
	{
		Graph* nw1=aIter->second;
		map<string,Graph*>::iterator bIter=aIter;
		bIter++;
		Vertex* u=speciesGraph.getVertex(aIter->first.c_str());
		for(;bIter!=extantNetworks.end();bIter++)
		{
			Graph* nw2=bIter->second;
			Vertex* v=speciesGraph.getVertex(bIter->first.c_str());
			double pdist=u->getPathLength(v);
			string key(aIter->first.c_str());
			key.append("-");
			key.append(bIter->first.c_str());
			phyloDist[key]=pdist;
			double edist=getNwSim(nw1,nw2);
			nwDist[key]=edist;
		}
	}
	for(map<string,double>::iterator pIter=phyloDist.begin();pIter!=phyloDist.end();pIter++)
	{
		cout << pIter->first <<"\t" << pIter->second <<"\t" << nwDist[pIter->first] << endl;
	}
	return 0;
}

double
Framework::getNwSim(Graph* nw1, Graph* nw2)
{
	map<string,Vertex*>& vset1=nw1->getVertexSet();
	map<string,Vertex*>& vset2=nw2->getVertexSet();
	double edges1=0;
	double common=0;
	double edges2=0;
	for(map<string,Vertex*>::iterator vIter=vset1.begin();vIter!=vset1.end();vIter++)
	{
		Vertex* v1=vIter->second;
		map<string,Vertex*>& neighbors1=v1->getImmediateNeighbours();
		if(neighbors1.size()==0)
		{
			continue;
		}	
		edges1=edges1+neighbors1.size();
		if(vset2.find(vIter->first)==vset2.end())
		{
			continue;
		}
		Vertex* v2=vset2[vIter->first];
		map<string,Vertex*>& neighbors2=v2->getImmediateNeighbours();
		for(map<string,Vertex*>::iterator nIter=neighbors1.begin();nIter!=neighbors1.end();nIter++)
		{
			if(neighbors2.find(nIter->first)!=neighbors2.end())
			{
				common++;
			}
		}
	}

	for(map<string,Vertex*>::iterator vIter=vset2.begin();vIter!=vset2.end();vIter++)
	{
		Vertex* v=vIter->second;
		edges2=edges2+v->getImmediateNeighbours().size();
	}
	double prec=common/edges2;
	double recall=common/edges1;
	double fscore=(2*prec*recall)/(prec+recall);
	return fscore;
}

int
main(int argc, const char** argv)
{
	if(argc!=6)
	{
		cout << "simulateNw inputnetwork species_dist tflist outputdir proot" << endl;
		return 0;
	}
	Framework fw;
    double proot = atof(argv[5]);
    cout << "P(root=1)= " << proot << endl;
	fw.readAncientNetwork(argv[1]);
	fw.readSpeciesDistance(argv[2],proot);
	fw.readTFList(argv[3]);
	fw.generateNetworks();
	fw.dumpNetworks(argv[4]);
	fw.estimatePhyloTopoDist();
	return 0;
}
