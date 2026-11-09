#include <fstream>
#include <iostream>
#include <string.h>
#include <math.h>
#include "SpeciesDistance.H"

SpeciesDistance::SpeciesDistance()
{
    root=NULL;
    proot=0.5;
}

SpeciesDistance::~SpeciesDistance()
{
}

int
SpeciesDistance::readSpeciesTree(const char* aFName)
{
    ifstream inFile(aFName);
    cout << "SpeciesDistance::readSpeciesTree: " << aFName << endl;
    char buffer[1024];
    while(inFile.good())
    {
        inFile.getline(buffer,1023);
        if(strlen(buffer)<=0)
        {
            continue;
        }
        if(strchr(buffer,'#')!=NULL)
        {
            continue;
        }
        char* tok=strtok(buffer,"\t");
        int tokCnt=0;
        string childSpeciesName;
        string parentSpeciesName;
        //string childtype;
        double p_gain=0;
        double p_maintain_edge=0;
        while(tok!=NULL)
        {
            if(tokCnt==0)
            {
                childSpeciesName.append(tok);
            }
            else if(tokCnt==1)
            {
                parentSpeciesName.append(tok);
            }
            else if(tokCnt==2)
            {
                p_gain=atof(tok);
            }
            else if(tokCnt==3)
            {
                p_maintain_edge=atof(tok);
            }
            tok=strtok(NULL,"\t");
            tokCnt++;
        }
        SpeciesDistance::Species* childSpecies=NULL;
        SpeciesDistance::Species* parentSpecies=NULL;
        if(speciesSet.find(childSpeciesName)==speciesSet.end())
        {
            childSpecies=new SpeciesDistance::Species;
            childSpecies->name.append(childSpeciesName.c_str());
            childSpecies->parent=NULL;
            speciesSet[childSpeciesName]=childSpecies;
        }
        else
        {
            childSpecies=speciesSet[childSpeciesName];
        }
        if(speciesSet.find(parentSpeciesName)==speciesSet.end())
        {
            parentSpecies=new SpeciesDistance::Species;
            parentSpecies->name.append(parentSpeciesName.c_str());
            parentSpecies->parent=NULL;
            speciesSet[parentSpeciesName]=parentSpecies;
        }
        else
        {
            parentSpecies=speciesSet[parentSpeciesName];
        }
        childSpecies->parent=parentSpecies;
        if(parentSpecies->parent==NULL)
        {
            root=parentSpecies;
        }
        //Probabilities are for the child node
        childSpecies->p_maintain_edge=p_maintain_edge;
        childSpecies->p_gain=p_gain;
        parentSpecies->children.push_back(childSpecies);
        cout << childSpeciesName << " p_maintain_edge=" << p_maintain_edge << " p_gain=" << p_gain << endl;
    }
    cout <<"speciesSet.size()=" <<speciesSet.size() << " Root is " << root->name << endl;
    inFile.close();
    return 0;
}


//The probability that an edge is maintained in a child given that the edge is present in the ancestor
double
SpeciesDistance::getProbGainGain(string& spName)
{
    if(speciesSet.find(spName)==speciesSet.end())
    {
        cout <<"No species with name " << spName.c_str() << endl;
        exit(0);
    }
    Species* species=speciesSet[spName];
    double gain_gain=species->p_maintain_edge;
    return gain_gain;
}

//The probability that an edge is gained in a child given that the edge is absent in the ancestor
double
SpeciesDistance::getProbGainLoss(string& spName)
{
    if(speciesSet.find(spName)==speciesSet.end())
    {
        cout <<"No species with name " << spName.c_str() << endl;
        exit(0);
    }
    Species* species=speciesSet[spName];
    double gain_loss=species->p_gain;
    return gain_loss;
}

//The probability than an edge is lost in a child given that the edge was present in the ancestor
double
SpeciesDistance::getProbLossGain(string& spName)
{
    if(speciesSet.find(spName)==speciesSet.end())
    {
        cout <<"No species with name " << spName.c_str() << endl;
        exit(0);
    }
    Species* species=speciesSet[spName];
    double loss_gain=1-species->p_maintain_edge;
    return loss_gain;
}

//The  probability that an edge is lost in a child given the edge was absent in the ancestor
double
SpeciesDistance::getProbLossLoss(string& spName)
{
    if(speciesSet.find(spName)==speciesSet.end())
    {
        cout <<"No species with name " << spName.c_str() << endl;
        exit(0);
    }
    Species* species=speciesSet[spName];
    double loss_loss=1-species->p_gain;
    return loss_loss;
}

// new getEdgeStatusProb by Shilu
double
SpeciesDistance::getEdgeStatusProb(map<string,int>& edgeStatus)
{
    double edgeprior=0;
    int binkey=0;
    int binpos=0;
    for(map<string,int>::iterator eIter=edgeStatus.begin();eIter!=edgeStatus.end();eIter++)
    {
        if(eIter->second!=0)
        {
            int val=(int)pow(2,binpos);
            binkey=binkey+val;
        }
        binpos=binpos+1;
    }
    if(binEdgeKeyProbMap.find(binkey)!=binEdgeKeyProbMap.end())
    {
        return binEdgeKeyProbMap[binkey];
    }
    //Start with the root
    //changed 6/4 to test setting the prior probability of an edge in the root node = .2
    double score;
    if(edgeStatus[root->name]==0)
    {
        score = 1-proot;
    }
    else
    {
        score = proot;
    }
    edgeprior=score;
    cout << "P(G)=P(" << root->name <<"=" <<edgeStatus[root->name] <<")[" << score <<"]";
    for(int i=0;i<root->children.size();i++)
    {
        double childrenScore=getSubTree(edgeStatus[root->name],root->children[i],edgeStatus);
        edgeprior=edgeprior*childrenScore;
        cout << " * P(" << root->children[i]->name <<"=" <<edgeStatus[root->children[i]->name] <<"|" << root->name <<"=" <<edgeStatus[root->name] << ")[" <<  childrenScore << "]";
    }
    cout << "=" << edgeprior << endl;
    binEdgeKeyProbMap[binkey]=edgeprior;
    return edgeprior;
}


// new getSubTree function by Shilu
double
SpeciesDistance::getSubTree(bool parentEdge, Species* child, map<string,int>& edgeStatus)
{
    double score=0;
    int status=edgeStatus[child->name];
    if(child->children.empty())
    {
        //This is a leaf node
        if(parentEdge)
        {
            if(status==0)
            {
                score=getProbLossGain(child->name);
            }
            else
            {
                score=getProbGainGain(child->name);
            }
        }
        else
        {
            if(status==0)
            {
                score=getProbLossLoss(child->name);
            }
            else
            {
                score=getProbGainLoss(child->name);
            }
        }
    }
    else
    {
        //This is a parent node

        // its own score:
        if(edgeStatus[child->parent->name])
        {
            if(status==0)
            {
                score=getProbLossGain(child->name);
            }
            else
            {
                score=getProbGainGain(child->name);
            }
        }
        else
        {
            if(status==0)
            {
                score=getProbLossLoss(child->name);
            }
            else
            {
                score=getProbGainLoss(child->name);
            }
        }
        //cout << "Current " << child->name << " score=" << score << endl;
        //cout << " * P(" << child->name <<"=" <<edgeStatus[child->name] <<"|" << child->parent->name <<"=" <<edgeStatus[child->parent->name] << ")" ;
        // children score
        for(int i=0;i<child->children.size();i++)
        {
            // edgeStatus[child->name] status for current node, i.e. parent edge status for children
            double childrenScore=getSubTree(edgeStatus[child->name],child->children[i],edgeStatus);
            score=score*childrenScore;
            //cout << " * P(" << child->children[i]->name <<"=" <<edgeStatus[child->children[i]->name] <<"|" << child->name <<"=" <<edgeStatus[child->name] << ")" ;
        }
    }
    return score;
}

SpeciesDistance::Species*
SpeciesDistance::getRoot()
{
    return root;
}

// new genEdgeStatus function by Shilu
int
SpeciesDistance::genEdgeStatus(gsl_rng* r,int status, map<string,int>& edgeStatus)
{
    rndmizer=r;
    //Shilu: set edge status for root cell controlled by proot
    double pval=gsl_ran_flat(rndmizer,0,1);
    int rootstatus=-1;
    //similar to pval<=pselect, if p(root=1)>=cutoff(pval) rootstatus=1;
    if(pval<=proot)
    {
        rootstatus=1;
    }
    else
    {
        rootstatus=0;
    }
    edgeStatus[root->name]=rootstatus;
    //cout << root->name << " pval=" << pval << " proot=" << proot << " rootstatus=" << rootstatus << endl;
    for(int i=0;i<root->children.size();i++)
    {
        Species* child=root->children[i];
        sampleEdgeStatus_Subtree(child,rootstatus,rootstatus,edgeStatus);
    }
    return 0;
}

// new sampleEdgeStatus_Subtree function by Shilu
//If the node is a leaf parent then add to edgeStatus otherwise just sample
int
SpeciesDistance::sampleEdgeStatus_Subtree(Species* node, int rootstatus, int parentstatus, map<string,int>& edgeStatus)
{
    double pval=gsl_ran_flat(rndmizer,0,1);
    double pselect=-1;
    if(parentstatus==0)
    {
        pselect=getProbGainLoss(node->name);
    }
    else
    {
        pselect=getProbGainGain(node->name);
    }
    int childStatus=-1;
    if(pval<=pselect)
    {
        childStatus=1;
    }
    else
    {
        childStatus=0;
    }
    SpeciesDistance::EdgeStatusType etype=SpeciesDistance::ON_ON;
    if(rootstatus==0 && childStatus==0)
    {
        etype=OFF_OFF;
    }
    else if(rootstatus==0 && childStatus==1)
    {
        etype=OFF_ON;
    }
    else if(rootstatus==1 && childStatus==0)
    {
        etype=ON_OFF;
    }
    else
    {
        etype=ON_ON;
    }
    if(edgeTypeAnalysis.find(etype)==edgeTypeAnalysis.end())
    {
        edgeTypeAnalysis[etype]=1;
    }
    else
    {
        edgeTypeAnalysis[etype]=edgeTypeAnalysis[etype]+1;
    }
    //cout << node->name << " pval=" << pval << " pselect=[GL(parent=0) or GG(parent=1)]=" << pselect << " rootstatus=" << rootstatus <<" parentstatus=" << parentstatus << " childStatus=" << childStatus << endl;
    // leaf node
    if(node->children.empty())
    {
        edgeStatus[node->name]=childStatus;
    }
    else
    {
        //get the edge status for parent node
        edgeStatus[node->name]=childStatus;
        for(int i=0;i<node->children.size();i++)
        {
            sampleEdgeStatus_Subtree(node->children[i], rootstatus, childStatus, edgeStatus);
        }
    }
    return 0;
}

map<SpeciesDistance::EdgeStatusType,int>&
SpeciesDistance::getEdgeStatus()
{
    return edgeTypeAnalysis;
}

int
SpeciesDistance::setProot(double p)
{
    proot = p;
    return 0;
}

