#include <math.h>
#include <iostream>
#include <map>

#include "Distance.H"


Distance::Distance()
{
}

Distance::~Distance()
{
}

double
Distance::computeSymmKLDivergence(double m1, double v1, double m2, double v2)
{
	double kl1=computeKLDivergence(m1,v1,m2,v2);
	double kl2=computeKLDivergence(m2,v2,m1,v1);
	double klsim=(kl1+kl2)/2;
	return klsim;
}

double 
Distance::computeKLDivergence(double m1, double v1, double m2, double v2)
{
	double t1=log(v2/v1);
	double t2=((m1-m2)*(m1-m2))/(v2);
	double t3=v1/v2;
	double kld=(t1+t2+t3-1)/2;
	return kld;
}

double 
Distance::computeZstat(double m1,double v1,double m2, double v2,int sampleCnt)
{
	double temp1=m1-m2;
	double temp2=sqrt((v1+v2)/sampleCnt);
	double zstat=temp1/temp2;
	return zstat;
}
/*
double
Distance::computeCC(vector<double>& v1, vector<double>& v2)
{
	double cc=0;
	double m1=0;
	for(int i=0;i<v1.size();i++)
	{
		m1=m1+v1[i];
	}
	m1=m1/v1.size();

	double m2=0;
	for(int i=0;i<v2.size();i++)
	{
		m2=m2+v2[i];
	}
	m2=m2/v2.size();
	
	double xx=0;
	double yy=0;
	double xy=0;
	double oppRel=0;
	for(int i=0;i<v1.size();i++)
	{
		double diff1=v1[i]-m1;
		xx=xx + (diff1 * diff1);
		double diff2=v2[i]-m2;
		yy=yy+ (diff2 * diff2);
		xy=xy+(diff1*diff2);  	
		if( (diff1*diff2) < 0)
		{
			oppRel++;
		}
	}
	cc=sqrt((xy*xy)/(xx*yy));
	double threshold=v1.size()/2.0;

	if(oppRel > threshold)
	{
		cc=cc*(-1);
	}
	return cc;
}*/


double
Distance::computeCC(vector<double>& v1, vector<double>& v2)
{
	double cc=0;
	double m1=0;
	for(int i=0;i<v1.size();i++)
	{
		m1=m1+v1[i];
	}

	double m2=0;
	for(int i=0;i<v2.size();i++)
	{
		m2=m2+v2[i];
	}
	
	double sumx=m1;
	double sumy=m2;
	double xx=0;
	double yy=0;
	double xy=0;
	for(int i=0;i<v1.size();i++)
	{
		double xi=v1[i];
		xx=xx + (xi* xi);
		double yi=v2[i];
		yy=yy+ (yi*yi);
		xy=xy+(xi*yi);  	
	}
	double n=(double) v1.size();
	double num=n*xy-(sumx*sumy);
	double den1=sqrt(n*xx-(sumx*sumx));
	double den2=sqrt(n*yy-(sumy*sumy));
	cc=num/(den1*den2);
	return cc;
}


double
Distance::computeSpearmanCC(vector<double>& v1, vector<double>& v2)
{
	double cc=0;
	//First get the ranks
	vector<double> rank1(v1.size());
	vector<double> rank2(v2.size());
	getRanks(v1,rank1);
	getRanks(v2,rank2);
	cc=computeCC(rank1,rank2);
	return cc;
}

int
Distance::getRanks(vector<double>& values, vector<double>& ranks)
{
	vector<int> inds;
	for(int i=0;i<values.size();i++)
	{
		inds.push_back(i);
	}
	for(int i=0;i<values.size();i++)
	{
		for(int j=i+1;j<values.size();j++)
		{
			if(values[inds[i]]>values[inds[j]])
			{
				int oldind=inds[i];
				inds[i]=inds[j];
				inds[j]=oldind;
			}
		}
	}

	for(int i=0;i<inds.size();i++)
	{
		ranks[inds[i]]=i;	
	}
	//Resolve ties
	int i=0;
	while(i<(ranks.size()-1))
	{
		double ranksum=ranks[inds[i]];
		int j=i+1;
		while( (j<ranks.size()) && (values[inds[j]]==values[inds[i]]) )
		{
			ranksum=ranksum+ranks[inds[j]];
			j++;
		}
		double avgrank=ranksum/(j-i);
		if((j-i)>1)
		{
	//		cout <<"Found tie at " << i <<  " to " << j << " avgrank " << avgrank << endl;
		}
		for(int k=i;k<j;k++)
		{
			ranks[inds[k]]=avgrank;
		}
		i=j;
	}
	
	return 0;
}
