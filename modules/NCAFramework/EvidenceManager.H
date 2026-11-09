#ifndef _EVIDENCE_MANAGER
#define _EVIDENCE_MANAGER

#include <map>
#include <vector>
#include <fstream>
#include <sys/types.h>
#include <unistd.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>

#include "Error.H"
#include "VariableManager.H"
#include "Matrix.H"
#include "CommonTypes.H"

using namespace std;

class EvidenceManager
{
	public:
		EvidenceManager();
		~EvidenceManager();

		//This associates the variable information with the variable id
		int setVariableManager(VariableManager*);
		void skipfirstline();
		
		//Read the evidence file
		Error::ErrorCode loadEvidenceFromFile_Continuous(const char*);

		Error::ErrorCode loadEvidenceFromMultipleFiles_Continuous(const char*);

		int setDataMat(Matrix* mat);
		Matrix* getDataMat();
		Matrix* getDataMat(const string&);
		vector<string> getSpeciesNames();

		int removeVarsByID(vector<int>& ids);

		int writeEvidence(const char* oname);
		int updateEvidence(EvidenceManager* eMngr);
	private:
		int readExpFromFile(const char* inname, map<string,vector<double>* >& expMap);

		int loadEvidenceFromMultipleFiles_Continuous_helper(const string&, const char*, map<string,vector<double>* >&);

		VariableManager* vMgr;

		Matrix* dataMat;

		map<string, Matrix*> dataMats;

		vector<string> allcell;

		bool skip_firstline;
};

#endif //_EVIDENCE_MANAGER
