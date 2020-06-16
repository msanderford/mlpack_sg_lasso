#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>

using namespace std;

class alnData
{
	public:
		alnData(string speciesFile, string alnFileList);
		void readTraits(string fileName);
		void printTraits();
		void processFastaFileList(string fileName);
		void readAln(string fileName);
		void processAln();
		void generateResponseFile(string fileName);
		void generateFeatureFile(string fileName);
		void generateMappingFile(string fileName);
		void generateGroupIndicesFile(string baseName);
		void normalizeFeatures(bool normalize);
		void dropSingletons(bool ignoreSingletons);

	private:
		int featureIndex;
		int geneIndex;
		bool normalize;
		bool ignoreSingletons;
		string currentGene;
		vector<string> species;
		vector<string> groups;
		vector<vector<int>> groupIndices;
		map<string, int> geneGroupIndex;
		map<string, float> traits;
		map<string, string> seqs;
		map<int, string> featureMap;
		map<int, vector<int>> features;

};