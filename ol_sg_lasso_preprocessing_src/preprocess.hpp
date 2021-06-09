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
		alnData();
		void initialize(string speciesFile, string alnFileList);
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
		void setCountThreshold(int countThreshold);
		void setUpsampleBalance(bool upsampleBalance);
		void setDownsampleBalance(bool downsampleBalance);
		void setDiskCaching(bool useDiskCache);
		void balanceSample();

	private:
		int featureIndex;
		int geneIndex;
		bool normalize;
		bool ignoreSingletons;
		bool upsampleBalance;
		bool downsampleBalance;
		bool useDiskCache;
		int countThreshold;
		string currentGene;
		vector<string> species;
		vector<string> groups;
		vector<vector<int>> groupIndices;
		map<string, int> geneGroupIndex;
		map<string, float> traits;
		map<string, string> seqs;
		map<int, string> featureMap;
		map<int, vector<int>> features;
		map<int, string> featureCacheFiles;

};
