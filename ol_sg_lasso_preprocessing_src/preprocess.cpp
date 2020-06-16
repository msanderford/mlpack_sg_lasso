#include "preprocess.hpp"
#include <algorithm>
#include <sstream>

using namespace std;

void trim(string& s)
{
   size_t p = s.find_first_not_of(" \t\r\n");
   s.erase(0, p);

   p = s.find_last_not_of(" \t\r\n");
   if (string::npos != p)
      s.erase(p+1);
}


alnData::alnData(string speciesFile, string alnFileList)
{
	//Read species traits/response file
	this->readTraits(speciesFile);
	//this->printTraits();


	cout << "Processing FASTA files list: " << alnFileList << "..." << endl;


	//Read and process FASTA files
	this->featureIndex = 1;
	this->processFastaFileList(alnFileList);
	this->normalize = false;

}

void alnData::normalizeFeatures(bool normalize)
{
	this->normalize = normalize;
}

void alnData::dropSingletons(bool ignoreSingletons)
{
	this->ignoreSingletons = ignoreSingletons;
}

void alnData::readTraits(string speciesFile)
{
	string speciesTrait;
	ifstream speciesList (speciesFile);
	if (speciesList.is_open())
	{
		string delimiter = "	";
		string species;
		float trait;
		while (getline(speciesList, speciesTrait))
		{
			trim(speciesTrait);
			species = speciesTrait.substr(0, speciesTrait.find(delimiter));
			trait = stof(speciesTrait.substr(speciesTrait.find(delimiter), string::npos));
			if (trait != 0)
			{
				this->species.push_back(species);
				this->traits[species] = trait;
			}
		}
		speciesList.close();
	}
}

void alnData::printTraits()
{
	for (int i = 0; i < this->species.size(); i++)
	{
		cout << this->species[i] << "	" << this->traits[this->species[i]] << endl;
	}
}

void alnData::processFastaFileList(string alnFileList)
{
	string fastaFileGroup, fastaFileName;
	ifstream fileList (alnFileList);
	if (fileList.is_open())
	{
		while (getline(fileList,fastaFileGroup))
		{
			trim(fastaFileGroup);
			this->groups.push_back(fastaFileGroup);
			stringstream ss(fastaFileGroup);
			while(getline(ss, fastaFileName, ','))
			{
				if (this->geneGroupIndex.find(fastaFileName) == this->geneGroupIndex.end())
				{
					this->geneGroupIndex[fastaFileName] = this->geneGroupIndex.size();
					this->readAln(fastaFileName);
					fastaFileName.erase(fastaFileName.find(".fas"), string::npos);
					fastaFileName.erase(0, fastaFileName.find("/")+1);
					this->currentGene = fastaFileName;
					this->processAln();
					this->seqs.clear();
				}
			}
		}
		fileList.close();
	}
}

void alnData::readAln(string fastaFileName)
{
	string line;
	int seqlen;
	vector<string> tempSpecies;
	tempSpecies = this->species;
	ifstream fastaFile (fastaFileName);
	if (fastaFile.is_open())
	{
		cout << "Processing FASTA file: " << fastaFileName << "..." << endl;
		string seqid;
		while (getline(fastaFile, line))
		{
			trim(line);
			if (line[0] == '>')
			{
				seqid = line.substr(1,string::npos);
				getline(fastaFile,line);
				if (this->traits.find(seqid) != this->traits.end())
				{
					trim(line);
					this->seqs[seqid] = line;
					seqlen = line.length();
					tempSpecies.erase(find(tempSpecies.begin(), tempSpecies.end(), seqid));
				}
			}
		}
		fastaFile.close();
	}
	while (!tempSpecies.empty())
	{
		this->seqs[tempSpecies.at(0)] = string(seqlen, '-');
		tempSpecies.erase(tempSpecies.begin());
	}
}

void alnData::processAln()
{
	int seqLen = this->seqs[this->species[1]].size();
	int numSpecies = this->species.size();
	int groupStartIndex = this-> featureIndex;

	for (int i = 0; i < seqLen; i++)
	{
		//check if site is constant
		set<char> bases;
		bases.insert('-');
		for (int j = 0; j < numSpecies; j++)
		{
			bases.insert(this->seqs[this->species[j]][i]);
		}
		//process it into features if it isn't
		if (bases.size() > 2)
		{
			set<char>::iterator baseIter = bases.begin();
			for (int j = 0; j < bases.size(); j++)
			{
				if (*baseIter != '-')
				{
					vector<int> oneHot;
					string featureName = this->currentGene + "_" + to_string(i) + "_" + *baseIter;
					this->featureMap[this->featureIndex] = featureName;
					for (int k = 0; k < numSpecies; k++)
					{
						if (this->seqs[this->species[k]][i] == *baseIter)
						{
							oneHot.push_back(1);
						} else {
							oneHot.push_back(0);
						}
					}
					this->features[featureIndex] = oneHot;
					this->featureIndex++;
				}
				baseIter++;
			}
		}
	}
	this->groupIndices.push_back({groupStartIndex,this->featureIndex-1});
}

void alnData::generateResponseFile(string baseName)
{
	string responseFileName = "response_" + baseName + ".txt";
	ofstream responseFile (responseFileName);
	if (responseFile.is_open())
	{
		for (int i = 0; i < this->species.size(); i++)
		{
			responseFile << this->traits[this->species[i]] << endl;
		}
		responseFile.close();
	}
}

void alnData::generateFeatureFile(string baseName)
{
	string featuresFileName = "feature_" + baseName + ".txt";
	//transpose features first for efficiency
	vector<vector<float>> tFeatures;
	for (int i = 0; i < this->species.size(); i++)
	{
		tFeatures.push_back({});
	}
	for (int i = 0; i < this->features.size(); i++)
	{
		vector<int> oneHot = this->features[i];
		float sum = 0;
		float val;
		for (int j = 0; j < this->features[i].size(); j++)
		{
			sum += this->features[i][j];
		}
		val = 1.0;
		if (this->normalize)
		{
			val = 1.0/sum;
		}
		if (this->ignoreSingletons && sum == 1.0)
		{
			val = 0.0;
		}
		for (int j = 0; j < oneHot.size(); j++)
		{
			tFeatures[j].push_back(val * oneHot[j]);
		}
	}
	//write to file
	ofstream featuresFile (featuresFileName);
	if (featuresFile.is_open())
	{
		for (int i = 0; i < tFeatures.size(); i++)
		{
			string featureLine;
			for (int j = 0; j < tFeatures[i].size(); j++)
			{
				if (tFeatures[i][j] == 0.0)
				{
					featureLine.append(to_string(0) + '\t');
				}
				else
				{
					featureLine.append(to_string(tFeatures[i][j]) + '\t');
				}
			}
			featureLine.pop_back();
			featuresFile << featureLine << endl;
		}
		featuresFile.close();
	}

}

void alnData::generateMappingFile(string baseName)
{
	string mappingFileName = "feature_mapping_" + baseName + ".txt";
	ofstream mappingFile (mappingFileName);
	if (mappingFile.is_open())
	{
		for (int i = 0; i < this->featureMap.size(); i++)
		{
			mappingFile << i << "\t" << this->featureMap[i] << endl;
		}
		mappingFile.close();
	}
}

void alnData::generateGroupIndicesFile(string baseName)
{
	string groupIndicesFileName = "group_indices_" + baseName + ".txt";
	string fieldFileName = "field_" + baseName + ".txt";
	ofstream groupIndicesFile (groupIndicesFileName);
	ofstream fieldFile (fieldFileName);
	string indStarts, indEnds, weights, weightBuffer, gene;
	double weight;
	int geneStart, geneEnd;
	int groupStart = 1;
	int groupEnd = 0;
	if (groupIndicesFile.is_open() and fieldFile.is_open())
	{
		for (int i = 0; i < this->groups.size(); i++)
		{
			stringstream ss(this->groups[i]);
			while(getline(ss, gene, ','))
			{
				geneStart = groupIndices[geneGroupIndex[gene]][0];
				geneEnd = groupIndices[geneGroupIndex[gene]][1];
				for (int j = geneStart; j <= geneEnd; j++)
				{
					fieldFile << to_string(j) + "\t";
				}
				groupEnd += geneEnd-geneStart+1;
			}
			indStarts.append(to_string(groupStart) + "\t");
			indEnds.append(to_string(groupEnd) + "\t");
			weight = sqrt(1+(groupEnd - groupStart));
			weights.append(to_string(weight) + "\t");
			groupStart = groupEnd + 1;
		}
		indStarts.pop_back();
		indEnds.pop_back();
		weights.pop_back();

		groupIndicesFile << indStarts << endl;
		groupIndicesFile << indEnds << endl;
		groupIndicesFile << weights << endl;

		groupIndicesFile.close();
		fieldFile.close();
	}

}

