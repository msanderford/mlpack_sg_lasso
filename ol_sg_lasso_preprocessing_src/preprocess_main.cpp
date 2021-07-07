#include "preprocess.hpp"
#include <string.h>
#include <sys/stat.h>

using namespace std;

bool pathExists(const std::string &s)
{
	struct stat buffer;
	return (stat (s.c_str(), &buffer) == 0);
}

int main(int argc, char *argv[])
{
	cout << "You have entered " << argc - 1 << " arguments:" << "\n";

	if (argc < 4)
	{
		cout << "Preprocessor requires at least 3 arguments." << endl;
		return 0;
	}

	for (int i = 1; i < argc; ++i) cout << argv[i] << "\n";

	string basename = argv[3];

	if (!pathExists(basename))
	{
		system(("mkdir -p " + basename).c_str());
	}

	//Open species traits file and read everything with -1/1 into table

	//int testVal = processFasta(argv[2]);

	alnData* data = new alnData();

	if (argc > 4)
	{
		for (int i = 3; i < argc; i++)
		{
			if (strcmp(argv[i], "n") == 0)
			{
				data->normalizeFeatures(true);
				cout << "Normalizing one-hot features..." << endl;
			}
			if (strcmp(argv[i], "is") == 0)
			{
				data->dropSingletons(true);
				cout << "Ignoring singleton mutations..." << endl;
			}
			if (strcmp(argv[i], "ct") == 0)
                        {
                                cout << "Ignoring mutations observed fewer than " << argv[i+1] << "times..." << endl;
				data->setCountThreshold(std::stoi(argv[i+1]));
                        }
			if (strcmp(argv[i], "ub") == 0)
			{
				data->setUpsampleBalance(true);
			}
			if (strcmp(argv[i], "db") == 0)
			{
				data->setDownsampleBalance(true);
			}
			if (strcmp(argv[i], "useCaching") == 0)
			{
				data->setDiskCaching(true);
			}
		}
	}
	
	data->initialize(argv[1], argv[2]);
	
	cout << "Generating group indices file..." <<endl;
	data->generateGroupIndicesFile(basename);
	cout << "Generating mapping file..." <<endl;
	data->generateMappingFile(basename);
	cout << "Generating response file..." <<endl;
	data->generateResponseFile(basename);
	cout << "Generating feature file..." <<endl;
	data->generateFeatureFile(basename);
	cout << "Generating missing sequences list..." <<endl;
	data->generateMissingFile(basename);


	return 0;
}
