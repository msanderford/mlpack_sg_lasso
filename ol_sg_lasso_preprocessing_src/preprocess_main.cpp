#include "preprocess.hpp"
#include <string.h>

using namespace std;

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

	//Open species traits file and read everything with -1/1 into table

	//int testVal = processFasta(argv[2]);

	alnData* data = new alnData(argv[1], argv[2]);

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
		}
	}
	
	cout << "Generating group indices file..." <<endl;
	data->generateGroupIndicesFile(basename);
	cout << "Generating mapping file..." <<endl;
	data->generateMappingFile(basename);
	cout << "Generating response file..." <<endl;
	data->generateResponseFile(basename);
	cout << "Generating feature file..." <<endl;
	data->generateFeatureFile(basename);


	return 0;
}
