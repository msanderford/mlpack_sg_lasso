#include <argparse.hpp>
#include <armadillo>
#include "sg_lasso_leastr.hpp"

using namespace std;
using namespace arma;


map<string, string> processSlepOpts(string filename);
string writeModelToXMLStream(string);

int main(int argc, char *argv[]) {
  argparse::ArgumentParser program("sg_lasso_leastr");

  program.add_argument("-f", "--features")
    .required()
    .help("Specify the input features file.");

  program.add_argument("-n", "--groups")
    .required()
    .help("Specify the group indices file.");

  program.add_argument("-r", "--response")
    .required()
    .help("Specify the response file.");

  program.add_argument("-w", "--output")
    .required()
    .help("specify the output file.");

  program.add_argument("-s", "--slep")
    .default_value(std::string("-"))
    .help("Specify a file of key/value SLEP options.");

  program.add_argument("-z", "--lambda1")
    .default_value(0.1)
    .help("Specify individual feature sparsity.")
    .scan<'g', double>();

  program.add_argument("-y", "--lambda2")
    .default_value(0.1)
    .help("Specify group feature sparsity.")
    .scan<'g', double>();

  try {
    program.parse_args(argc, argv);
  }
  catch (const std::runtime_error& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    std::exit(1);
  }

  //auto input = program.get<int>("square");
  //std::cout << (input * input) << std::endl;

  double lambda[2] = {program.get<double>("lambda1"), program.get<double>("lambda2")};

  mat features;
  mat opts_ind;
  rowvec responses;

  SGLassoLeastR* sgl;

  //features.load(csv_name(program.get<std::string>("features"),csv_opts::semicolon));
  features.load(csv_name(program.get<std::string>("features"),csv_opts::trans));

  //responses.load(csv_name(program.get<std::string>("response"),csv_opts::semicolon));
  responses.load(csv_name(program.get<std::string>("response"),csv_opts::trans));





  if (responses.n_cols != features.n_cols)
  {
    //Log::Fatal << "The responses must have the same number of columns as the feature set." << endl;
    throw std::invalid_argument("\nThe responses must have the same number of columns as the feature set.\n");
  }

  opts_ind.load(csv_name(program.get<std::string>("groups"),csv_opts::trans));

  sgl = new SGLassoLeastR(features, responses, opts_ind, lambda, processSlepOpts(program.get<std::string>("slep")));

  //std::cout << sgl->writeModelToXMLStream();
  ofstream fileStream(program.get<std::string>("output") + ".xml");
  if (fileStream.is_open())
  {
    //fileStream << sgl->writeModelToXMLStream();
    sgl->writeModelToXMLStream(fileStream);
    fileStream.close();
  } else {
    std::cout << "Could not open output file for writing." << std::endl;
  }

  return 0;

}


map<string, string> processSlepOpts(string filename)
{
  map<string, string> slep_opts;
  std::cout << "Processing SLEP options file: " << filename << "..." << std::endl;

  string line;
  ifstream optsFile (filename);

  if (optsFile.is_open())
  {

    int splitpos;
    string opt_key;
    while (getline(optsFile, line))
    {
      splitpos = line.find("\t");
      if (splitpos != std::string::npos)
      {
        opt_key = line.substr(0, line.find("\t"));
        slep_opts[opt_key] = line.substr(line.find("\t")+1, std::string::npos);
      }
    }
    optsFile.close();
  }

  return slep_opts;
}




