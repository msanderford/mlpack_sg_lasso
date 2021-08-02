/**
 * @file sg_lasso_main.cpp
 * @author James Cline
 *
 * Main function for least-squares linear regression.
 *
 * mlpack is free software; you may redistribute it and/or modify it under the
 * terms of the 3-clause BSD license.  You should have received a copy of the
 * 3-clause BSD license along with mlpack.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#include <mlpack/prereqs.hpp>
#include <mlpack/core/util/cli.hpp>
#include <mlpack/core/util/mlpack_main.hpp>

#include "sg_lasso.hpp"

using namespace mlpack;
using namespace mlpack::regression;
using namespace mlpack::util;
using namespace arma;
using namespace std;

PROGRAM_INFO("Simple Linear Regression and Prediction",
    // Short description.
    "An implementation of sparse group lasso using "
    "logistic loss.  Given a dataset and responses, a model can be "
    "trained and saved for later use.",
    "An implementation of sparse group lasso using "
    "logistic loss.  Given a dataset and responses, a model can be "
    "trained and saved for later use.",
//        "output_model", "lr_model") +
//    "\n\n"
//    "Then, to use " + PRINT_MODEL("lr_model") + " to predict responses for a "
//    "test set " + PRINT_DATASET("X_test") + ", saving the predictions to " +
//    PRINT_DATASET("X_test_responses") + ", the following command could be "
//    "used:"
//    "\n\n" +
//    PRINT_CALL("sg_lasso", "input_model", "lr_model", "test", "X_test",
//        "output_predictions", "X_test_responses"),
    SEE_ALSO("Linear/ridge regression tutorial", "@doxygen/lrtutorial.html"),
    SEE_ALSO("@lars", "#lars"),
    SEE_ALSO("Linear regression on Wikipedia",
        "https://en.wikipedia.org/wiki/Linear_regression"),
    SEE_ALSO("mlpack::regression::SGLasso C++ class documentation",
        "@doxygen/classmlpack_1_1regression_1_1SGLasso.html"));

PARAM_MATRIX_IN("features", "Matrix containing feature set A.", "f");
PARAM_MATRIX_IN("opts_ind", "Matrix of indices defining overlapping group information.", "nodes");
PARAM_ROW_IN("responses", "Vector containing responses y ", "r");
PARAM_DOUBLE_IN("lambda1", "Feature regularization parameter (z1 >=0).", "z", 0.0);
PARAM_DOUBLE_IN("lambda2", "Group regularization parameter (z2 >=0).", "y", 0.0);
PARAM_STRING_IN("fastaFile", "FASTA alignment to be converted to one-hot encoded features file.", "seqs", "none");

// This is the future name of the parameter.
PARAM_MODEL_OUT(SGLasso, "feature_weights", "File of learned feature weights.", "weights");

map<string, string> processFasta(string filename);

static void mlpackMain()
{
  const double lambda1 = CLI::GetParam<double>("lambda1");
  const double lambda2 = CLI::GetParam<double>("lambda2");

  double lambda[2] = {lambda1, lambda2};

  mat features;
  mat opts_ind;
  rowvec responses;

  SGLasso* sgl;

  const bool computeModel = true;


  // An input file was given and we need to generate the model.
  if (computeModel)
  {
	string fastaFileName = CLI::GetParam<string>("fastaFile");
	//int testVal = processFasta(fastaFileName);

	//std::cout << "processFasta exit status: " << testVal << std::endl;

    Timer::Start("load_features");
    features = std::move(CLI::GetParam<mat>("features"));
    Timer::Stop("load_features");


    // The initial predictors for y, Nx1.
    Timer::Start("load_responses");
    responses = CLI::GetParam<rowvec>("responses");
    Timer::Stop("load_responses");

    if (responses.n_cols != features.n_cols)
    {
      Log::Fatal << "The responses must have the same number of columns "
          "as the feature set." << endl;
    }

    Timer::Start("load_opts_ind");
    opts_ind = std::move(CLI::GetParam<mat>("opts_ind"));
    Timer::Stop("load_opts_ind");

    Timer::Start("sparse_group_lasso");
    sgl = new SGLasso(features, responses, opts_ind, lambda, processFasta(fastaFileName));
    Timer::Stop("sparse_group_lasso");
  }


  // Save the model if needed.
  //CLI::GetParam<SGLasso*>("feature_weights") = sgl;
  CLI::GetParam<SGLasso*>("feature_weights") = sgl;
}


map<string, string> processFasta(string filename)
{
  map<string, string> slep_opts;
  std::cout << "Processing FASTA file: " << filename << "..." << std::endl;

  string line;
  ifstream fastaFile (filename);

  if (fastaFile.is_open())
  {
  /*
     int i = 0;
     while (getline(fastaFile, line))
     {
        i++;
	 }
	 fastaFile.close();
	 fastaFile.open(filename);
	 i = i/2;
	 string speciesIds[i];
	 string seqs[i];
     i = 0;
     while (getline(fastaFile, line))
     {
        if (line[0] == '>')
        {
           speciesIds[i] = line.substr(1,line.length()-1);
		} else {
           seqs[i] = line.substr(0,line.length());
           i = i + 1;
        }
     }
  */
    int splitpos;
    string opt_key;
    while (getline(fastaFile, line))
    {
      splitpos = line.find("\t");
      if (splitpos != std::string::npos)
      {
        opt_key = line.substr(0, line.find("\t"));
        slep_opts[opt_key] = line.substr(line.find("\t"), std::string::npos);
      }
    }	
  fastaFile.close();


//     for (int j = 0; j < i; j = j + 1)
//     {
//        std::cout << speciesIds[j] << ":" << seqs[j] << std::endl;
//     }

  }

  return slep_opts;
}
