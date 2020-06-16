/**
 * @file sg_lasso_leastr.cpp
 * @author James Cline
 * @author Michael Fox
 *
 * Implementation of simple linear regression.
 *
 * mlpack is free software; you may redistribute it and/or modify it under the
 * terms of the 3-clause BSD license.  You should have received a copy of the
 * 3-clause BSD license along with mlpack.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#include "overlapping_sg_lasso_leastr.hpp"
#include <mlpack/core/util/log.hpp>
#include "overlapping.hpp"
#include <algorithm>

using namespace mlpack;
using namespace mlpack::regression;

//OLSGLassoLeastR::OLSGLassoLeastR(const arma::mat& features,
//                                   const arma::rowvec& responses,
//                                   const double lambda,
//                                   const bool intercept) :
//    OLSGLassoLeastR(features, responses, arma::rowvec(), lambda, intercept)
//{}

OLSGLassoLeastR::OLSGLassoLeastR(const arma::mat& features,
                                   const arma::rowvec& responses,
                                   const arma::mat& weights,
                                   const arma::rowvec& field,
                                   double* lambda,
                                   const bool intercept) :
    lambda(lambda),
    intercept(intercept)
{
  Train(features, responses, weights, field, intercept);
}

//double OLSGLassoLeastR::Train(const arma::mat& features,
//                               const arma::rowvec& responses,
//                               const bool intercept)
//{
//  return Train(features, responses, arma::rowvec(), intercept);
//}

arma::rowvec& OLSGLassoLeastR::Train(const arma::mat& features,
                               const arma::rowvec& responses,
                               const arma::mat& weights,
                               const arma::rowvec& field,
                               const bool intercept)
{
  this->intercept = intercept;

  //Set all optional parameters to defaults
  int opts_maxIter = 5000;
  int opts_init = 2;
  int opts_tFlag = 3;
  int opts_nFlag = 0;
  int opts_rFlag = 1;
  int opts_mFlag = 0;
  int opts_rStartNum = 100;
  double opts_tol = 0.00001;
  arma::mat opts_ind = weights;
  opts_ind.cols(0,1) = opts_ind.cols(0,1) - 1;
  arma::rowvec opts_field = field - 1;

  //Set overlapping specific parameters to defaults

  int opts_maxIter2 = 1000;
  double opts_tol2 = 0.00000001;
  int opts_flag2 = 2;

  /*
   * We want to calculate the a_i coefficients of:
   * \sum_{i=0}^n (a_i * x_i^i)
   * In order to get the intercept value, we will add a row of ones.
   */

  // We store the number of rows and columns of the features.
  // Reminder: Armadillo stores the data transposed from how we think of it,
  //           that is, columns are actually rows (see: column major order).
  const size_t nCols = features.n_cols;

//  arma::mat p = features;
//  arma::rowvec r = responses;


  arma::mat A = features.t();
  arma::mat& ind = opts_ind;
  arma::colvec y = responses.t();
  double* z;
  z = this->Lambda();
  double lambda2_max;
  const size_t m = A.n_rows;
  const size_t n = A.n_cols;
  arma::colvec m_ones(m, arma::fill::ones);
  arma::colvec m_zeros(m, arma::fill::zeros);
  arma::colvec n_zeros(n, arma::fill::zeros);

  double lambda1 = z[0];
  double lambda2 = z[1];

  double *gap;
  gap = (double*) malloc(sizeof(double));
  double penalty2 [5];

  if (lambda1<0 || lambda2<0)
  {
	Log::Fatal << "\n z should be nonnegative!\n" << std::endl;
  }

  if(ind.n_cols != 3)
  {
	Log::Fatal << "\n Check opts_ind, expected 3 cols" << std::endl;
  }

  int groupNum = ind.n_rows;
  arma::colvec Y(opts_field.n_cols, arma::fill::zeros);

  //sgLogisticR.m:177-195
  //arma::colvec sample_weights(m);
  //sample_weights.fill(1.0/m);

//std::cout << "1..." << std::endl;

  //sgLogisticR.m:200-202
  //arma::uvec p_flag = arma::find(y == 1);
  //arma::uvec not_p_flag = arma::find(y != 1);
  //arma::colvec b(m);
//std::cout << "p_flag.n_elem:" << p_flag.n_elem << std::endl;
  //double m1 = static_cast<double>(p_flag.n_elem) / (double)m;
  //double m2 = 1 - m1;

  //sgLeastR.m:168-175

  arma::mat ATy = A.t() * y;


  //sgLeastR.m:178-201
  double* lambda;

  if (opts_rFlag == 0)
  {
	  lambda = z;
  } else {
	 if (lambda1<0 || lambda1>1 || lambda2<0 || lambda2>1)
	 {
		Log::Fatal << "\n opts.rFlag=1, and z should be in [0,1]" << std::endl;
	 }

	 arma::mat temp = arma::abs(ATy);

	 double lambda1_max = arma::as_scalar(arma::max(temp));

	 lambda1 = lambda1 * lambda1_max;

	 temp = arma::max(temp - lambda1, n_zeros);

	 //lambda2_max = computeLambda2Max(temp, n, ind, ind.n_rows);
	 lambda2_max = 1;

	 lambda2 = lambda2 * lambda2_max;
  }

  //sgLogisticR.m:203-216
  arma::colvec x(n, arma::fill::zeros);   //x.fill(0);

  if (opts_init != 2){
	  x = ATy;
  }
//std::cout << "4..." << std::endl;


//std::cout << "A_cols:" << A.n_cols << " A_rows:" << A.n_rows << " x_cols:" << x.n_cols << " x_rows:" << x.n_rows << std::endl;

  //sgLeastR.m:219-226
  arma::mat Ax = A * x;

//std::cout << "Ax_cols:" << Ax.n_cols << " Ax_rows:" << Ax.n_rows << std::endl;

//std::cout << "2..." << std::endl;

  //sgLeastR.m:228-253
  int bFlag = 0;
  double L = 1.0;

  //arma::colvec weighty = sample_weights % y;

  arma::colvec xp = x;   arma::colvec Axp = Ax;   arma::colvec xxp(n, arma::fill::zeros);
  //double cp = c;   double ccp = 0;

  double alphap = 0;   double alpha = 1;


  //sgLeastR.m:255-382
  double beta, sc, gc, fun_s, fun_x, l_sum, r_sum, tree_norm;
//  arma::mat As;
  arma::colvec aa, bb, prob, s, v;
  arma::rowvec ValueL(opts_maxIter);
  arma::rowvec funVal(opts_maxIter);
  //arma::mat ind_work(ind.n_rows + 1, ind.n_cols);

//std::cout << "3..." << std::endl;

std::cout << "m:" << m << " n:" << n << std::endl;

opts_ind = opts_ind.t();

  for (int iterStep = 1; iterStep <= opts_maxIter; iterStep = iterStep + 1)
  {
    //sgLogisticR.m:293-304
    beta = (alphap - 1)/alpha;   s = x + (xxp * beta);
//std::cout << "sc:" << sc << " c:" << c << " beta:" << beta << " ccp:" << ccp << " m1:" << m1 << " m2:" << m2 << std::endl;
//std::cout << "As_elems:" << As.n_elem << "Ax_elems:" << Ax.n_elem << " Axp_elems:" << Axp.n_elem << " beta:" << beta << std::endl;
    arma::mat As = Ax + ((Ax - Axp) * beta);

    //sgLogisticR.m:318-329
//std::cout << "A_rows:" << A.n_rows << " A_cols:" << A.n_cols << " b_elems:" << b.n_elem << std::endl;
    arma::mat ATAs = A.t() * As;
    arma::mat g = ATAs - ATy;
    xp = x;   Axp = Ax;

    int firstFlag = 1;
    //sgLogisticR.m:331-378
    while (true)
    {
      v = s - g/L;
      //sgLogisticR.m:337-338
//      arma::mat ind_work = ind;
//      ind_work.col(2) = ind_work.col(2) * (lambda2/L);
//      arma::rowvec first_row = {-1, -1, lambda1/L};
//      ind_work.insert_rows(0, first_row);
      //ind_work(0,0) = -1;
      //ind_work(0,1) = -1;
      //ind_work(0,2) = lambda1/L;
      //ind_work.submat(1,0,ind_work.n_rows,1) = ind.cols(0,1);
      //ind_work(arma::span(1,ind_work.n_rows),2) = ind.col(2) * (lambda2/L);

      //sgLogisticR.m:340-342
      //x = altra(v, n, ind_work, ind_work.n_rows);

//std::cout << "x:" << std::endl;
//std::cout << x << std::endl;
//std::cout << "v:" << std::endl;
//std::cout << v << std::endl;
//std::cout << "Y:" << std::endl;
//std::cout << Y << std::endl;
//std::cout << "opts_tol2:" << std::endl;
//std::cout << opts_tol2 << std::endl;

      overlapping(x.memptr(), gap, penalty2, v.memptr(),  n, groupNum, lambda1/L, lambda2/L, opts_ind.memptr(), opts_field.memptr(), Y.memptr(), opts_maxIter2, opts_flag2, opts_tol2);

//std::cout << "gap:" << std::endl;
//std::cout << gap << std::endl;
//std::cout << "penalty2:" << std::endl;
//std::cout << penalty2 << std::endl;
//std::cout << "opts_field:" << std::endl;
//std::cout << opts_field << std::endl;
//std::cout << "opts_ind:" << std::endl;
//std::cout << opts_ind << std::endl;
//std::cout << "x:" << std::endl;
//std::cout << x << std::endl;
//std::cout << "v:" << std::endl;
//std::cout << v << std::endl;

//      arma::colvec temp_x = altra(v, n, ind_work, ind_work.n_rows);
//std::cout << "temp_x_rows:" << temp_x.n_rows << " temp_x_cols:" << temp_x.n_cols << " x_rows:" << x.n_rows << " x_cols:" << x.n_cols << std::endl;
//std::cout << "9..." << std::endl;
      v = x - s;
      int nonzero_x_count = 0;
      for (int i = 0; i <= x.n_rows; i = i + 1)
      {
		 //std::cout << "aa["<< i << "]:" << aa(i) << " bb[" << i << "]:" << bb(i) << " Ax[" << i << "]:" << Ax(i) << std::endl;
		 if (x(i) != 0) { nonzero_x_count = nonzero_x_count + 1;}
	  }
//std::cout << "nonzero_x_count:" << nonzero_x_count << std::endl;

      //sgLogisticR.m:345-353
      Ax = A * x;
      //sgLogisticR.m:356-377
      arma::mat Av = Ax - As;
      r_sum = arma::as_scalar(v.t() * v); l_sum = arma::as_scalar(Av.t() * Av);
//std::cout << "r_sum:" << r_sum << " l_sum:" << l_sum << " L:" << L << std::endl;
//throw 20;

      if (r_sum <= std::pow(0.1, 20))
      {
	     bFlag = 1;
	     break;
	  }

	  if (l_sum <= r_sum * L)
	  {
	     break;
	  } else {
	     L = std::max(2*L, l_sum/r_sum);
	  }
    }
    //sgLogisticR.m:382-401
    alphap = alpha;   alpha = (1 + std::pow(4 * alpha * alpha + 1.0, 0.5))/2.0;

    xxp = x - xp; arma::mat Axy = Ax - y;

    ValueL(iterStep) = L;
    //funVal(iterStep) = fun_x;

    //ind_work(0,0) = -1;
    //ind_work(0,1) = -1;
    //ind_work(0,2) = lambda1;
    //ind_work.submat(1,0,ind_work.n_rows,1) = ind.cols(0,1);
    //ind_work(arma::span(1,ind_work.n_rows),2) = ind.col(2) * lambda2;
//    arma::mat ind_work = ind;
//    ind_work.col(2) = ind_work.col(2) * (lambda2);
//    arma::rowvec first_row = {-1, -1, lambda1};
//    ind_work.insert_rows(0, first_row);

//    tree_norm = treeNorm(x, n, ind_work, ind_work.n_rows);

    funVal(iterStep) = arma::as_scalar(Axy.t() * Axy) + lambda1 * sum(abs(x)) + lambda2 *penalty2[0];

    if (bFlag) {break;}
    double norm_xp, norm_xxp;
    norm_xxp = sqrt(arma::as_scalar(xxp.t() * xxp));

    //sgLogisticR.m:403-435
    switch (opts_tFlag)
    {
	  case 0:
        if (iterStep >=2)
        {
	      if (std::abs(funVal(iterStep) - funVal(iterStep - 1)) <= opts_tol * funVal(iterStep - 1))
	      {
	        bFlag = 1;
	      }
	    }
	    break;
	  case 1:
	    if (iterStep >=2)
	    {
		  if (abs(funVal(iterStep) - funVal(iterStep - 1)) <= opts_tol * funVal(iterStep - 1))
	      {
			bFlag = 1;
	      }
		}
		break;
	  case 2:
	    if (funVal(iterStep) <= opts_tol)
	    {
		  bFlag = 1;
	    }
	    break;
	  case 3:
	    norm_xxp = sqrt(arma::as_scalar(xxp.t() * xxp));
	    if (norm_xxp <= opts_tol)
	    {
		  bFlag = 1;
	    }
	    break;
	  case 4:
	    norm_xp = sqrt(arma::as_scalar(xp.t() * xp)); norm_xxp = sqrt(arma::as_scalar(xxp.t() * xxp));
	    if (norm_xxp <= opts_tol * std::max(norm_xp, 1.0))
	    {
		  bFlag = 1;
	    }
	    break;
	  case 5:
        if (iterStep >= opts_maxIter)
        {
          bFlag = 1;
        }
        break;
    }
    if (bFlag) {break;}

	//sgLogisticR.m:438-441
	if (iterStep % opts_rStartNum == 0)
	{
	  alphap = 0;   alpha = 1;   xp = x;   Axp = Ax;   xxp = n_zeros;   L = L/2;
	}

  }
  arma::rowvec x_row = x.t();

  parameters = x_row;

  free(gap);
  return x_row;



  // Here we add the row of ones to the features.
  // The intercept is not penalized. Add an "all ones" row to design and set
  // intercept = false to get a penalized intercept.
//  if (intercept)
//  {
//    p.insert_rows(0, arma::ones<arma::mat>(1, nCols));
//  }

//  if (weights.n_elem > 0)
//  {
//    p = p * diagmat(sqrt(weights));
//    r = sqrt(weights) % responses;
//  }

  // Convert to this form:
  // a * (X X^T) = y X^T.
  // Then we'll use Armadillo to solve it.
  // The total runtime of this should be O(d^2 N) + O(d^3) + O(dN).
  // (assuming the SVD is used to solve it)
//  arma::mat cov = p * p.t() +
//      lambda * arma::eye<arma::mat>(p.n_rows, p.n_rows);

//  parameters = arma::solve(cov, p * r.t());
//  return ComputeError(features, responses);

}

void OLSGLassoLeastR::Predict(const arma::mat& points,
    arma::rowvec& predictions) const
{
  if (intercept)
  {
    // We want to be sure we have the correct number of dimensions in the
    // dataset.
    Log::Assert(points.n_rows == parameters.n_rows - 1);
    // Get the predictions, but this ignores the intercept value
    // (parameters[0]).
    predictions = arma::trans(parameters.subvec(1, parameters.n_elem - 1))
        * points;
    // Now add the intercept.
    predictions += parameters(0);
  }
  else
  {
    // We want to be sure we have the correct number of dimensions in
    // the dataset.
    Log::Assert(points.n_rows == parameters.n_rows);
    predictions = arma::trans(parameters) * points;
  }
}

double OLSGLassoLeastR::ComputeError(const arma::mat& features,
                                      const arma::rowvec& responses) const
{
  // Get the number of columns and rows of the dataset.
  const size_t nCols = features.n_cols;
  const size_t nRows = features.n_rows;

  // Calculate the differences between actual responses and predicted responses.
  // We must also add the intercept (parameters(0)) to the predictions.
  arma::rowvec temp;
  if (intercept)
  {
    // Ensure that we have the correct number of dimensions in the dataset.
    if (nRows != parameters.n_rows - 1)
    {
      Log::Fatal << "The test data must have the same number of columns as the "
          "training file." << std::endl;
    }
    temp = responses - (parameters(0) +
        arma::trans(parameters.subvec(1, parameters.n_elem - 1)) * features);
  }
  else
  {
    // Ensure that we have the correct number of dimensions in the dataset.
    if (nRows != parameters.n_rows)
    {
      Log::Fatal << "The test data must have the same number of columns as the "
          "training file." << std::endl;
    }
    temp = responses - arma::trans(parameters) * features;
  }
  const double cost = arma::dot(temp, temp) / nCols;

  return cost;
}

const arma::colvec OLSGLassoLeastR::altra(const arma::colvec& v_in,
                            const int n,
                            const arma::mat& ind_mat,
                            const int nodes) const
{
	double *x;
	x = (double*) malloc(n*sizeof(double));
	const double* v = v_in.memptr();
    int i, j, m;
    double lambda,twoNorm, ratio;
    double ind[ind_mat.n_cols * ind_mat.n_rows];

    for(int k=0;k<ind_mat.n_cols;k++)
    {
	    for(int l=0;l<ind_mat.n_rows;l++)
	    {
		    ind[(l*3)+k] = ind_mat(l,k);
		}
	}

//	for(int k=0;k<ind_mat.n_cols*ind_mat.n_rows;k++)
//	{
//	    std::cout << "ind[" << k << "]:" << ind[k] << std::endl;
//	}


//std::cout << "Altra 2..." << std::endl;
    /*
     * test whether the first node is special
     */
    if ((int) ind[0]==-1){

        /*
         *Recheck whether ind[1] equals to zero
         */
        if ((int) ind[1]!=-1){
            printf("\n Error! \n Check ind");
            exit(1);
        }

        lambda=ind[2];

        for(j=0;j<n;j++){
            if (v[j]>lambda)
                x[j]=v[j]-lambda;
            else
                if (v[j]<-lambda)
                    x[j]=v[j]+lambda;
                else
                    x[j]=0;
        }

        i=1;
    }
    else{
        memcpy(x, v, sizeof(double) * n);
        i=0;
    }
//std::cout << "Altra 3..." << std::endl;
    /*
     * sequentially process each node
     *
     */
	for(;i < nodes; i++){
        /*
         * compute the L2 norm of this group
         */
		twoNorm=0;
		for(j=(int) ind[3*i]-1;j< (int) ind[3*i+1];j++)
			twoNorm += x[j] * x[j];
        twoNorm=sqrt(twoNorm);

        lambda=ind[3*i+2];
        if (twoNorm>lambda){
            ratio=(twoNorm-lambda)/twoNorm;

            /*
             * shrinkage this group by ratio
             */
            for(j=(int) ind[3*i]-1;j<(int) ind[3*i+1];j++)
                x[j]*=ratio;
        }
        else{
            /*
             * threshold this group to zero
             */
            for(j=(int) ind[3*i]-1;j<(int) ind[3*i+1];j++)
                x[j]=0;
        }
	}
//std::cout << "Altra 4..." << std::endl;
	arma::colvec x_col(&x[0], n);
	free(x);
//std::cout << "Altra 5..." << std::endl;
	return x_col;
}


const double OLSGLassoLeastR::treeNorm(const arma::rowvec& x_in,
                            const int n,
                            const arma::mat& ind_mat,
                            const int nodes) const
{
	double tree_norm;
	const double* x = x_in.memptr();
    int i, j, m;
    double twoNorm, lambda;
    double ind[ind_mat.n_cols * ind_mat.n_rows];

    for(int k=0;k<ind_mat.n_cols;k++)
    {
	    for(int l=0;l<ind_mat.n_rows;l++)
	    {
		    ind[(l*3)+k] = ind_mat(l,k);
		}
	}

//	for(int k=0;k<ind_mat.n_cols*ind_mat.n_rows;k++)
//	{
//	    std::cout << "ind[" << k << "]:" << ind[k] << std::endl;
//	}

    tree_norm=0;

    /*
     * test whether the first node is special
     */
    if ((int) ind[0]==-1){

        /*
         *Recheck whether ind[1] equals to zero
         */
        if ((int) ind[1]!=-1){
            printf("\n Error! \n Check ind");
            exit(1);
        }

        lambda=ind[2];

        for(j=0;j<n;j++){
            tree_norm+=fabs(x[j]);
        }

        tree_norm=tree_norm * lambda;

        i=1;
    }
    else{
        i=0;
    }

    /*
     * sequentially process each node
     *
     */
	for(;i < nodes; i++){
        /*
         * compute the L2 norm of this group
         */
		twoNorm=0;
		for(j=(int) ind[3*i]-1;j< (int) ind[3*i+1];j++)
			twoNorm += x[j] * x[j];
        twoNorm=sqrt(twoNorm);

        lambda=ind[3*i+2];

        tree_norm=tree_norm + lambda*twoNorm;
	}

	return tree_norm;
}


const double OLSGLassoLeastR::computeLambda2Max(const arma::rowvec& x_in,
                            const int n,
                            const arma::mat& ind_mat,
                            const int nodes) const
{
    int i, j, m;
    double lambda,twoNorm;
    const double* x = x_in.memptr();
    double ind[ind_mat.n_cols * ind_mat.n_rows];

    double lambda2_max = 0;

    for(int k=0;k<ind_mat.n_cols;k++)
    {
	    for(int l=0;l<ind_mat.n_rows;l++)
	    {
		    ind[(l*3)+k] = ind_mat(l,k);
		}
	}

    for(i=0;i < nodes; i++){
        /*
         * compute the L2 norm of this group
         */
		twoNorm=0;
		for(j=(int) ind[3*i]-1;j< (int) ind[3*i+1];j++)
			twoNorm += x[j] * x[j];
        twoNorm=sqrt(twoNorm);

        twoNorm=twoNorm/ind[3*i+2];

        if (twoNorm >lambda2_max )
            lambda2_max=twoNorm;
	}

	return lambda2_max;
}




