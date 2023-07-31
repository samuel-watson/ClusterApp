#ifndef RANDOMEFFECTS_HPP
#define RANDOMEFFECTS_HPP

#include "general.h"
#include "covariance.hpp"
#include "modelbits.hpp"
#include "maths.h"
#include "sparseext.h"
#include "calculator.hpp"

namespace glmmr {

using namespace Eigen;

class RandomEffects{
public:
  sparse ZL;
  MatrixXd u_;
  MatrixXd zu_;
  glmmr::ModelBits& model;
  RandomEffects(glmmr::ModelBits& model_) : 
    ZL(model_.n(),model_.covariance.Q()),
    u_(MatrixXd::Zero(model_.covariance.Q(),1)),
    zu_(model_.n(),1), model(model_) { if(model.covariance.parameters_.size()>0)ZL = model.covariance.ZL_sparse();};
  MatrixXd Zu(){return zu_;};
  MatrixXd u(bool scaled = true);
  vector_matrix predict_re(const ArrayXXd& newdata_,const ArrayXd& newoffset_);
};

}

inline MatrixXd glmmr::RandomEffects::u(bool scaled){
  if(scaled){
    return model.covariance.Lu(u_);
  } else {
    return u_;
  }
}

inline vector_matrix glmmr::RandomEffects::predict_re(const ArrayXXd& newdata_,
                                                      const ArrayXd& newoffset_){
  // generate the merged data
  int nnew = newdata_.rows();
  ArrayXXd mergedata(model.n()+nnew,model.covariance.data_.cols());
  mergedata.topRows(model.n()) = model.covariance.data_;
  mergedata.bottomRows(nnew) = newdata_;
  ArrayXd mergeoffset(model.n()+nnew);
  mergeoffset.head(model.n()) = model.data.offset;
  mergeoffset.tail(nnew) = newoffset_;
  glmmr::Covariance covariancenew(model.formula.formula_,
                                   mergedata,
                                   model.covariance.colnames_,
                                   model.covariance.parameters_);
  glmmr::Covariance covariancenewnew(model.formula.formula_,
                                      newdata_,
                                      model.covariance.colnames_,
                                      model.covariance.parameters_);
  glmmr::LinearPredictor newlinpred(model.formula,
                                     mergedata,
                                     model.linear_predictor.colnames(),
                                     model.linear_predictor.parameters);
  // //generate sigma
  int newQ = covariancenewnew.Q();
  vector_matrix result(newQ);
  result.vec.setZero();
  result.mat.setZero();
  MatrixXd D = covariancenew.D(false,false);
  result.mat = D.block(model.covariance.Q(),model.covariance.Q(),newQ,newQ);
  MatrixXd D22 = D.block(0,0,model.covariance.Q(),model.covariance.Q());
  D22 = D22.llt().solve(MatrixXd::Identity(model.covariance.Q(),model.covariance.Q()));
  MatrixXd D12 = D.block(model.covariance.Q(),0,newQ,model.covariance.Q());
  MatrixXd Lu = model.covariance.Lu(u(false));
  MatrixXd SSV = D12 * D22 * Lu;
  result.vec = SSV.rowwise().mean();
  MatrixXd D121 = D12 * D22 * D12.transpose();
  result.mat -= D12 * D22 * D12.transpose();
  return result;
}

#endif