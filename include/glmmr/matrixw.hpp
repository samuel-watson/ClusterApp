#pragma once

#include "general.h"
#include "modelbits.hpp"
#include "openmpheader.h"
#include "maths.h"

namespace glmmr {

using namespace Eigen;

template<typename modeltype>
class MatrixW{
public:
  bool        attenuated = false;
  VectorXd    W_ = VectorXd::Constant(1,1.0);
  modeltype&  model;
  MatrixW(modeltype& model_): model(model_) { update(); };
  VectorXd    W() const;
  void        update();
};
}

template<typename modeltype>
inline VectorXd glmmr::MatrixW<modeltype>::W() const{
  return W_;
}

template<typename modeltype>
inline void glmmr::MatrixW<modeltype>::update(){
  if(W_.size() != model.n())W_.conservativeResize(model.n());
  ArrayXd nvar_par(model.n());
  ArrayXd xb(model.n());
  switch(model.family.family){
  case Fam::gaussian:
    nvar_par = model.data.variance;
    break;
  case Fam::gamma:
    nvar_par = model.data.variance.inverse();
    break;
  case Fam::beta:
    nvar_par = (1+model.data.variance);
    break;
  case Fam::binomial:
    nvar_par = model.data.variance.inverse();
    break;
  default:
    nvar_par.setConstant(1.0);
  }
  
  if(attenuated){
    xb = glmmr::maths::attenuted_xb(model.xb(),model.covariance.Z(),model.covariance.D(),model.family.link);
  } else {
    xb = model.xb();
  }
  W_ = glmmr::maths::dhdmu(xb,model.family);
  W_ = (W_.array()*nvar_par).matrix();
  W_ = ((W_.array().inverse()) * model.data.weights).matrix();
}
