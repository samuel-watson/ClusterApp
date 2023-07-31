#ifndef MATRIXW_HPP
#define MATRIXW_HPP

#include "general.h"
#include "modelbits.hpp"
#include "openmpheader.h"
#include "maths.h"

namespace glmmr {

using namespace Eigen;

class MatrixW{
public:
  bool attenuated = false;
  VectorXd W_ = VectorXd::Constant(1,1.0);
  glmmr::ModelBits& model;
  MatrixW(glmmr::ModelBits& model_): model(model_) { update(); };
  VectorXd W();
  void update();
};

}

inline VectorXd glmmr::MatrixW::W(){
  return W_;
}

inline void glmmr::MatrixW::update(){
  if(W_.size() != model.n())W_.conservativeResize(model.n());
  ArrayXd nvar_par(model.n());
  ArrayXd xb(model.n());
  if(model.family.family=="gaussian"){
    nvar_par = model.data.variance;
  } else if(model.family.family=="Gamma"){
    nvar_par = model.data.variance.inverse();
  } else if(model.family.family=="beta"){
    nvar_par = (1+model.data.variance);
  } else if(model.family.family=="binomial"){
    nvar_par = model.data.variance.inverse();
  } else {
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

#endif