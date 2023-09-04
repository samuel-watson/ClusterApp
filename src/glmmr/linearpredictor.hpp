#pragma once

#include "general.h"
#include "interpreter.h"
#include "calculator.hpp"
#include "formula.hpp"

namespace glmmr {

class LinearPredictor {
public:
  dblvec parameters;
  glmmr::calculator calc;
  MatrixXd Xdata;
  LinearPredictor(glmmr::Formula& form_,
                  const Eigen::ArrayXXd &data_,
                  const strvec& colnames_) :
    Xdata(data_.rows(),1),
    colnames_vec(colnames_),  
    form(form_),
    n_(data_.rows()),
    X_(MatrixXd::Zero(n_,1))
    {
      form.calculate_linear_predictor(calc,data_,colnames_,Xdata);
      P_ = calc.parameter_names.size();
      parameters.resize(P_);
      std::fill(parameters.begin(),parameters.end(),0.0);
      X_.conservativeResize(n_,P_);
      if(!calc.any_nonlinear){
        X_ = calc.jacobian(parameters,Xdata);
      } else {
        X_.setZero();
      }
    };

  LinearPredictor(glmmr::Formula& form_,
             const Eigen::ArrayXXd &data_,
             const strvec& colnames_,
             const dblvec& parameters_) :
    Xdata(data_.rows(),1),
    colnames_vec(colnames_), 
    form(form_),
    n_(data_.rows()),
    X_(MatrixXd::Zero(n_,1))
     {
      form.calculate_linear_predictor(calc,data_,colnames_,Xdata);
      update_parameters(parameters);
      P_ = calc.parameter_names.size();
      X_.conservativeResize(n_,P_);
      X_ = calc.jacobian(parameters,Xdata);
      x_set = true;
    };

  LinearPredictor(glmmr::Formula& form_,
             const Eigen::ArrayXXd &data_,
             const strvec& colnames_,
             const Eigen::ArrayXd& parameters_) :
    Xdata(data_.rows(),1),
    colnames_vec(colnames_), 
    form(form_),
    n_(data_.rows()),
    X_(MatrixXd::Zero(n_,1))
     {
      form.calculate_linear_predictor(calc,data_,colnames_,Xdata);
      update_parameters(parameters);
      P_ = calc.parameter_names.size();
      X_.conservativeResize(n_,P_);
      X_ = calc.jacobian(parameters,Xdata);
      x_set = true;
    };

  void update_parameters(const dblvec& parameters_){
    parameters = parameters_;
    if(!x_set){
      X_ = calc.jacobian(parameters,Xdata);
      x_set = true;
    }
  };

  void update_parameters(const Eigen::ArrayXd& parameters_){
    dblvec new_parameters(parameters_.data(),parameters_.data()+parameters_.size());
    update_parameters(new_parameters);
  };

  int P(){
    return P_;
  }
  
  int n(){
    return n_;
  }
  
  strvec colnames(){
    return colnames_vec;
  }

  VectorXd xb(){
    VectorXd xb(n());
    if(calc.any_nonlinear){
      xb = calc.linear_predictor(parameters,Xdata);
    } else {
      Map<VectorXd> beta(parameters.data(),parameters.size());
      xb = X_ * beta;
    }
    
    return xb;
  }

  MatrixXd X(){
    if(calc.any_nonlinear){
      X_ = calc.jacobian(parameters,Xdata);
    }
    return X_;
  }
  
  strvec parameter_names(){
    return calc.parameter_names;
  }
  
  VectorXd parameter_vector(){
    VectorXd pars = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(parameters.data(),parameters.size());
    return pars;
  }
  
  bool any_nonlinear(){
    return calc.any_nonlinear;
  }
  
  VectorXd predict_xb(const ArrayXXd& newdata_,
             const ArrayXd& newoffset_){
    glmmr::LinearPredictor newlinpred(form,
                                       newdata_,
                                       colnames(),
                                       parameters);
    VectorXd xb = newlinpred.xb() + newoffset_.matrix();
    return xb;
  }

private:
  strvec colnames_vec;
  glmmr::Formula& form;
  int P_;
  int n_;
  intvec x_cols;
  MatrixXd X_;
  bool x_set = false;
};
}
