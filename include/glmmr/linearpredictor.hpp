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
  glmmr::Formula& form;
  LinearPredictor(glmmr::Formula& form_,const Eigen::ArrayXXd &data_,const strvec& colnames_);
  LinearPredictor(glmmr::Formula& form_,const Eigen::ArrayXXd &data_,const strvec& colnames_,const dblvec& parameters_);
  LinearPredictor(glmmr::Formula& form_,const Eigen::ArrayXXd &data_,const strvec& colnames_,const Eigen::ArrayXd& parameters_);
  LinearPredictor(const glmmr::LinearPredictor& linpred);
  virtual void update_parameters(const dblvec& parameters_);
  virtual void update_parameters(const Eigen::ArrayXd& parameters_);
  int P() const;
  int n() const;
  strvec colnames() const;
  virtual VectorXd xb();
  virtual MatrixXd X();
  strvec parameter_names() const;
  VectorXd parameter_vector();
  bool any_nonlinear() const;
  virtual VectorXd predict_xb(const ArrayXXd& newdata_,const ArrayXd& newoffset_);
protected:
  strvec colnames_vec;
  int P_;
  int n_;
  intvec x_cols;
  MatrixXd X_;
  bool x_set = false;
};
}

inline glmmr::LinearPredictor::LinearPredictor(glmmr::Formula& form_,
                const Eigen::ArrayXXd &data_,
                const strvec& colnames_) :
  Xdata(data_.rows(),1),
  form(form_),
  colnames_vec(colnames_),  
  n_(data_.rows()),
  X_(MatrixXd::Zero(n_,1))
{
  form.calculate_linear_predictor(calc,data_,colnames_,Xdata);
  P_ = calc.parameter_names.size();
  parameters.resize(P_);
  if(!calc.any_nonlinear){
    std::fill(parameters.begin(),parameters.end(),0.0);
  } else {
    // to avoid problems like dividing by zero
    std::fill(parameters.begin(),parameters.end(),1.0);
  }
  X_.conservativeResize(n_,P_);
  if(!calc.any_nonlinear){
    X_ = calc.jacobian(parameters,Xdata);
  } else {
    X_.setZero();
  }
};

inline glmmr::LinearPredictor::LinearPredictor(glmmr::Formula& form_,
                const Eigen::ArrayXXd &data_,
                const strvec& colnames_,
                const dblvec& parameters_) :
  Xdata(data_.rows(),1),
  form(form_),
  colnames_vec(colnames_), 
  n_(data_.rows()),
  X_(MatrixXd::Zero(n_,1))
{
  form.calculate_linear_predictor(calc,data_,colnames_,Xdata);
  P_ = calc.parameter_names.size();
  update_parameters(parameters);
  X_.conservativeResize(n_,P_);
  X_ = calc.jacobian(parameters,Xdata);
  x_set = true;
};

inline glmmr::LinearPredictor::LinearPredictor(glmmr::Formula& form_,
                const Eigen::ArrayXXd &data_,
                const strvec& colnames_,
                const Eigen::ArrayXd& parameters_) :
  Xdata(data_.rows(),1),
  form(form_),
  colnames_vec(colnames_), 
  n_(data_.rows()),
  X_(MatrixXd::Zero(n_,1))
{
  form.calculate_linear_predictor(calc,data_,colnames_,Xdata);
  P_ = calc.parameter_names.size();
  update_parameters(parameters);
  X_.conservativeResize(n_,P_);
  X_ = calc.jacobian(parameters,Xdata);
  x_set = true;
};

inline glmmr::LinearPredictor::LinearPredictor(const glmmr::LinearPredictor& linpred) :
  Xdata(linpred.Xdata.rows(),1),
  form(linpred.form),
  colnames_vec(linpred.colnames_vec), 
  n_(linpred.Xdata.rows()),
  X_(MatrixXd::Zero(n_,1))
{
  form.calculate_linear_predictor(calc,linpred.Xdata.array(),linpred.colnames_vec,Xdata);
  P_ = calc.parameter_names.size();
  update_parameters(linpred.parameters);
  X_.conservativeResize(n_,P_);
  X_ = calc.jacobian(parameters,Xdata);
  x_set = true;
};

inline void glmmr::LinearPredictor::update_parameters(const dblvec& parameters_){
#if defined(R_BUILD) && defined(ENABLE_DEBUG)
  Rcpp::Rcout << "\nUpdating parameters... Old parameters: ";
  for(const auto& i: parameters)Rcpp::Rcout << i << " ";
  Rcpp::Rcout << "\nNew parameters: ";
  for(const auto& i: parameters_)Rcpp::Rcout << i << " ";
#endif
  
  #ifdef R_BUILD
  if(parameters_.size()!=(unsigned)P())Rcpp::stop(std::to_string(parameters_.size())+" parameters provided, "+std::to_string(P())+" required");
  if(parameters_.size()!=(unsigned)calc.parameter_count)Rcpp::stop(std::to_string(parameters_.size())+" parameters provided, "+std::to_string(calc.parameter_count)+" required");
  #endif
  
  
  parameters = parameters_;
  if(!x_set){
    X_ = calc.jacobian(parameters,Xdata);
    x_set = true;
  }
};

inline void glmmr::LinearPredictor::update_parameters(const Eigen::ArrayXd& parameters_){
  dblvec new_parameters(parameters_.data(),parameters_.data()+parameters_.size());
  update_parameters(new_parameters);
};

inline int glmmr::LinearPredictor::P() const{
  return P_;
}

inline int glmmr::LinearPredictor::n() const{
  return n_;
}

inline strvec glmmr::LinearPredictor::colnames() const{
  return colnames_vec;
}

inline VectorXd glmmr::LinearPredictor::xb(){
  VectorXd xb(n());
  if(calc.any_nonlinear){
    xb = calc.linear_predictor(parameters,Xdata);
  } else {
    Map<VectorXd> beta(parameters.data(),parameters.size());
    xb = X_ * beta;
  }
  return xb;
}

inline MatrixXd glmmr::LinearPredictor::X(){
  if(calc.any_nonlinear){
    X_ = calc.jacobian(parameters,Xdata);
  }
  return X_;
}

inline strvec glmmr::LinearPredictor::parameter_names() const{
  return calc.parameter_names;
}

inline VectorXd glmmr::LinearPredictor::parameter_vector(){
  VectorXd pars = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(parameters.data(),parameters.size());
  return pars;
}

inline bool glmmr::LinearPredictor::any_nonlinear() const{
  return calc.any_nonlinear;
}

inline VectorXd glmmr::LinearPredictor::predict_xb(const ArrayXXd& newdata_,
                    const ArrayXd& newoffset_) {
  LinearPredictor newlinpred(form,newdata_,colnames());
  newlinpred.update_parameters(parameters);
  VectorXd xb = newlinpred.xb() + newoffset_.matrix();
  return xb;
}




