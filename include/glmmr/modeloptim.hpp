#pragma once

#include "general.h"
#include "modelbits.hpp"
#include "randomeffects.hpp"
#include "modelmatrix.hpp"
#include "openmpheader.h"
#include "maths.h"
#include "algo.h"
#include "sparse.h"
#include "calculator.hpp"
#include "optim/optim.h"

namespace glmmr {

using namespace Eigen;

template<typename modeltype>
class ModelOptim {
  
public:
  // members 
  modeltype&                        model;
  glmmr::ModelMatrix<modeltype>&    matrix;
  glmmr::RandomEffects<modeltype>&  re;
  int                               trace = 0;
  // ArrayXXd                          ll_previous; // log likelihood values for all u samples
  ArrayXXd                          ll_current;
  std::pair<double,double>          current_ll_values = {0.0,0.0};
  std::pair<double,double>          previous_ll_values = {0.0,0.0};
  std::pair<double,double>          current_ll_var = {0.0,0.0};
  std::pair<double,double>          previous_ll_var = {0.0,0.0};
  std::pair<int,int>                fn_counter = {0,0};
  
  // constructor
  ModelOptim(modeltype& model_, glmmr::ModelMatrix<modeltype>& matrix_,glmmr::RandomEffects<modeltype>& re_) ;

  // control parameters for the optimisers - direct will be removed as its useless.
  struct OptimControl {
    int     npt = 0;
    double  rhobeg = 0;
    double  rhoend = 0;
    bool    direct = false;
    double  direct_range_beta = 3.0;
    int     max_iter_direct = 100; 
    double  epsilon = 1e-4; 
    bool    select_one = true; 
    bool    trisect_once = false; 
    int     max_eval = 0; 
    bool    mrdirect = false;
    double  g_epsilon = 1e-8;
    int     past = 3;
    double  delta = 1e-8;
    int     max_linesearch = 64;
    double  alpha = 0.8;
    bool    saem = false;
    bool    pr_average = true;
    bool    reml = true;
  } control;
  
  // functions
  virtual void    update_beta(const dblvec &beta);
  virtual void    update_beta(const VectorXd &beta);
  virtual void    update_theta(const dblvec &theta);
  virtual void    update_theta(const VectorXd &theta);
  virtual void    update_u(const MatrixXd& u_, bool append); // two versions needed so CRAN won't complain with linked package
  virtual void    update_u(const MatrixXd& u_); 
  virtual double  log_likelihood(bool beta);
  virtual double  log_likelihood();
  virtual double  full_log_likelihood();
  virtual void    nr_beta();
  virtual void    laplace_nr_beta_u();
  virtual void    update_var_par(const double& v);
  virtual void    update_var_par(const ArrayXd& v);
  template<class algo, typename = std::enable_if_t<std::is_base_of<optim_algo, algo>::value> >
  void            ml_beta();
  template<class algo, typename = std::enable_if_t<std::is_base_of<optim_algo, algo>::value> >
  void            ml_theta();
  template<class algo, typename = std::enable_if_t<std::is_base_of<optim_algo, algo>::value> >
  void            ml_all();
  template<class algo, typename = std::enable_if_t<std::is_base_of<optim_algo, algo>::value> >
  void            laplace_ml_beta_u();
  template<class algo, typename = std::enable_if_t<std::is_base_of<optim_algo, algo>::value> >
  void            laplace_ml_theta();
  template<class algo, typename = std::enable_if_t<std::is_base_of<optim_algo, algo>::value> >
  void            laplace_ml_beta_theta();
  virtual double  aic();
  virtual ArrayXd optimum_weights(double N, VectorXd C, double tol = 1e-5, int max_iter = 501);
  void            set_bobyqa_control(int npt_, double rhobeg_, double rhoend_);
  void            set_direct_control(bool direct = false, double direct_range_beta = 3.0, int max_iter = 100, double epsilon = 1e-4, bool select_one = true, bool trisect_once = false, 
                                      int max_eval = 0, bool mrdirect = false);
  void            set_lbfgs_control(double g_epsilon = 1e-8, int past = 3, double delta = 1e-8, int max_linesearch = 64);
  void            set_bound(const dblvec& bound, bool lower = true);
  void            set_theta_bound(const dblvec& bound, bool lower = true);
  void            use_reml(bool reml);
  int             P() const;
  int             Q() const;
  double          ll_diff_variance(bool beta = true, bool theta = true);
  std::pair<double,double>  current_likelihood_values();
  std::pair<double,double>  u_diagnostic();
  void            reset_fn_counter();
  void            set_quantile(const double& q);
  // functions to optimise
  double          log_likelihood_beta(const dblvec &beta);
  double          log_likelihood_beta_with_gradient(const VectorXd &beta, VectorXd& g);
  double          log_likelihood_theta(const dblvec &theta);
  double          log_likelihood_theta_with_gradient(const VectorXd& theta, VectorXd& g);
  double          log_likelihood_all(const dblvec &par);
  double          log_likelihood_laplace_beta_u(const dblvec &par);
  double          log_likelihood_laplace_beta_u_with_gradient(const VectorXd& theta, VectorXd& g);
  double          log_likelihood_laplace_theta(const dblvec &par);
  double          log_likelihood_laplace_beta_theta(const dblvec &par);
  
protected:
// objects
  dblvec    lower_bound;
  dblvec    upper_bound; // bounds for beta
  dblvec    lower_bound_theta;
  dblvec    upper_bound_theta; // bounds for beta
  bool      beta_bounded = false;
  double    quantile = 0.5;
  
  // functions
  void            calculate_var_par();
  dblvec          get_start_values(bool beta, bool theta, bool var = true);
  dblvec          get_lower_values(bool beta, bool theta, bool var = true, bool u = false);
  dblvec          get_upper_values(bool beta, bool theta, bool var = true, bool u = false);
  void            set_direct_control(directd& op);
  void            set_bobyqa_control(bobyqad& op);
  void            set_newuoa_control(newuoad& op);
  void            set_lbfgs_control(lbfgsd& op);
  
private:
  // used for REML
  void        generate_czz();
  MatrixXd    CZZ = MatrixXd::Zero(1,1);
};

}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::set_direct_control(directd& op){
    op.control.max_iter = control.max_iter_direct;
    op.control.epsilon = control.epsilon;
    op.control.select_one = control.select_one;
    op.control.trisect_once = control.trisect_once;
    op.control.trace = trace;
    op.control.mrdirect = control.mrdirect;
    op.control.max_eval = control.max_eval;
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::set_bobyqa_control(bobyqad& op){
  op.control.trace = trace;
  op.control.rhobeg = control.rhobeg;
  op.control.rhoend = control.rhoend;
  op.control.npt = control.npt;
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::set_newuoa_control(newuoad& op){
  op.control.trace = trace;
  op.control.rhobeg = control.rhobeg;
  op.control.rhoend = control.rhoend;
  op.control.npt = control.npt;
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::set_lbfgs_control(lbfgsd& op){
  op.control.trace = trace;
  op.control.g_epsilon = control.g_epsilon;
  op.control.past = control.past;
  op.control.delta = control.delta;
  op.control.max_linesearch = control.max_linesearch;
}

template<typename modeltype>
template<class algo, typename>
inline void glmmr::ModelOptim<modeltype>::ml_beta(){
  dblvec start = get_start_values(true,false,false);
  // store previous log likelihood values for convergence calculations
  previous_ll_values.first = current_ll_values.first;
  previous_ll_var.first = current_ll_var.first;
  
  if constexpr (std::is_same_v<algo,LBFGS>){
    VectorXd start_vec = Map<VectorXd>(start.data(),start.size());
    optim<double(const VectorXd&, VectorXd&),algo> op(start_vec);
    set_lbfgs_control(op);
    if(beta_bounded) op.set_bounds(lower_bound,upper_bound);
      if constexpr (std::is_same_v<modeltype,bits>)
    {
      op.template fn<&glmmr::ModelOptim<bits>::log_likelihood_beta_with_gradient, glmmr::ModelOptim<bits> >(this);
    } else if constexpr (std::is_same_v<modeltype,bits_nngp>) {
      op.template fn<&glmmr::ModelOptim<bits_nngp>::log_likelihood_beta_with_gradient, glmmr::ModelOptim<bits_nngp> >(this);
    } else if constexpr (std::is_same_v<modeltype,bits_hsgp>){
      op.template fn<&glmmr::ModelOptim<bits_hsgp>::log_likelihood_beta_with_gradient, glmmr::ModelOptim<bits_hsgp> >(this);
    }
    op.minimise();
  } else {
    optim<double(const std::vector<double>&),algo> op(start);
    if constexpr (std::is_same_v<algo,DIRECT>) {
      op.set_bounds(start,dblvec(start.size(),control.direct_range_beta),true);
      set_direct_control(op);
    } else if constexpr (std::is_same_v<algo,BOBYQA>) {
      set_bobyqa_control(op);
    } else if constexpr (std::is_same_v<algo,NEWUOA>) {
      set_newuoa_control(op);
    }
    if(beta_bounded) op.set_bounds(lower_bound,upper_bound);
      if constexpr (std::is_same_v<modeltype,bits>)
    {
      op.template fn<&glmmr::ModelOptim<bits>::log_likelihood_beta, glmmr::ModelOptim<bits> >(this);
    } else if constexpr (std::is_same_v<modeltype,bits_nngp>) {
      op.template fn<&glmmr::ModelOptim<bits_nngp>::log_likelihood_beta, glmmr::ModelOptim<bits_nngp> >(this);
    } else if constexpr (std::is_same_v<modeltype,bits_hsgp>){
      op.template fn<&glmmr::ModelOptim<bits_hsgp>::log_likelihood_beta, glmmr::ModelOptim<bits_hsgp> >(this);
    }
    op.minimise();
  }
  int eval_size = control.saem ? re.mcmc_block_size : ll_current.rows();
  current_ll_values.first = ll_current.col(0).tail(eval_size).mean();
  current_ll_var.first = (ll_current.col(0).tail(eval_size) - ll_current.col(0).tail(eval_size).mean()).square().sum() / (eval_size - 1);
}

template<typename modeltype>
template<class algo, typename>
inline void glmmr::ModelOptim<modeltype>::ml_theta(){  
  dblvec start = get_start_values(false,true,false);  
  dblvec lower = get_lower_values(false,true,false);
  dblvec upper = get_upper_values(false,true,false);
  // store previous log likelihood values for convergence calculations
  previous_ll_values.second = current_ll_values.second;
  previous_ll_var.second = current_ll_var.second;
  if(re.scaled_u_.cols() != re.u_.cols())re.scaled_u_.resize(NoChange,re.u_.cols());
  re.scaled_u_ = model.covariance.Lu(re.u_);  
  if(control.reml) generate_czz();
  // optimisation
  if constexpr (std::is_same_v<algo,LBFGS>){
    VectorXd start_vec = Map<VectorXd>(start.data(),start.size());
    optim<double(const VectorXd&, VectorXd&),algo> op(start_vec); 
    op.set_bounds(lower,upper);
    set_lbfgs_control(op);
    if constexpr (std::is_same_v<modeltype,bits>)
    {
      op.template fn<&glmmr::ModelOptim<bits>::log_likelihood_theta_with_gradient, glmmr::ModelOptim<bits> >(this);
    } else if constexpr (std::is_same_v<modeltype,bits_nngp>) {
      op.template fn<&glmmr::ModelOptim<bits_nngp>::log_likelihood_theta_with_gradient, glmmr::ModelOptim<bits_nngp> >(this);
    } else if constexpr (std::is_same_v<modeltype,bits_hsgp>){
      op.template fn<&glmmr::ModelOptim<bits_hsgp>::log_likelihood_theta_with_gradient, glmmr::ModelOptim<bits_hsgp> >(this);
    }
    op.minimise();
  } else {
    optim<double(const std::vector<double>&),algo> op(start);
    if constexpr (std::is_same_v<algo,DIRECT>) {      
      dblvec upper2(lower.size());
      std::fill(upper2.begin(),upper2.end(),1.0);
      op.set_bounds(lower,upper2,false);
      set_direct_control(op);
    } else if constexpr (std::is_same_v<algo,BOBYQA>) {
      set_bobyqa_control(op);
      op.set_bounds(lower,upper);
    } else if constexpr (std::is_same_v<algo,NEWUOA>) {
      set_newuoa_control(op);
      op.set_bounds(lower,upper);
    }
      if constexpr (std::is_same_v<modeltype,bits>)
    {
      op.template fn<&glmmr::ModelOptim<bits>::log_likelihood_theta, glmmr::ModelOptim<bits> >(this);
    } else if constexpr (std::is_same_v<modeltype,bits_nngp>) {
      op.template fn<&glmmr::ModelOptim<bits_nngp>::log_likelihood_theta, glmmr::ModelOptim<bits_nngp> >(this);
    } else if constexpr (std::is_same_v<modeltype,bits_hsgp>){
      op.template fn<&glmmr::ModelOptim<bits_hsgp>::log_likelihood_theta, glmmr::ModelOptim<bits_hsgp> >(this);
    }
    op.minimise();
  }
  int eval_size = control.saem ? re.mcmc_block_size : ll_current.rows();
  current_ll_values.second = ll_current.col(1).tail(eval_size).mean();
  current_ll_var.second = (ll_current.col(1).tail(eval_size) - ll_current.col(1).tail(eval_size).mean()).square().sum() / (eval_size - 1);
  calculate_var_par();
}

// THIS FUNCTION BELOW NEEDS UPDATING - IT IS NOT CURRENTLY USED AND LACKS FUNCTIONALITY
// FOR U DIAGNOSTICS 
template<typename modeltype>
template<class algo, typename>
inline void glmmr::ModelOptim<modeltype>::ml_all(){
  if(re.scaled_u_.cols() != re.u_.cols())re.scaled_u_.resize(NoChange,re.u_.cols());
  re.scaled_u_ = model.covariance.Lu(re.u_);
  dblvec start = get_start_values(true,true,false);
  optim<double(const std::vector<double>&),algo> op(start);
  if constexpr (std::is_same_v<algo,DIRECT>) {
    op.set_bounds(start,dblvec(start.size(),control.direct_range_beta),true);
    set_direct_control(op);
  } else if constexpr (std::is_same_v<algo,BOBYQA>) {
    set_bobyqa_control(op);
  } else if constexpr (std::is_same_v<algo,NEWUOA>) {
    set_newuoa_control(op);
  } else if constexpr (std::is_same_v<algo,LBFGS>) {
    throw std::runtime_error("L-BFGS not available for beta & theta optimisation");
  }
  dblvec lower = get_lower_values(true,true,false);
  dblvec upper = get_upper_values(true,true,false);
  op.set_bounds(lower,upper);
  if constexpr (std::is_same_v<modeltype,bits>)
  {
    op.template fn<&glmmr::ModelOptim<bits>::log_likelihood_all, glmmr::ModelOptim<bits> >(this);
  } else if constexpr (std::is_same_v<modeltype,bits_nngp>) {
    op.template fn<&glmmr::ModelOptim<bits_nngp>::log_likelihood_all, glmmr::ModelOptim<bits_nngp> >(this);
  } else if constexpr (std::is_same_v<modeltype,bits_hsgp>){
    op.template fn<&glmmr::ModelOptim<bits_hsgp>::log_likelihood_all, glmmr::ModelOptim<bits_hsgp> >(this);
  }
  op.minimise();
  calculate_var_par();
}

template<typename modeltype>
template<class algo, typename>
inline void glmmr::ModelOptim<modeltype>::laplace_ml_beta_u(){
  dblvec start = get_start_values(true,false,false);
  for(int i = 0; i< Q(); i++) start.push_back(re.u_(i,0));
  if constexpr (std::is_same_v<algo,LBFGS>){
    VectorXd start_vec = Map<VectorXd>(start.data(),start.size());
    optim<double(const VectorXd&, VectorXd&),algo> op(start_vec);
    set_lbfgs_control(op);
    if(lower_bound.size()==P())
    {
      dblvec lower = get_lower_values(true,false,false,true);
      dblvec upper = get_upper_values(true,false,false,true);
      op.set_bounds(lower,upper);
    }
      if constexpr (std::is_same_v<modeltype,bits>)
    {
      op.template fn<&glmmr::ModelOptim<bits>::log_likelihood_laplace_beta_u_with_gradient, glmmr::ModelOptim<bits> >(this);
    } else if constexpr (std::is_same_v<modeltype,bits_nngp>) {
      op.template fn<&glmmr::ModelOptim<bits_nngp>::log_likelihood_laplace_beta_u_with_gradient, glmmr::ModelOptim<bits_nngp> >(this);
    } else if constexpr (std::is_same_v<modeltype,bits_hsgp>){
      op.template fn<&glmmr::ModelOptim<bits_hsgp>::log_likelihood_laplace_beta_u_with_gradient, glmmr::ModelOptim<bits_hsgp> >(this);
    }
    op.minimise();
  } else {
    optim<double(const std::vector<double>&),algo> op(start);
    if constexpr (std::is_same_v<algo,DIRECT>) {
      op.set_bounds(start,dblvec(start.size(),control.direct_range_beta),true);
      set_direct_control(op);
    } else if constexpr (std::is_same_v<algo,BOBYQA>) {
      set_bobyqa_control(op);
    } else if constexpr (std::is_same_v<algo,NEWUOA>) {
      set_newuoa_control(op);
    }
    if(lower_bound.size()==P())
    {
      dblvec lower = get_lower_values(true,false,false,true);
      dblvec upper = get_upper_values(true,false,false,true);
      op.set_bounds(lower,upper);
    }
    if constexpr (std::is_same_v<modeltype,bits>)
    {
      op.template fn<&glmmr::ModelOptim<bits>::log_likelihood_laplace_beta_u, glmmr::ModelOptim<bits> >(this);
    } else if constexpr (std::is_same_v<modeltype,bits_nngp>) {
      op.template fn<&glmmr::ModelOptim<bits_nngp>::log_likelihood_laplace_beta_u, glmmr::ModelOptim<bits_nngp> >(this);
    } else if constexpr (std::is_same_v<modeltype,bits_hsgp>){
      op.template fn<&glmmr::ModelOptim<bits_hsgp>::log_likelihood_laplace_beta_u, glmmr::ModelOptim<bits_hsgp> >(this);
    }
    op.minimise();
  }
  calculate_var_par();
}

template<typename modeltype>
template<class algo, typename>
inline void glmmr::ModelOptim<modeltype>::laplace_ml_theta()
{
  
  dblvec start = get_start_values(false,true,false);  
  dblvec lower = get_lower_values(false,true,false);
  dblvec upper = get_upper_values(false,true,false);
  
  if(re.scaled_u_.cols() != re.u_.cols())re.scaled_u_.conservativeResize(NoChange,re.u_.cols());
  re.scaled_u_ = model.covariance.Lu(re.u_);
  if(control.reml) generate_czz();
  if constexpr (std::is_same_v<algo,LBFGS>){
    VectorXd start_vec = Map<VectorXd>(start.data(),start.size());
    optim<double(const VectorXd&, VectorXd&),algo> op(start_vec);
    op.set_bounds(lower,upper);
    set_lbfgs_control(op);
      if constexpr (std::is_same_v<modeltype,bits>)
    {
      op.template fn<&glmmr::ModelOptim<bits>::log_likelihood_theta_with_gradient, glmmr::ModelOptim<bits> >(this);
    } else {
      throw std::runtime_error("L-BFGS not available for approximate covariance");
    }
    op.minimise();
  } else {
     optim<double(const std::vector<double>&),algo> op(start);
    if constexpr (std::is_same_v<algo,DIRECT>) {
      op.set_bounds(start,dblvec(start.size(),1.0),true);
      set_direct_control(op);
    } else if constexpr (std::is_same_v<algo,BOBYQA>) {
      set_bobyqa_control(op);
    } else if constexpr (std::is_same_v<algo,NEWUOA>) {
      set_newuoa_control(op);
    } else if constexpr (std::is_same_v<algo,LBFGS>) {
      throw std::runtime_error("L-BFGS not available for Laplace theta optimisation");
    }
    op.set_bounds(lower,upper);
    if constexpr (std::is_same_v<modeltype,bits>)
    {
      op.template fn<&glmmr::ModelOptim<bits>::log_likelihood_laplace_theta, glmmr::ModelOptim<bits> >(this);
    } else if constexpr (std::is_same_v<modeltype,bits_nngp>) {
      op.template fn<&glmmr::ModelOptim<bits_nngp>::log_likelihood_laplace_theta, glmmr::ModelOptim<bits_nngp> >(this);
    } else if constexpr (std::is_same_v<modeltype,bits_hsgp>){
      op.template fn<&glmmr::ModelOptim<bits_hsgp>::log_likelihood_laplace_theta, glmmr::ModelOptim<bits_hsgp> >(this);
    }
    op.minimise();
  }

 
}

template<typename modeltype>
template<class algo, typename>
inline void glmmr::ModelOptim<modeltype>::laplace_ml_beta_theta(){
  if(re.scaled_u_.cols() != re.u_.cols())re.scaled_u_.conservativeResize(NoChange,re.u_.cols());
  re.scaled_u_ = model.covariance.Lu(re.u_);
  dblvec start = get_start_values(true,true,false);  
  dblvec lower = get_lower_values(true,true,false);
  dblvec upper = get_upper_values(true,true,false);
  optim<double(const std::vector<double>&),algo> op(start);
  if constexpr (std::is_same_v<algo,DIRECT>) {
    op.set_bounds(start,dblvec(start.size(),control.direct_range_beta),true);
    set_direct_control(op);
  } else if constexpr (std::is_same_v<algo,BOBYQA>) {
    set_bobyqa_control(op);
  } else if constexpr (std::is_same_v<algo,NEWUOA>) {
    set_newuoa_control(op);
  } else if constexpr (std::is_same_v<algo,LBFGS>) {
    throw std::runtime_error("L-BFGS not available for Laplace beta-theta optimisation");
  }
  op.set_bounds(lower,upper);
  if constexpr (std::is_same_v<modeltype,bits>)
  {
    op.template fn<&glmmr::ModelOptim<bits>::log_likelihood_laplace_beta_theta, glmmr::ModelOptim<bits> >(this);
  } else if constexpr (std::is_same_v<modeltype,bits_nngp>) {
    op.template fn<&glmmr::ModelOptim<bits_nngp>::log_likelihood_laplace_beta_theta, glmmr::ModelOptim<bits_nngp> >(this);
  } else if constexpr (std::is_same_v<modeltype,bits_hsgp>){
    op.template fn<&glmmr::ModelOptim<bits_hsgp>::log_likelihood_laplace_beta_theta, glmmr::ModelOptim<bits_hsgp> >(this);
  }
  op.minimise();
  calculate_var_par();
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::set_direct_control(bool direct, double direct_range_beta, int max_iter, double epsilon, bool select_one, bool trisect_once, 
                                      int max_eval, bool mrdirect)
                                      {
  control.direct = direct;
  control.max_iter_direct = max_iter;
  control.epsilon = epsilon;
  control.select_one = select_one;
  control.trisect_once = trisect_once;
  control.max_eval = max_eval;
  control.mrdirect = mrdirect;
  control.direct_range_beta = direct_range_beta;
 }

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::reset_fn_counter()
{
  fn_counter.first = 0;
  fn_counter.second = 0;
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::set_quantile(const double& q)
{
  if(q <= 0 || q >= 1)throw std::runtime_error("q !in [0,1]");
  quantile = q;
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::ll_diff_variance(bool beta, bool theta)
{
  
  // int eval_size = ll_previous.rows() < ll_current.rows() ? ll_previous.rows() : ll_current.rows();
  // 
  // ArrayXd diff = ArrayXd::Zero(eval_size); 
  // if(beta) diff += ll_current.col(0).head(eval_size) - ll_previous.col(0).head(eval_size);
  // if(theta) diff += ll_current.col(1).head(eval_size) - ll_previous.col(1).head(eval_size);
  double var = 0;
  if(beta) var += current_ll_var.first + previous_ll_var.first;
  if(theta) var += current_ll_var.second + previous_ll_var.second;
  return var ; //(diff - diff.mean()).square().sum() / (diff.size() - 1);
}

template<typename modeltype>
inline std::pair<double,double> glmmr::ModelOptim<modeltype>::current_likelihood_values()
{
  return current_ll_values;
}

template<typename modeltype>
inline std::pair<double,double> glmmr::ModelOptim<modeltype>::u_diagnostic()
{
  std::pair<double, double> ll;
  ll.first = current_ll_values.first - previous_ll_values.first;
  ll.second = current_ll_values.second - previous_ll_values.second;
  return ll;
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::set_lbfgs_control(double g_epsilon, int past, double delta, int max_linesearch)
{
  control.g_epsilon = g_epsilon;
  control.past = past;
  control.delta = delta;
  control.max_linesearch = max_linesearch;
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::log_likelihood_all(const dblvec &par)
{
  int G = model.covariance.npar();
  auto first = par.begin();
  auto last1 = par.begin() + P();
  auto last2 = par.begin() + P() + G;
  dblvec beta(first,last1);
  dblvec theta(last1,last2);
  model.linear_predictor.update_parameters(beta);
  update_theta(theta);
  if(model.family.family==Fam::gaussian || model.family.family==Fam::gamma || model.family.family==Fam::beta)update_var_par(par[P()+G]);
  double ll = full_log_likelihood();
  return -1.0*ll;
  // need to fix importance sampling - this isn't right
  // if(importance){
  //   return -1.0 * log(exp(ll)/ exp(denomD));
  // } else {
  //   return -1.0*ll;
  // }
}


template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::log_likelihood_laplace_beta_u(const dblvec &par)
{
  auto start = par.begin();
  auto end = par.begin() + P();
  dblvec beta(start,end);
  MatrixXd v(Q(),1);
  for(int i = 0; i < Q(); i++)v(i,0) = par[P() + i];
  model.linear_predictor.update_parameters(beta);
  update_u(v);
  double logl = v.col(0).transpose()*v.col(0);
  double ll = log_likelihood();
  matrix.W.update();
  MatrixXd LZWZL = model.covariance.LZWZL(matrix.W.W());
  double LZWdet = glmmr::maths::logdet(LZWZL);
  return -1.0*(ll - 0.5*logl - 0.5*LZWdet);
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::log_likelihood_laplace_theta(const dblvec &par)
{
  // update_theta(par);
  // matrix.W.update();
  // double logl = re.u_.col(0).transpose() * re.u_.col(0);
  // double ll = log_likelihood();
  // MatrixXd LZWZL = model.covariance.LZWZL(matrix.W.W());
  // double LZWdet = glmmr::maths::logdet(LZWZL);
  // return -1*(ll - 0.5*logl - 0.5*LZWdet);
  
  model.covariance.update_parameters(par);
  double ll =  model.covariance.log_likelihood(re.scaled_u_.col(0));
  
  if(control.reml){
    // REML correction to the log-likelihood
    MatrixXd D = model.covariance.D().llt().solve(MatrixXd::Identity(Q(),Q()));
    // if(D.rows()!=Q() || D.cols()!=Q()) throw std::runtime_error("D wrong size");
    // if(CZZ.rows()!=Q() || CZZ.cols()!=Q()) throw std::runtime_error("CZZ wrong size");
    MatrixXd CZZD = CZZ + D;
    CZZD = CZZD.llt().solve(MatrixXd::Identity(Q(),Q()));
    double trCZZ = 0;
    for(int i = 0; i < Q(); i++){
      for(int j = 0; j < Q(); j++){
        trCZZ += D(i,j)*CZZD(j,i);
      }
    }
    ll += -0.5*trCZZ;
  }
  
  return -1.0*ll;
  
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::log_likelihood_laplace_beta_theta(const dblvec &par)
{
  auto start = par.begin();
  auto end1 = par.begin() + P();
  auto end2 = par.begin() + P() + model.covariance.npar();
  dblvec beta(start,end1);
  dblvec theta(end1,end2);
  model.linear_predictor.update_parameters(beta);
  update_theta(theta);
  double ll = log_likelihood();
  double logl = re.u_.col(0).transpose() * re.u_.col(0);
  matrix.W.update();
  MatrixXd LZWZL = model.covariance.LZWZL(matrix.W.W());
  double LZWdet = glmmr::maths::logdet(LZWZL);
  return -1*(ll - 0.5*logl - 0.5*LZWdet);
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::log_likelihood_beta(const dblvec& beta)
{
  model.linear_predictor.update_parameters(beta);
  double ll = log_likelihood();
  fn_counter.first += re.scaled_u_.cols();
  if(control.saem)
  {
    int     iteration = std::max((int)re.zu_.cols() / re.mcmc_block_size, 1);
    double  gamma = pow(1.0/iteration,control.alpha);
    double  ll_t = 0;
    double  ll_pr = 0;
    for(int i = 0; i < iteration; i++){
      int lower_range = i * re.mcmc_block_size;
      int upper_range = (i + 1) * re.mcmc_block_size;
      if(i == (iteration - 1) && iteration > 1){
        double ll_t_c = ll_t;
        double ll_pr_c = ll_pr;
        ll_t = ll_t + gamma*(ll_current.col(0).segment(lower_range, re.mcmc_block_size).mean() - ll_t);
        if(control.pr_average) ll_pr += ll_t;
        for(int j = lower_range; j < upper_range; j++)
        {
          ll_current(j,0) = ll_t_c + gamma*(ll_current(j,0) - ll_t_c);
          if(control.pr_average) ll_current(j,0) = (ll_current(j,0) + ll_pr_c)/((double)iteration);
        }
      } else {
        ll_t = ll_t + gamma*(ll_current.col(0).segment(lower_range, re.mcmc_block_size).mean() - ll_t);
        if(control.pr_average) ll_pr += ll_t;
      }
    }
    if(control.pr_average){
      ll = ll_pr / (double)iteration;
    } else {
      ll = ll_t;
    }
  } else {
    ll = log_likelihood();
  }
  return -1*ll;
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::log_likelihood_beta_with_gradient(const VectorXd& beta, VectorXd& g)
{
  model.linear_predictor.update_parameters(beta.array());
  g.setZero();
  double ll = log_likelihood();
  fn_counter.first += re.scaled_u_.cols();
  if(control.saem)
  {
    throw std::runtime_error("L-BFGS-B not currently available with SAEM");
  } else {
    for(int i = 0; i < re.u_.cols(); i++) g += matrix.log_gradient(re.u_.col(i),true);
    g.array() *= -1.0 / (double) re.u_.cols();
    ll = log_likelihood();
  }
  return -1*ll;
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::log_likelihood_theta(const dblvec& theta)
{
  model.covariance.update_parameters(theta);
  fn_counter.second += re.scaled_u_.cols();
// #pragma omp parallel
  for(int i = 0; i < re.scaled_u_.cols(); i++)
  {
    ll_current(i,1) = model.covariance.log_likelihood(re.scaled_u_.col(i));
  }
  
  if(control.reml){
    // REML correction to the log-likelihood
    MatrixXd D = model.covariance.D().llt().solve(MatrixXd::Identity(Q(),Q()));
    MatrixXd CZZD = CZZ + D;
    CZZD = CZZD.llt().solve(MatrixXd::Identity(Q(),Q()));
    double trCZZ = 0;
    for(int i = 0; i < Q(); i++){
      for(int j = 0; j < Q(); j++){
        trCZZ += D(i,j)*CZZD(j,i);
      }
    }
    ll_current.col(1).array() -= 0.5*trCZZ;
  }
  
  double ll = 0;
  if(control.saem)
  {
    int     iteration = std::max((int)re.zu_.cols() / re.mcmc_block_size, 1);
    double  gamma = pow(1.0/iteration,control.alpha);
    double  ll_t = 0;
    double  ll_pr = 0;
    for(int i = 0; i < iteration; i++){
      int lower_range = i * re.mcmc_block_size;
      int upper_range = (i + 1) * re.mcmc_block_size;
      if(i == (iteration - 1) && iteration > 1){
        double ll_t_c = ll_t;
        double ll_pr_c = ll_pr;
        ll_t = ll_t + gamma*(ll_current.col(1).segment(lower_range, re.mcmc_block_size).mean() - ll_t);
        if(control.pr_average) ll_pr += ll_t;
        for(int j = lower_range; j < upper_range; j++)
        {
          ll_current(j,1) = ll_t_c + gamma*(ll_current(j,1) - ll_t_c);
          if(control.pr_average) ll_current(j,1) = (ll_current(j,1) + ll_pr_c)/((double)iteration);
        }
      } else {
        ll_t = ll_t + gamma*(ll_current.col(1).segment(lower_range, re.mcmc_block_size).mean() - ll_t);
        if(control.pr_average) ll_pr += ll_t;
      }
    }
    if(control.pr_average){
      ll = ll_pr / (double)iteration;
    } else {
      ll = ll_t;
    }
  } else {
    ll = ll_current.col(1).mean();
  }
  
  return -1*ll;
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::log_likelihood_theta_with_gradient(const VectorXd& theta, VectorXd& g)
{
  model.covariance.update_parameters(theta.array());
  fn_counter.second += re.scaled_u_.cols();
  double ll = 0;
  if(control.saem)
  {
    throw std::runtime_error("L-BFGS-B not currently available with SAEM");
  } else {
    g = model.covariance.log_gradient(re.scaled_u_, ll);
  }
  
  if(control.reml){
    // REML correction to the log-likelihood
    std::vector<MatrixXd> derivs;
    model.covariance.derivatives(derivs,1);
    int npars = derivs.size()-1;
    MatrixXd D = model.covariance.D().llt().solve(MatrixXd::Identity(Q(),Q()));
    MatrixXd CZZD = CZZ + D;
    CZZD = CZZD.llt().solve(MatrixXd::Identity(Q(),Q()));
    double trCZZD, trCZZ;
    trCZZ = 0;
    for(int i = 0; i < Q(); i++){
      for(int j = 0; j < Q(); j++){
        trCZZ += D(i,j)*CZZD(j,i);
      }
    }
    ll -= 0.5*trCZZ;
    for(int k = 0; k < npars; k++){
      trCZZD = 0;
      MatrixXd Dderiv = D * derivs[k+1] * D;
      for(int i = 0; i < Q(); i++){
        for(int j = 0; j < Q(); j++){
          trCZZD += D(i,j)*CZZD(j,i);
        }
      }
      g(k) += -0.5*trCZZ;
    }
  }
  
  return -1*ll;
}

template<>
inline double glmmr::ModelOptim<bits_nngp>::log_likelihood_theta_with_gradient(const VectorXd& theta, VectorXd& g)
{
  if(control.reml) throw std::runtime_error("REML not currently available with gradient based NNGP optimisation");
  model.covariance.update_parameters_d(theta.array());
  fn_counter.second += re.scaled_u_.cols();
  double ll = 0;
  if(control.saem)
  {
    throw std::runtime_error("L-BFGS-B not currently available with SAEM");
  } else {
    g = model.covariance.log_gradient(re.scaled_u_, ll);
  }
  return -1*ll;
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::log_likelihood_laplace_beta_u_with_gradient(const VectorXd& x, VectorXd& g)
{
  MatrixXd v(Q(),1);
  v.col(0) = x.tail(Q());
  model.linear_predictor.update_parameters(x.head(P()).array());
  update_u(v);
  double logl = v.col(0).transpose()*v.col(0);
  double ll = log_likelihood();
  matrix.W.update();
  MatrixXd LZWZL = model.covariance.LZWZL(matrix.W.W());
  double LZWdet = glmmr::maths::logdet(LZWZL);
  g.head(P()) = matrix.log_gradient(v.col(0),true);
  g.tail(Q()) = matrix.log_gradient(v.col(0),false);
  g.array() *= -1.0;
  return -1.0*(ll - 0.5*logl - 0.5*LZWdet);
}

template<>
inline double glmmr::ModelOptim<bits_hsgp>::log_likelihood_theta_with_gradient(const VectorXd& theta, VectorXd& g)
{
// THIS NEEDS REWORKING - L-BFGS-B THETA IS NOT INCLUDED IN v0.6.1
throw std::runtime_error("No L-BFGS-B with HSGP");
//   model.covariance.update_parameters(theta.array());
//   double ll = log_likelihood();
//   ArrayXd xb = model.xb();
//   MatrixXd grad(2,re.u_.cols());
//   MatrixXd ZLd0 = model.covariance.ZL_deriv(0);
//   MatrixXd ZLd1 = model.covariance.ZL_deriv(1);
//   int niter = re.u_.cols();
// #pragma omp parallel for if(niter > 30)
//   for(int i = 0; i < niter; i++)
//   { 
//     ArrayXd mu = xb + re.zu_.col(i).array();
//     mu = mu.exp();    
//     grad(0,i) = ( (model.data.y-mu.matrix()) * (re.u_.col(i).transpose()) * ZLd0.transpose()).trace();
//     grad(1,i) = ( (model.data.y-mu.matrix()) * (re.u_.col(i).transpose()) * ZLd1.transpose()).trace();
//   }  
//   g = grad.rowwise().mean();
//   g.array() *= -1.0;
//   return -1.0*ll;
}

template<>
inline double glmmr::ModelOptim<bits_hsgp>::log_likelihood_theta(const dblvec& theta){
  if(control.reml) throw std::runtime_error("REML not currently available with HSGP");
  model.covariance.update_parameters(theta);
  re.zu_ = model.covariance.ZLu(re.u_);
  double ll = log_likelihood(false);
  fn_counter.first += re.scaled_u_.cols();
  if(control.saem)
  {
    int     iteration = std::max((int)re.zu_.cols() / re.mcmc_block_size, 1);
    double  gamma = pow(1.0/iteration,control.alpha);
    double  ll_t = 0;
    double  ll_pr = 0;
    for(int i = 0; i < iteration; i++){
      int lower_range = i * re.mcmc_block_size;
      int upper_range = (i + 1) * re.mcmc_block_size;
      if(i == (iteration - 1) && iteration > 1){
        double ll_t_c = ll_t;
        double ll_pr_c = ll_pr;
        ll_t = ll_t + gamma*(ll_current.col(1).segment(lower_range, re.mcmc_block_size).mean() - ll_t);
        if(control.pr_average) ll_pr += ll_t;
        for(int j = lower_range; j < upper_range; j++)
        {
          ll_current(j,1) = ll_t_c + gamma*(ll_current(j,1) - ll_t_c);
          if(control.pr_average) ll_current(j,1) = (ll_current(j,1) + ll_pr_c)/((double)iteration);
        }
      } else {
        ll_t = ll_t + gamma*(ll_current.col(1).segment(lower_range, re.mcmc_block_size).mean() - ll_t);
        if(control.pr_average) ll_pr += ll_t;
      }
    }
    if(control.pr_average){
      ll = ll_pr / (double)iteration;
    } else {
      ll = ll_t;
    }
  } else {
    ll = log_likelihood(false);
  }
  return -1*ll;
}

template<typename modeltype>
inline glmmr::ModelOptim<modeltype>::ModelOptim(modeltype& model_, 
                                                glmmr::ModelMatrix<modeltype>& matrix_,
                                                glmmr::RandomEffects<modeltype>& re_) : model(model_), matrix(matrix_), re(re_), ll_current(ArrayXXd::Zero(re_.mcmc_block_size,2)) {}; //ll_previous(ArrayXXd::Zero(re_.mcmc_block_size,2)), 

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::set_bobyqa_control(int npt_, double rhobeg_, double rhoend_){
  control.npt = npt_;
  control.rhobeg = rhobeg_;
  control.rhoend = rhoend_;
}


template<typename modeltype>
inline int glmmr::ModelOptim<modeltype>::P() const {
  return model.linear_predictor.P();
}


template<typename modeltype>
inline int glmmr::ModelOptim<modeltype>::Q() const {
  return model.covariance.Q();
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::update_beta(const dblvec &beta){
  bool update = true;
  if(beta_bounded)
  {
    for(int i = 0; i < beta.size(); i++)
    {
      if(beta[i] < lower_bound[i] || beta[i] > upper_bound[i]) 
      {
        update = false;
        throw std::runtime_error("beta out of bounds");
      }
    }
  }
  if(update) model.linear_predictor.update_parameters(beta);
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::update_beta(const VectorXd &beta){
  bool update = true;
  if(beta_bounded)
  {
    for(int i = 0; i < beta.size(); i++)
    {
      if(beta(i) < lower_bound[i] || beta(i) > upper_bound[i]) 
      {
        update = false;
        throw std::runtime_error("beta out of bounds");
      }
    }
  }
  if(update)model.linear_predictor.update_parameters(beta.array());
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::update_theta(const dblvec &theta){
  model.covariance.update_parameters(theta);
  re.zu_ = model.covariance.ZLu(re.u_);
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::update_theta(const VectorXd &theta){
  model.covariance.update_parameters(theta.array());
  re.zu_ = model.covariance.ZLu(re.u_);
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::update_u(const MatrixXd &u_, bool append){
  bool action_append = append;
  // if HSGP then check and update the size of u is m has changed
  if constexpr (std::is_same_v<modeltype,bits_hsgp>){
    if(model.covariance.Q() != re.u_.rows()){
      re.u_.resize(model.covariance.Q(),1);
      re.u_.setZero();
    }
  }
  
  int newcolsize = u_.cols();
  int currcolsize = re.u_.cols();
  // check if the existing samples are a single column of zeros - if so remove them
  if(append && re.u_.cols() == 1 && re.u_.col(0).isZero()) action_append = false;
  // update stored ll values 
  // if(ll_previous.rows() != ll_current.rows()) ll_previous.resize(ll_current.rows(),NoChange);
  // ll_previous = ll_current;
  
  if(action_append){
    re.u_.conservativeResize(NoChange,currcolsize + newcolsize);
    re.zu_.conservativeResize(NoChange,currcolsize + newcolsize);
    re.u_.rightCols(newcolsize) = u_;
    ll_current.resize(currcolsize + newcolsize,NoChange);
  } else {
    if(u_.cols()!=re.u_.cols()){
      re.u_.resize(NoChange,newcolsize);
      re.zu_.resize(NoChange,newcolsize);
    }
    re.u_ = u_;
    if(re.u_.cols() != ll_current.rows()) ll_current.resize(newcolsize,NoChange);
  }
  re.zu_ = model.covariance.ZLu(re.u_);
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::update_u(const MatrixXd &u_){
  int newcolsize = u_.cols();  
  if(u_.cols()!=re.u_.cols()){
      re.u_.resize(NoChange,newcolsize);
      re.zu_.resize(NoChange,newcolsize);
    }
  re.u_ = u_;
  if(newcolsize != ll_current.rows()) ll_current.resize(newcolsize,NoChange);
  re.zu_ = model.covariance.ZLu(re.u_);
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::log_likelihood(bool beta) {
  ArrayXd xb(model.xb());
  int llcol = beta ? 0 : 1;
  ll_current.col(llcol).setZero();
  
  if(model.weighted){
    if(model.family.family==Fam::gaussian){
#pragma omp parallel for
      for(int j= 0; j< re.zu_.cols() ; j++){
        for(int i = 0; i<model.n(); i++){
          ll_current(j,llcol ) += glmmr::maths::log_likelihood(model.data.y(i),xb(i) + re.zu_(i,j),
                                             model.data.variance(i)/model.data.weights(i),
                                             model.family);
        }
      }
    } else {
#pragma omp parallel for
      for(int j=0; j< re.zu_.cols() ; j++){
        for(int i = 0; i<model.n(); i++){
          ll_current(j,llcol) += model.data.weights(i)*glmmr::maths::log_likelihood(model.data.y(i),xb(i) + re.zu_(i,j),
                                   model.data.variance(i),model.family);
        }
      }
      ll_current.col(llcol) *= model.data.weights.sum()/model.n();
    }
  } else {
#pragma omp parallel for
    for(int j= 0; j< re.zu_.cols() ; j++){
      for(int i = 0; i<model.n(); i++){
        ll_current(j,llcol) += glmmr::maths::log_likelihood(model.data.y(i),xb(i) + re.zu_(i,j),
                                           model.data.variance(i),model.family);
      }
    }
  }
  return ll_current.col(llcol).mean();
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::log_likelihood(){
  double ll = log_likelihood(true);
  return ll;
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::full_log_likelihood(){
  double ll = log_likelihood();
  double logl = 0;
  MatrixXd Lu = model.covariance.Lu(re.u(false));
#pragma omp parallel for reduction (+:logl)
  for(int i = 0; i < Lu.cols(); i++){
    logl += model.covariance.log_likelihood(Lu.col(i));
  }
  logl *= 1/Lu.cols();
  return ll+logl;
}

template<typename modeltype>
inline dblvec glmmr::ModelOptim<modeltype>::get_start_values(bool beta, bool theta, bool var)
{
  dblvec start;
  if(beta){
    for(const auto& i: model.linear_predictor.parameters)start.push_back(i);
    if(theta)for(const auto& j: model.covariance.parameters_)start.push_back(j);
  } else {
    start = model.covariance.parameters_;
  }
  if(var && (model.family.family==Fam::gaussian||model.family.family==Fam::gamma||model.family.family==Fam::beta)){
    start.push_back(model.data.var_par);
  }
  return start;
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::set_bound(const dblvec& bound, bool lower)
{
  if(static_cast<int>(bound.size())!=P())throw std::runtime_error("Bound not equal to number of parameters");
  if(lower){
    if(lower_bound.size() != bound.size())lower_bound.resize(P());
    lower_bound = bound; 
  } else {
    if(upper_bound.size() != bound.size())upper_bound.resize(P());
    upper_bound = bound;
  }
  beta_bounded = true;
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::set_theta_bound(const dblvec& bound, bool lower)
{
  if(lower){
    lower_bound_theta = bound; 
  } else {
    upper_bound_theta = bound;
  }
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::use_reml(bool reml)
{
  control.reml = reml;
}


template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::generate_czz()
{
  CZZ.resize(Q(),Q());
  matrix.W.update();
  VectorXd w = matrix.W.W();
  MatrixXd Winv = w.array().inverse().matrix().asDiagonal().toDenseMatrix();
  MatrixXd X = model.linear_predictor.X();
  MatrixXd XtWX = X.transpose() * Winv * X;
  XtWX = XtWX.llt().solve(MatrixXd::Identity(P(),P()));
  MatrixXd S = Winv - Winv * (X * (XtWX * X.transpose())) * Winv;
  CZZ = model.covariance.Z().transpose() * S * model.covariance.Z();
}

template<typename modeltype>
inline dblvec glmmr::ModelOptim<modeltype>::get_lower_values(bool beta, bool theta, bool var, bool u)
{
#ifndef R_BUILD
  double R_NegInf = -1.0 * std::numeric_limits<double>::infinity();
#endif
  dblvec lower;
  if(beta){
    if(lower_bound.size()==0)
    {
      for(int i = 0; i< P(); i++)lower.push_back(R_NegInf);
    } else {
      lower = lower_bound;
    }
  } 
  if(theta)
  {
    if(lower_bound_theta.size()==0)
    {
      for(int i=0; i< model.covariance.npar(); i++)lower.push_back(1e-6);
    } else {
      for(const auto& par: lower_bound_theta)lower.push_back(par);
    }
  }
  if(var && (model.family.family==Fam::gaussian||model.family.family==Fam::gamma||model.family.family==Fam::beta))
  {
    lower.push_back(0.0);
  }
  if(u)
  {
    for(int i = 0; i< Q(); i++) lower.push_back(R_NegInf);
  }
  return lower;
}

template<typename modeltype>
inline dblvec glmmr::ModelOptim<modeltype>::get_upper_values(bool beta, bool theta, bool var, bool u){
  dblvec upper;
#ifndef R_BUILD
  double R_PosInf = std::numeric_limits<double>::infinity();
#endif
  if(beta)
  {
    if(upper_bound.size()==0){
      for(int i = 0; i< P(); i++)upper.push_back(R_PosInf);
    } else {
      upper = upper_bound;
    }
  } 
  if(theta)
  {
    if(upper_bound_theta.size()==0){
      for(int i=0; i< model.covariance.npar(); i++)upper.push_back(R_PosInf);
    } else {
      for(const auto& par: upper_bound_theta)upper.push_back(par);
    }
  }
  if(var && (model.family.family==Fam::gaussian||model.family.family==Fam::gamma||model.family.family==Fam::beta))
  {
    upper.push_back(R_PosInf);
  }
  if(u)
  {
    for(int i = 0; i < Q(); i++) upper.push_back(R_PosInf);
  }
  return upper;
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::nr_beta(){
  // save the old likelihood values
  previous_ll_values.first = current_ll_values.first;
  previous_ll_var.first = current_ll_var.first;
  // if(ll_previous.rows() != ll_current.rows()) ll_previous.resize(ll_current.rows(),NoChange);
  // ll_previous.col(0) = ll_current.col(0);
  
  int niter = re.u(false).cols();
  MatrixXd zd = matrix.linpred();
  ArrayXd sigmas(niter);
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
  Rcpp::Rcout << "\nNR Beta: XtWX";
#endif
  MatrixXd XtXW = MatrixXd::Zero(P()*niter,P());
  MatrixXd Wu = MatrixXd::Zero(model.n(),niter);
  MatrixXd X = model.linear_predictor.X();
  ArrayXd nvar_par(model.n());
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
  
#pragma omp parallel for
  for(int i = 0; i < niter; ++i){
    VectorXd w = glmmr::maths::dhdmu(zd.col(i),model.family);
    w = ((w.array() *nvar_par).inverse() * model.data.weights).matrix();
    ArrayXd resid(model.n());
    matrix.gradient_eta(re.u_.col(i),resid);
    XtXW.block(P()*i, 0, P(), P()) = X.transpose() * w.asDiagonal() * X;
    w = w.cwiseProduct(resid.matrix());
    Wu.col(i) = w;
  }
  XtXW *= (double)1.0/niter;
  MatrixXd XtWXm = XtXW.block(0,0,P(),P());
  for(int i = 1; i<niter; i++) XtWXm += XtXW.block(P()*i,0,P(),P());
  XtWXm = XtWXm.llt().solve(MatrixXd::Identity(P(),P()));
  if(model.linear_predictor.calc.any_nonlinear)
  {
    MatrixXd A = matrix.hessian_nonlinear_correction();
    glmmr::Eigen_ext::near_semi_pd(A);
    XtWXm += A;
  }
  VectorXd Wum = Wu.rowwise().mean();
  VectorXd bincr = XtWXm * X.transpose() * Wum;
  update_beta(model.linear_predictor.parameter_vector() + bincr);
  calculate_var_par();
  
  // repopulate loglikelihood history - assumes NR can 
  // only be used with MCEM for the theta parameters - TO DO: allow NR for beta and SAEM for theta
  current_ll_values.first = log_likelihood();
  current_ll_var.first = (ll_current.col(0) - ll_current.col(0).mean()).square().sum() / (ll_current.col(0).size() - 1);
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::laplace_nr_beta_u(){
  matrix.W.update();
  MatrixXd infomat = matrix.template observed_information_matrix<IM::EIM>();
  infomat = infomat.llt().solve(MatrixXd::Identity(P()+Q(),P()+Q()));
  ArrayXd resid(model.n());
  matrix.gradient_eta(re.u_.col(0),resid);
  VectorXd w = matrix.W.W();
  w = w.cwiseProduct(resid.matrix());
  VectorXd params(P()+Q());
  params.head(P()) = model.linear_predictor.parameter_vector();
  params.tail(Q()) = re.u_.col(0);
  VectorXd pderiv(P()+Q());
  pderiv.head(P()) = (model.linear_predictor.X()).transpose() * w;
  pderiv.tail(Q()) = matrix.log_gradient(re.u_.col(0));
  params += infomat*pderiv;
  update_beta(params.head(P()));
  update_u(params.tail(Q()));
  calculate_var_par();
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::update_var_par(const double& v){
  model.data.var_par = v;
  model.data.variance.setConstant(v);
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::update_var_par(const ArrayXd& v){
  model.data.variance = v;
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::calculate_var_par(){
  if(model.family.family==Fam::gaussian || model.family.family==Fam::quantile_scaled){
    // revise this for beta and Gamma re residuals
    if(control.reml){
      // reml update of marginal variance parameter
      MatrixXd SigmaInv = matrix.Sigma(true);
      MatrixXd X = model.linear_predictor.X();
      MatrixXd Z = model.covariance.Z();
      
      VectorXd xb = model.linear_predictor.xb();
      ArrayXd resid = (model.data.y.array() - xb.array());
      resid *= model.data.weights.sqrt();
      
      MatrixXd M = matrix.template observed_information_matrix<IM::EIM2>();
      M = M.llt().solve(MatrixXd::Identity(P()+Q(),P()+Q()));
      // MatrixXd SX = SigmaInv * X;
      // MatrixXd PP = SigmaInv - SX * M.topLeftCorner(P(),P()) * SX.transpose();
      // VectorXd e = model.data.var_par * model.data.weights.inverse().matrix().asDiagonal() * PP * model.data.y;
      VectorXd e = model.data.var_par * SigmaInv * resid.matrix();
      M.topLeftCorner(P(),P()) = X * M.topLeftCorner(P(),P()) * X.transpose();
      M.bottomRightCorner(Q(),Q()) = Z * M.bottomRightCorner(Q(),Q()) * Z.transpose();
      double new_var_par;
      if((model.data.weights != 1).any()){
        new_var_par = (1.0/model.n()) * (e.transpose() * e + (model.data.weights.inverse().matrix().asDiagonal()*M).trace());
      } else {
        new_var_par = (1.0/model.n()) * (e.transpose()* e + M.trace());
      }
      update_var_par(new_var_par);
    } else {
      int niter = re.u(false).cols();
      ArrayXd sigmas(niter);
      sigmas.setZero();
      MatrixXd zd = matrix.linpred();
#pragma omp parallel for
      for(int i = 0; i < niter; ++i){
        VectorXd zdu = glmmr::maths::mod_inv_func(zd.col(i), model.family.link);
        ArrayXd resid = (model.data.y - zdu);
        resid *= model.data.weights.sqrt();
        if(model.family.family==Fam::gaussian){
          sigmas(i) = (resid - resid.mean()).square().sum()/(resid.size()-1.0);
        } else {
          for(int j = 0; j < resid.size(); j++){
            if(resid(j) < 0){
              sigmas(i) += resid(j)*(model.family.quantile - 1.0);
            } else {
              sigmas(i) += resid(j)*model.family.quantile;
            }
          }
          sigmas(i) *= 1.0/resid.size();
        }
      }
      update_var_par(sigmas.mean());
    }
  }
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::aic(){
  MatrixXd Lu = re.u();
  int dof = P() + model.covariance.npar();
  double logl = 0;
#pragma omp parallel for reduction (+:logl)
  for(int i = 0; i < Lu.cols(); i++){
    logl += model.covariance.log_likelihood(Lu.col(i));
  }
  double ll = log_likelihood();
  
  return (-2*( ll + logl ) + 2*dof); 
}

template<typename modeltype>
inline ArrayXd glmmr::ModelOptim<modeltype>::optimum_weights(double N, 
                                                             VectorXd C,
                                                             double tol,
                                                             int max_iter){
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
  if(C.size()!=P())throw std::runtime_error("C is wrong size");
#endif 
  
  VectorXd Cvec(C);
  ArrayXd weights = ArrayXd::Constant(model.n(),1.0*model.n());
  VectorXd holder(model.n());
  weights = weights.inverse();
  ArrayXd weightsnew(weights);
  ArrayXd w = (matrix.W.W()).array().inverse();
  std::vector<MatrixXd> ZDZ;
  std::vector<MatrixXd> Sigmas;
  std::vector<MatrixXd> Xs;
  std::vector<glmmr::SigmaBlock> SB(matrix.get_sigma_blocks());
#ifdef R_BUILD
  Rcpp::Rcout << "\n### Preparing data ###";
  Rcpp::Rcout << "\nThere are " << SB.size() << " independent blocks and " << model.n() << " cells.";
#endif
  int maxprint = model.n() < 10 ? model.n() : 10;
  for(auto& sb: SB){
    sparse ZLs = submat_sparse(model.covariance.ZL_sparse(),sb.RowIndexes);
    MatrixXd ZL = sparse_to_dense(ZLs,false);
    MatrixXd S = ZL * ZL.transpose();
    ZDZ.push_back(S);
    Sigmas.push_back(S);
    ArrayXi rows = Map<ArrayXi,Unaligned>(sb.RowIndexes.data(),sb.RowIndexes.size());
    MatrixXd X = glmmr::Eigen_ext::submat(model.linear_predictor.X(),rows,ArrayXi::LinSpaced(P(),0,P()-1));
    Xs.push_back(X);
  }
  
  double diff = 1;
  int block_size;
  MatrixXd M(P(),P());
  int iter = 0;
#ifdef R_BUILD
  Rcpp::Rcout << "\n### Starting optimisation ###";
#endif
  while(diff > tol && iter < max_iter){
    iter++;
#ifdef R_BUILD
    Rcpp::Rcout << "\nIteration " << iter << "\n------------\nweights: [" << weights.segment(0,maxprint).transpose() << " ...]";
#endif
    //add check to remove weights that are below a certain threshold
    if((weights < 1e-8).any()){
      for(int i = 0 ; i < SB.size(); i++){
        auto it = SB[i].RowIndexes.begin();
        while(it != SB[i].RowIndexes.end()){
          if(weights(*it) < 1e-8){
            weights(*it) = 0;
            int idx = it - SB[i].RowIndexes.begin();
            glmmr::Eigen_ext::removeRow(Xs[i],idx);
            glmmr::Eigen_ext::removeRow(ZDZ[i],idx);
            glmmr::Eigen_ext::removeColumn(ZDZ[i],idx);
            Sigmas[i].conservativeResize(ZDZ[i].rows(),ZDZ[i].cols());
            it = SB[i].RowIndexes.erase(it);
#ifdef R_BUILD
            Rcpp::Rcout << "\n Removing point " << idx << " in block " << i;
#endif
          } else {
            it++;
          }
        }
      }
    }
    
    M.setZero();
    for(int i = 0 ; i < static_cast<int>(SB.size()); i++){
      Sigmas[i] = ZDZ[i];
      for(int j = 0; j < Sigmas[i].rows(); j++){
        // sigma_sq
        Sigmas[i](j,j) += w(SB[i].RowIndexes[j])/(N*weights(SB[i].RowIndexes[j]));
      }
      Sigmas[i] = Sigmas[i].llt().solve(MatrixXd::Identity(Sigmas[i].rows(),Sigmas[i].cols()));
      M += Xs[i].transpose() * Sigmas[i] * Xs[i];
    }
    
    //check if positive definite, if not remove the offending column(s)
    bool isspd = glmmr::Eigen_ext::issympd(M);
    if(isspd){
#ifdef R_BUILD
      Rcpp::Rcout << "\n Information matrix not postive definite: ";
#endif
      ArrayXd M_row_sums = M.rowwise().sum();
      int fake_it = 0;
      int countZero = 0;
      for(int j = 0; j < M_row_sums.size(); j++){
        if(M_row_sums(j) == 0){
#ifdef R_BUILD
          Rcpp::Rcout << "\n   Removing column " << fake_it;
#endif
          for(int k = 0; k < Xs.size(); k++){
            glmmr::Eigen_ext::removeColumn(Xs[k],fake_it);
          }
          glmmr::Eigen_ext::removeElement(Cvec,fake_it);
          countZero++;
        } else {
          fake_it++;
        }
      }
      M.conservativeResize(M.rows()-countZero,M.cols()-countZero);
      M.setZero();
      for(int k = 0; k < static_cast<int>(SB.size()); k++){
        M += Xs[k].transpose() * Sigmas[k] * Xs[k];
      }
    }
    M = M.llt().solve(MatrixXd::Identity(M.rows(),M.cols()));
    VectorXd Mc = M*Cvec;
    weightsnew.setZero();
    for(int i = 0 ; i < SB.size(); i++){
      block_size = SB[i].RowIndexes.size();
      holder.segment(0,block_size) = Sigmas[i] * Xs[i] * Mc;
      for(int j = 0; j < block_size; j++){
        weightsnew(SB[i].RowIndexes[j]) = holder(j);
      }
    }
    weightsnew = weightsnew.abs();
    weightsnew *= 1/weightsnew.sum();
    diff = ((weights-weightsnew).abs()).maxCoeff();
    weights = weightsnew;
#ifdef R_BUILD
    Rcpp::Rcout << "\n(Max. diff: " << diff << ")\n";
#endif
  }
#ifdef R_BUILD
  if(iter<max_iter){
    Rcpp::Rcout << "\n### CONVERGED Final weights: [" << weights.segment(0,maxprint).transpose() << "...]";
  } else {
    Rcpp::Rcout << "\n### NOT CONVERGED Reached maximum iterations";
  }
#endif
  return weights;
}