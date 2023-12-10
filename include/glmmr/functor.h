// Copyright (c) 2022 Yi Pan <ypan1988@gmail.com>
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  A copy of the GNU General Public License is available at
//  https://www.R-project.org/Licenses/

#ifndef FUNCTOR_H_
#define FUNCTOR_H_

namespace rminqa {

template<typename T>
class Functor {
public:
  struct OptStruct {
    bool has_grad_ = false;
    bool has_hess_ = false;
    T ndeps_;       // tolerances for numerical derivatives
    double fnscale_ = 1.0;  // scaling for objective
    T parscale_;    // scaling for parameters
    int usebounds_ = 0;
    T lower_, upper_;
    bool sann_use_custom_function_ = false;
  } os;

public:
  int feval = 0;

  Functor() {}
  virtual ~Functor() {}
  virtual double operator()(const std::vector<double> &par) = 0;
  
  // do we need virtual gradient and hessian functions? will they ever be overridden?
  void Gradient(const T &par, T &grad){
    if (os.parscale_.empty()) {
      os.parscale_ = T(par.size());
      for(int i = 0; i< par.size(); i++)os.parscale_[i] = 1.0;
    } 
    if (os.ndeps_.empty()){
      os.ndeps_ = T(par.size());
      for(int i = 0; i< par.size(); i++)os.ndeps_[i] = 1e-6;
    }
    
    grad = T(par.size(),0.0);;// arma::zeros<arma::vec>(par.size());
    T x(par.size());
    
    for(int i = 0; i< (int)par.size(); i++)x[i] = par[i]*os.parscale_[i];
    
    if (os.usebounds_ == 0) {
      for (int i = 0; i != par.size(); ++i) {
        double eps = os.ndeps_[i];
        
        x[i] = (par[i] + eps) * os.parscale_[i];
        double val1 = operator()(x) / os.fnscale_;
        
        x[i] = (par[i] - eps) * os.parscale_[i];
        double val2 = operator()(x) / os.fnscale_;
        
        grad[i] = (val1 - val2) / (2 * eps);
        
        x[i] = par[i] * os.parscale_[i];
      }
    } else {  // use bounds
      for (int i = 0; i != par.size(); ++i) {
        double epsused = os.ndeps_[i];
        double eps = os.ndeps_[i];
        
        double tmp = par[i] + eps;
        if (tmp > os.upper_[i]) {
          tmp = os.upper_[i];
          epsused = tmp - par[i];
        }
        
        x[i] = tmp * os.parscale_[i];
        double val1 = operator()(x) / os.fnscale_;
        
        tmp = par[i] - eps;
        if (tmp < os.lower_[i]) {
          tmp = os.lower_[i];
          eps = par[i] - tmp;
        }
        
        x[i] = tmp * os.parscale_[i];
        double val2 = operator()(x) / os.fnscale_;
        
        grad[i] = (val1 - val2) / (epsused + eps);
        
        x[i] = par[i] * os.parscale_[i];
      }
    }
  }
  
  // currently hess needs to be a type which the elements can be access using [n] in column-major
  // order. For Eigen library this means either passing a vector and then resizing it, or passing
  // it an std::vector with the data 
  template<typename T2>
  void Hessian(const T &par, T2 &hess){
    int n = par.size();
    if (os.parscale_.empty()) {
      os.parscale_ = T(par.size());
      for(int i = 0; i< par.size(); i++)os.parscale_[i] = 1.0;
    } 
    if (os.ndeps_.empty()){
      os.ndeps_ = T(par.size());
      for(int i = 0; i< par.size(); i++)os.ndeps_[i] = 1e-3;
    }
    
    //hess = std::vector<double>(n*n,0.0);// arma::zeros<arma::mat>(par.size(), par.size());
    for(int i = 0; i< hess.size(); i++)hess[i] = 0.0;
    std::vector<double> dpar(n);// = par / os.parscale_;
    for(int i = 0; i < n; i++)dpar[i] = par[i]/os.parscale_[i];
    std::vector<double> df1(n,0.0);
    std::vector<double> df2(n,0.0);
    
    for (int i = 0; i != n; ++i) {
      double eps = os.ndeps_[i] / os.parscale_[i];
      dpar[i] += eps;
      Gradient(dpar, df1);
      dpar[i] -= 2 * eps;
      Gradient(dpar, df2);
      for (int j = 0; j < n; ++j)
        hess[i*n + j] = os.fnscale_ * (df1[j] - df2[j]) /
          (2 * eps * os.parscale_[i] * os.parscale_[j]); //this was in row-major order but have changed
      dpar[i] = dpar[i] + eps;
    }
    
    // now symmetrize
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        double tmp = 0.5 * (hess[i*n + j] + hess[j*n + i]);
        hess[i*n + j] = tmp;
        hess[j*n + i] = tmp;
      }
    }
    
    //this will need casting back to a matrix or other type
  }
  
};

inline double minqa_objfun(long n, const double *x, void *data) {
  std::vector<double> par;
  par.assign(x,x+n);
  ++(static_cast<Functor<std::vector<double> > *>(data)->feval);
  return static_cast<Functor<std::vector<double> > *>(data)->operator()(par);
}

} // namespace rminqa

#endif
