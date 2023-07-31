#ifndef MATHS_H
#define MATHS_H

#define _USE_MATH_DEFINES

#include <cmath> 
#include <unordered_map>
#include <Eigen/Core>
#include <Eigen/LU>
#include <boost/math/distributions/normal.hpp>
#include "algo.h"
#include "general.h"
#include "family.hpp"

namespace glmmr {

    template<class T>
    T randomGaussian(T generator,
        VectorXd& res)
    {
        for (size_t i = 0; i < res.size(); ++i)
            res(i) = generator();
        // Note the generator is returned back
        return  generator;
    }

namespace maths {

inline double gaussian_cdf(double value)
{
  return 0.5 * erfc(-value * 0.707106781186547524401);
}

inline Eigen::VectorXd gaussian_cdf_vec(const Eigen::VectorXd& v) {
  Eigen::VectorXd res(v.size());
  for (int i = 0; i < v.size(); ++i)
    res(i) = gaussian_cdf(v(i));
  return res;
}

template <typename T>
inline T gaussian_pdf(T x)
{
  static const T inv_sqrt_2pi = 0.3989422804014327;
  return inv_sqrt_2pi * std::exp(-T(0.5) * x * x);
}

inline Eigen::VectorXd gaussian_pdf_vec(const Eigen::VectorXd& v)
{
  Eigen::VectorXd res(v.size());
  for (int i = 0; i < v.size(); ++i)
    res(i) = gaussian_pdf(v(i));
  return res;
}

inline Eigen::VectorXd exp_vec(const Eigen::VectorXd& x,
                               bool logit = false)
{
  Eigen::VectorXd z(x.size());
  for (int i = 0; i < x.size(); i++)
  {
    z(i) = logit ? std::exp(x(i)) / (1 + std::exp(x(i))) : std::exp(x(i));
  }
  return z;
}

inline Eigen::VectorXd mod_inv_func(Eigen::VectorXd mu,
                                    std::string link)
{
  const static std::unordered_map<std::string, int> string_to_case{
    {"logit",1},
    {"log",2},
    {"probit",3},
    {"identity",4},
    {"inverse",5}
  };
  switch (string_to_case.at(link)) {
  case 1:
    mu = exp_vec(mu, true);
    break;
  case 2:
    mu = exp_vec(mu);
    break;
  case 3:
    mu = gaussian_cdf_vec(mu);
    break;
  case 4:
    break;
  case 5:
    mu = mu.array().inverse().matrix();
    break;
  }
  return mu;
}

inline Eigen::VectorXd dhdmu(const Eigen::VectorXd& xb,
                             const glmmr::Family& family) {
  
  Eigen::VectorXd wdiag(xb.size());
  Eigen::ArrayXd p(xb.size());
  
  switch (family.flink) {
  case 1:
    wdiag = exp_vec(-1.0 * xb);
    break;
  case 2:
    wdiag = exp_vec(xb);
    break;
  case 3: case 13:
    p = mod_inv_func(xb, "logit");
    for(int i =0; i< xb.size(); i++){
      wdiag(i) = 1/(p(i)*(1.0 - p(i)));
    }
    break;
  case 4: case 14:
    p = mod_inv_func(xb, "log");
    for(int i =0; i< xb.size(); i++){
      wdiag(i) = (1.0 - p(i))/p(i);
    }
    break;
  case 5: case 15:
    p = mod_inv_func(xb, "identity");
    wdiag = (p * (1 - p)).matrix();
    break;
  case 6: case 16:
    {
      p = mod_inv_func(xb, "probit");
      Eigen::ArrayXd pinv = gaussian_pdf_vec(xb);
      for(int i =0; i< xb.size(); i++){
        wdiag(i) = (p(i) * (1-p(i)))/pinv(i);
      }
      break;
    }
  case 7:
    for(int i =0; i< xb.size(); i++){
      wdiag(i) = 1.0;
    }
    break;
  case 8:
    for(int i =0; i< xb.size(); i++){
      wdiag(i) = 1/exp(xb(i));
    }
    break;
  case 9:
    for(int i =0; i< xb.size(); i++){
      wdiag(i) = 1.0;
    }
    break;
  case 10:
    for(int i =0; i< xb.size(); i++){
      wdiag(i) = 1/(xb(i)*xb(i));
    }
    break;
  case 11:
    for(int i =0; i< xb.size(); i++){
      wdiag(i) = (xb(i)*xb(i));
    }
    break;
  case 12:
    p = mod_inv_func(xb, "logit");
    for(int i =0; i< xb.size(); i++){
      wdiag(i) = 1/(p(i)*(1.0 - p(i)));
    }
    break;
  }
  return wdiag;
}

inline Eigen::VectorXd detadmu(const Eigen::VectorXd& xb,
                               std::string link) {
  
  Eigen::VectorXd wdiag(xb.size());
  Eigen::VectorXd p(xb.size());
  const static std::unordered_map<std::string, int> string_to_case{
    {"log",1},
    {"identity",2},
    {"logit",3},
    {"probit",4},
    {"inverse",5}
  };
  
  switch (string_to_case.at(link)) {
  case 1:
    wdiag = glmmr::maths::exp_vec(-1.0 * xb);
    break;
  case 2:
    for(int i =0; i< xb.size(); i++){
      wdiag(i) = 1.0;
    }
    break;
  case 3:
    p = glmmr::maths::mod_inv_func(xb, "logit");
    for(int i =0; i< xb.size(); i++){
      wdiag(i) = 1/(p(i)*(1.0 - p(i)));
    }
    break;
  case 4:
    {
      Eigen::ArrayXd pinv = gaussian_pdf_vec(xb);
      wdiag = (pinv.inverse()).matrix();
      break;
    }
  case 5:
    for(int i =0; i< xb.size(); i++){
      wdiag(i) = -1.0 * xb(i) * xb(i);
    }
    break;
    
  }
  
  return wdiag;
}

inline Eigen::VectorXd attenuted_xb(const Eigen::VectorXd& xb,
                                    const Eigen::MatrixXd& Z,
                                    const Eigen::MatrixXd& D,
                                    const std::string& link){
  Eigen::ArrayXd xbnew(xb.array());
  int n = xb.size();
  if(link=="log"){
    for(int i=0; i<n; i++){
      xbnew(i) += (Z.row(i)*D*Z.row(i).transpose())(0)/2;
    }
  } else if(link=="probit"){
    Eigen::ArrayXd zprod(n);
    Eigen::MatrixXd Dzz(D.rows(),D.cols());
    Eigen::PartialPivLU<Eigen::MatrixXd> pluDzz;
    for(int i=0; i<n; i++){
      Dzz = D*Z.row(i).transpose()*Z.row(i) + Eigen::MatrixXd::Identity(D.rows(),D.cols());
      pluDzz = Eigen::PartialPivLU<Eigen::MatrixXd>(Dzz);
      zprod(i) = pluDzz.determinant();
    }
    zprod = zprod.inverse().sqrt();
    xbnew *= zprod;
  } else if(link=="logit"){
    double c = 0.5880842;
    Eigen::ArrayXd zprod(n);
    Eigen::MatrixXd Dzz(D.rows(),D.cols());
    Eigen::PartialPivLU<Eigen::MatrixXd> pluDzz;
    for(int i=0; i<n; i++){
      Dzz = c*D*Z.row(i).transpose()*Z.row(i) + Eigen::MatrixXd::Identity(D.rows(),D.cols());
      pluDzz = Eigen::PartialPivLU<Eigen::MatrixXd>(Dzz);
      zprod(i) = pluDzz.determinant();
    }
    zprod = zprod.inverse().sqrt();
    xbnew *= zprod;
  }
  
  return xbnew.matrix();
}

inline Eigen::VectorXd marginal_var(const Eigen::VectorXd& mu,
                                    const std::string& family,
                                    double var_par = 1.0){
  Eigen::ArrayXd wdiag(mu.size());
  const static std::unordered_map<std::string, int> string_to_case{
    {"gaussian",1},
    {"binomial",2},
    {"poisson",3},
    {"Gamma",4},
    {"beta",5}
  };
  
  switch (string_to_case.at(family)) {
  case 1:
    wdiag.setConstant(var_par);
    break;
  case 2:
    wdiag = mu.array()*(1-mu.array());
    break;
  case 3:
    wdiag = mu.array();
    break;
  case 4:
    wdiag = mu.array().square();
    break;
  case 5:
    wdiag = mu.array()*(1-mu.array())/(var_par+1);
    break;
  }
  
  return wdiag.matrix();
}

//ramanujans approximation
inline double log_factorial_approx(double n){
  double ans;
  if(n==0){
    ans = 0;
  } else {
    ans = n*log(n) - n + log(n*(1+4*n*(1+2*n)))/6 + log(3.141593)/2;
  }
  return ans;
}

inline double log_likelihood(double y,
                             double mu,
                             double var_par,
                             int flink) {
  double logl = 0;
  
  switch (flink){
  case 1:
  {
    double lf1 = glmmr::maths::log_factorial_approx(y);
    logl = y * mu - exp(mu) - lf1;
    break;
  }
  case 2:
  {
    double lf1 = log_factorial_approx(y);
    logl = y*log(mu) - mu-lf1;
    break;
  }
  case 3:
    if(y==1){
      logl = log(1/(1+exp(-1.0*mu)));
    } else {
      logl = log(1 - 1/(1+exp(-1.0*mu)));
    }
    break;
  case 4:
    if(y==1){
      logl = mu;
    } else {
      logl = log(1 - exp(mu));
    }
    break;
  case 5:
    if(y==1){
      logl = log(mu);
    } else {
      logl = log(1 - mu);
    }
    break;
  case 6:
  {
      boost::math::normal norm(0,1);
      if (y == 1) {
          logl = (double)cdf(norm,mu);
      }
      else {
          logl = log(1 - (double)cdf(norm, mu));
      }
      break;
  }
      
  case 7:
    logl = -0.5*log(var_par) -0.5*log(2*3.141593) -
      0.5*(y - mu)*(y - mu)/var_par;
    break;
  case 8:
    logl = -0.5*log(var_par) -0.5*log(2*3.141593) -
      0.5*(log(y) - mu)*(log(y) - mu)/var_par;
    break;
  case 9:
      {
        double ymu = var_par*y/exp(mu);
        logl = log(1/(tgamma(var_par)*y)) + var_par*log(ymu) - ymu;
        break;
      }
  case 10:
      {
        double ymu = var_par*y*mu;
        logl = log(1/(tgamma(var_par)*y)) + var_par*log(ymu) - ymu;
        break;
      }
  case 11:
    logl = log(1/(tgamma(var_par)*y)) + var_par*log(var_par*y/mu) - var_par*y/mu;
    break;
  case 12:
    logl = (mu*var_par - 1)*log(y) + ((1-mu)*var_par - 1)*log(1-y) - lgamma(mu*var_par) - lgamma((1-mu)*var_par) + lgamma(var_par);
  case 13:
    {
      double lfk = glmmr::maths::log_factorial_approx(y);
      double lfn = glmmr::maths::log_factorial_approx(var_par);
      double lfnk = glmmr::maths::log_factorial_approx(var_par - y);
      logl = lfn - lfk - lfnk + y*log(1/(1+exp(-1.0*mu))) + (var_par - y)*log(1 - 1/(1+exp(-1.0*mu)));
      break;
    }
  case 14:
    {
      double lfk = glmmr::maths::log_factorial_approx(y);
      double lfn = glmmr::maths::log_factorial_approx(var_par);
      double lfnk = glmmr::maths::log_factorial_approx(var_par - y);
      logl = lfn - lfk - lfnk + y*mu + (var_par - y)*log(1 - exp(mu));
      break;
    }
  case 15:
    {
      double lfk = glmmr::maths::log_factorial_approx(y);
      double lfn = glmmr::maths::log_factorial_approx(var_par);
      double lfnk = glmmr::maths::log_factorial_approx(var_par - y);
      logl = lfn - lfk - lfnk + y*log(mu) + (var_par - y)*log(1 - mu);
      break;
    }
  case 16:
    {
      boost::math::normal norm(0, 1);
      double lfk = glmmr::maths::log_factorial_approx(y);
      double lfn = glmmr::maths::log_factorial_approx(var_par);
      double lfnk = glmmr::maths::log_factorial_approx(var_par - y);
      logl = lfn - lfk - lfnk + y*((double)cdf(norm,mu)) + (var_par - y)*log(1 - (double)cdf(norm,mu));
      break;
    }
  }
  return logl;
}

inline double logdet(const Eigen::MatrixXd& M) {
  double ld = 0;
  Eigen::LLT<Eigen::MatrixXd> chol(M);
  auto& U = chol.matrixL();
  for (unsigned i = 0; i < M.rows(); ++i)
    ld += log(U(i,i));
  ld *= 2;
  return ld;
}

inline MatrixXd sample_MVN(const vector_matrix& mu,
                           int m){
  int n = mu.vec.size();
  MatrixXd L = mu.mat.llt().matrixL();
  boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
      generator(boost::mt19937(time(0)),
          boost::normal_distribution<>());
  VectorXd z(n);
  MatrixXd samps(n,m);
  for(int i = 0; i < m; i++){
      randomGaussian(generator, z);
    samps.col(i) = z;
    samps.col(i) += mu.vec;
  }
  return samps;
}
  
}

namespace tools {
inline std::vector<int> linseq(int start, int end) {
  std::vector<int> idx;
  for (int i = start; i <= end; i++) {
    idx.push_back(i);
  }
  return idx;
}
}
}


#endif