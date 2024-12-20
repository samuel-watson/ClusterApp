#pragma once

#include <boost/math/distributions/normal.hpp>
#include <boost/random.hpp>
#include "general.h"
#include "algo.h"
#include "family.hpp"

// [[Rcpp::depends(RcppEigen)]]

struct VectorMatrix {
public:
  VectorXd vec;
  MatrixXd mat;
  VectorMatrix(int n): vec(n), mat(n,n) {};
  VectorMatrix(const VectorMatrix& x) : vec(x.vec), mat(x.mat) {};
  VectorMatrix& operator=(VectorMatrix x){
    vec = x.vec;
    mat = x.mat;
    return *this;
  };
};

struct MatrixMatrix {
public:
  MatrixXd mat1;
  MatrixXd mat2;
  double a = 0;
  double b = 0;
  MatrixMatrix(int n1, int m1, int n2, int m2): mat1(n1,m1), mat2(n2,m2) {};
  MatrixMatrix(const MatrixMatrix& x) : mat1(x.mat1), mat2(x.mat2) {};
  MatrixMatrix& operator=(MatrixMatrix x){
    mat1 = x.mat1;
    mat2 = x.mat2;
    a = x.a;
    b = x.b;
    return *this;
  };
};

namespace glmmr {

template<class T>
inline T randomGaussian(T generator,
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

inline VectorXd mod_inv_func(const VectorXd& muin,
                                    Link link)
{
  VectorXd mu(muin);
  switch (link) {
  case Link::logit:
    mu = exp_vec(mu, true);
    break;
  case Link::loglink:
    mu = exp_vec(mu);
    break;
  case Link::probit:
    mu = gaussian_cdf_vec(mu);
    break;
  case Link::identity:
    break;
  case Link::inverse:
    mu = mu.array().inverse().matrix();
    break;
  }
  return mu;
}

inline double mod_inv_func(const double& muin,
                             Link link)
{
  double mu(muin);
  switch (link) {
  case Link::logit:
    mu = 1/(1+std::exp(-1.0*mu));
    break;
  case Link::loglink:
    mu = std::exp(mu);
    break;
  case Link::probit:
    mu = gaussian_cdf(mu);
    break;
  case Link::identity:
    break;
  case Link::inverse:
    mu = 1/mu;
    break;
  }
  return mu;
}



inline Eigen::VectorXd dhdmu(const Eigen::VectorXd& xb,
                             const glmmr::Family& family) {
  Eigen::VectorXd wdiag(xb.size());
  Eigen::ArrayXd p(xb.size());
  
  switch(family.family){
    case Fam::poisson:
      {
        switch(family.link){
          case Link::identity:
            wdiag = exp_vec(xb);
            break;
          default:
            wdiag = exp_vec(-1.0 * xb);
            break;
          }
      break;
      }
    case Fam::bernoulli: case Fam::binomial:
      {
        switch(family.link){
          case Link::loglink:
            p = mod_inv_func(xb, family.link);
            for(int i =0; i< xb.size(); i++){
              wdiag(i) = (1.0 - p(i))/p(i);
            }
            break;
          case Link::identity:
            p = mod_inv_func(xb, family.link);
            wdiag = (p * (1 - p)).matrix();
            break;
          case Link::probit:
          {
            p = mod_inv_func(xb, family.link);
            Eigen::ArrayXd pinv = gaussian_pdf_vec(xb);
            for(int i =0; i< xb.size(); i++){
              wdiag(i) = (p(i) * (1-p(i)))/pinv(i);
            }
            break;
          }
          default:
            p = mod_inv_func(xb, family.link);
            for(int i =0; i< xb.size(); i++){
              wdiag(i) = 1/(p(i)*(1.0 - p(i)));
            }
            break;
        }
      break;
      }
    case Fam::gaussian:
      {
        switch(family.link){
          case Link::loglink:
            for(int i =0; i< xb.size(); i++){
              wdiag(i) = 1/exp(xb(i));
            }
            break;
          default:
            //identity
            wdiag.setConstant(1.0);
            break;
        }
      break;
      }
    case Fam::gamma:
      {
        switch(family.link){
          case Link::inverse:
            for(int i =0; i< xb.size(); i++){
              wdiag(i) = 1/(xb(i)*xb(i));
            }
            break;
          case Link::identity:
            for(int i =0; i< xb.size(); i++){
              wdiag(i) = (xb(i)*xb(i));
            }
            break;
          default:
            //log
            wdiag.setConstant(1.0);
            break;
        }
      break;
      }
    case Fam::beta:
      {
        //only logit currently
        p = mod_inv_func(xb, family.link);
        for(int i =0; i< xb.size(); i++){
          wdiag(i) = 1/(p(i)*(1.0 - p(i)));
        }
        break;
      }
    case Fam::quantile: case Fam::quantile_scaled:
    {
      switch(family.link){
        case Link::loglink:
          p = mod_inv_func(xb, family.link);
          p = p.square().inverse();
          break;
        case Link::identity:
          wdiag.setConstant(1.0);
          break;
        case Link::probit:
        {
          p = gaussian_pdf_vec(xb);
          p = p.inverse();
          break;
        }
        case Link::inverse:
        {
          p = xb.array().square().inverse();
          break;
        }
        case Link::logit:
          p = mod_inv_func(xb, family.link);
          p = p*(1-p);
          p = p.square().inverse();
          break;
        }
      
      // p *= (family.quantile * family.quantile)/(1 + pow(family.quantile,4));
      break;
    }
  }
  return wdiag;
}

inline Eigen::VectorXd detadmu(const Eigen::VectorXd& xb,
                               const Link link) {
  Eigen::VectorXd wdiag(xb.size());
  Eigen::VectorXd p(xb.size());
  
  switch (link) {
  case Link::loglink:
    wdiag = glmmr::maths::exp_vec(-1.0 * xb);
    break;
  case Link::identity:
    for(int i =0; i< xb.size(); i++){
      wdiag(i) = 1.0;
    }
    break;
  case Link::logit:
    p = glmmr::maths::mod_inv_func(xb, link);
    for(int i =0; i< xb.size(); i++){
      wdiag(i) = 1/(p(i)*(1.0 - p(i)));
    }
    break;
  case Link::probit:
    {
      Eigen::ArrayXd pinv = gaussian_pdf_vec(xb);
      wdiag = (pinv.inverse()).matrix();
      break;
    }
  case Link::inverse:
    for(int i =0; i< xb.size(); i++){
      wdiag(i) = -1.0 * xb(i) * xb(i);
    }
    break;
    
  }
  return wdiag;
}

// inline double log_likelihood_derivative(const double mu, 
//                                         const double y,
//                                         int order, 
//                                         const Fam family,
//                                         double var_par = 1.0)
// {
//   if(order < 1 || order > 2)throw std::runtime_error("Only orders 1 and 2 supported for derivatives");
//   double result;
//   
//   if(order == 1){
//     switch (family) {
//     case Fam::gaussian:
//       result = (y - mu)/var_par;
//       break;
//     case Fam::bernoulli: 
//       result = y/mu - (1-y)/(1-mu);
//       break;
//     case Fam::binomial:
//       result = y/mu - (var_par-y)/(1-mu);
//       break;
//     case Fam::poisson:
//       result = y/mu - 1;
//       break;
//     case Fam::gamma:
//       result = 0;
//       break;
//     case Fam::beta:
//       result = 0;
//       break;
//     }
//   } else {
//     switch (family) {
//     case Fam::gaussian:
//       result = -1.0/var_par;
//       break;
//     case Fam::bernoulli: 
//       result = -1.0*y/(mu*mu) - (1-y)/((1-mu)*(1-mu));
//       break;
//     case Fam::binomial:
//       result = -1.0*y/(mu*mu) - (var_par-y)/((1-mu)*(1-mu));
//       break;
//     case Fam::poisson:
//       result = -1.0*y/(mu*mu);
//       break;
//     case Fam::gamma:
//       result = 0;
//       break;
//     case Fam::beta:
//       result = 0;
//       break;
//     }
//   }
//   return result;
// }
// 
// inline double link_derivative(const double eta, 
//                               int order, 
//                               const Link link)
// {
//   if(order < 1 || order > 2)throw std::runtime_error("Only orders 1 and 2 supported for derivatives");
//   double result;
//   
//   if(order == 1){
//     switch (link) {
//     case Link::identity:
//       result = 1.0;
//       break;
//     case Link::loglink:
//       result = exp(eta);
//       break;
//     case Link::logit:
//       {
//         double exp_eta = exp(eta);
//         result = exp_eta/((1+exp_eta)*(1+exp_eta));
//         break;
//       }
//     default:
//       result = 0;
//       break;
//     }
//   } else {
//     switch (link) {
//     case Link::identity:
//       result = 0.0;
//       break;
//     case Link::loglink:
//       result = exp(eta);
//       break;
//     case Link::logit:
//     {
//       double exp_eta = exp(eta);
//       result = exp_eta*(1+exp_eta)/((1+exp_eta)*(1+exp_eta)*(1+exp_eta));
//       break;
//     }
//     default:
//       result = 0;
//       break;
//     }
//   }
//   return result;
// }



inline Eigen::VectorXd attenuted_xb(const Eigen::VectorXd& xb,
                                    const Eigen::MatrixXd& Z,
                                    const Eigen::MatrixXd& D,
                                    const Link link){
  Eigen::ArrayXd xbnew(xb.array());
  int n = xb.size();
  if(link==Link::loglink){
    for(int i=0; i<n; i++){
      xbnew(i) += (Z.row(i)*D*Z.row(i).transpose())(0)/2;
    }
  } else if(link==Link::probit){
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
  } else if(link==Link::logit){
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
                                    const Fam family,
                                    double var_par = 1.0){
  Eigen::ArrayXd wdiag(mu.size());
  
  switch (family) {
  case Fam::gaussian: case Fam::quantile: case Fam::quantile_scaled:
    wdiag.setConstant(var_par);
    break;
  case Fam::bernoulli: case Fam::binomial:
    wdiag = mu.array()*(1-mu.array());
    break;
  case Fam::poisson:
    wdiag = mu.array();
    break;
  case Fam::gamma:
    wdiag = mu.array().square();
    break;
  case Fam::beta:
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

inline double log_likelihood(const double y,
                             const double mu,
                             const double var_par,
                             const glmmr::Family& family) {
  double logl = 0;
  switch(family.family){
  case Fam::poisson:
  {
    switch(family.link){
  case Link::identity:
  {
    double lf1 = log_factorial_approx(y);
    logl = y*log(mu) - mu-lf1;
    break;
  }
  default:
    {
      double lf1 = glmmr::maths::log_factorial_approx(y);
      logl = y * mu - exp(mu) - lf1;
      break;
    }
  }
    break;
  }
  case Fam::bernoulli:
  {
    switch(family.link){
  case Link::loglink:
    if(y==1){
      logl = mu;
    } else {
      logl = log(1 - exp(mu));
    }
    break;
  case Link::identity:
    if(y==1){
      logl = log(mu);
    } else {
      logl = log(1 - mu);
    }
    break;
  case Link::probit:
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
  default:
    //logit
    if(y==1){
      logl = log(1/(1+exp(-1.0*mu)));
    } else {
      logl = log(1 - 1/(1+exp(-1.0*mu)));
    }
    break;
  }
    break;
  }
  case Fam::binomial:
  {
    switch(family.link){
  case Link::loglink:
  {
    double lfk = glmmr::maths::log_factorial_approx(y);
    double lfn = glmmr::maths::log_factorial_approx(var_par);
    double lfnk = glmmr::maths::log_factorial_approx(var_par - y);
    logl = lfn - lfk - lfnk + y*mu + (var_par - y)*log(1 - exp(mu));
    break;
  }
  case Link::identity:
  {
    double lfk = glmmr::maths::log_factorial_approx(y);
    double lfn = glmmr::maths::log_factorial_approx(var_par);
    double lfnk = glmmr::maths::log_factorial_approx(var_par - y);
    logl = lfn - lfk - lfnk + y*log(mu) + (var_par - y)*log(1 - mu);
    break;
  }
  case Link::probit:
  {
    boost::math::normal norm(0, 1);
    double lfk = glmmr::maths::log_factorial_approx(y);
    double lfn = glmmr::maths::log_factorial_approx(var_par);
    double lfnk = glmmr::maths::log_factorial_approx(var_par - y);
    logl = lfn - lfk - lfnk + y*((double)cdf(norm,mu)) + (var_par - y)*log(1 - (double)cdf(norm,mu));
    break;
  }
  default:
  {
    double lfk = glmmr::maths::log_factorial_approx(y);
    double lfn = glmmr::maths::log_factorial_approx(var_par);
    double lfnk = glmmr::maths::log_factorial_approx(var_par - y);
    logl = lfn - lfk - lfnk + y*log(1/(1+exp(-1.0*mu))) + (var_par - y)*log(1 - 1/(1+exp(-1.0*mu)));
    break;
  }
  }
    break;
  }
  case Fam::gaussian:
  {
    switch(family.link){
  case Link::loglink:
    logl = -0.5*log(var_par) -0.5*log(2*3.141593) -
      0.5*(log(y) - mu)*(log(y) - mu)/var_par;
    break;
  default:
    //identity
    logl = -0.5*log(var_par) -0.5*log(2*3.141593) -
      0.5*(y - mu)*(y - mu)/var_par;
    break;
  }
    break;
  }
  case Fam::gamma:
  {
    switch(family.link){
  case Link::inverse:
  {
    double ymu = var_par*y*mu;
    logl = log(1/(tgamma(var_par)*y)) + var_par*log(ymu) - ymu;
    break;
  }
  case Link::identity:
    logl = log(1/(tgamma(var_par)*y)) + var_par*log(var_par*y/mu) - var_par*y/mu;
    break;
  default:
    //log
    {
      double ymu = var_par*y/exp(mu);
      logl = log(1/(tgamma(var_par)*y)) + var_par*log(ymu) - ymu;
      break;
    }
  }
    break;
  }
  case Fam::beta:
  {
    //only logit currently
    logl = (mu*var_par - 1)*log(y) + ((1-mu)*var_par - 1)*log(1-y) - lgamma(mu*var_par) - lgamma((1-mu)*var_par) + lgamma(var_par);
    break;
  }
  case Fam::quantile: case Fam::quantile_scaled:
  {
    double resid = y - mod_inv_func(mu,family.link);
    // if(family.family == Fam::quantile_scaled) resid *= (1.0/var_par);
    logl = resid <= 0 ? resid*(1.0 - family.quantile) : -1.0*resid*family.quantile;
    break;
  }
  }
  return logl;
}

inline double logdet(const Eigen::MatrixXd& M) {
  double ld = 0;
  Eigen::LLT<Eigen::MatrixXd> chol(M);
  auto& U = chol.matrixL();
  for (int i = 0; i < M.rows(); ++i)
    ld += log(U(i,i));
  ld *= 2;
  return ld;
}

inline MatrixXd sample_MVN(const VectorMatrix& mu,
                           int m) {
  int n = mu.vec.size();
  MatrixXd L = mu.mat.llt().matrixL();
  boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
    generator(boost::mt19937(time(0)),
              boost::normal_distribution<>());
  VectorXd z(n);
  MatrixXd samps(n, m);
  for (int i = 0; i < m; i++) {
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
