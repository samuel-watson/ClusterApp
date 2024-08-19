#pragma once

//defines

#define _USE_MATH_DEFINES
// #define ENABLE_DEBUG // COMMENT/UNCOMMENT FOR DEBUG - currently only useful in R builds as uses R print, will add more general error logging
// #define R_BUILD //Uncomment to build for R with RCPP

#ifdef R_BUILD
#include <RcppEigen.h>
#else
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#endif

#ifdef __clang__
#define EIGEN_HAS_STD_RESULT_OF 0 // This has no effect with RcppEigen as it has Eigen <0.3.4
#endif
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS 
#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int
// includes

#include <vector>
#include <array>
#include <string>
#include <cstring>
#include <sstream>
#include <regex>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <stack>
#include <variant>
#include <set>
#include <unordered_map>
#include <queue>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <random>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/random.hpp>
#include "sparsematrix.h"
#include "operators.h"
#include "sparsechol.h"

using namespace Eigen;
using namespace SparseOperators;

typedef std::string str;
typedef std::vector<str> strvec;
typedef std::vector<int> intvec;
typedef std::vector<double> dblvec;
typedef std::vector<strvec> strvec2d;
typedef std::vector<dblvec> dblvec2d;
typedef std::vector<intvec> intvec2d;
typedef std::vector<dblvec2d> dblvec3d;
typedef std::vector<intvec2d> intvec3d;
typedef std::pair<double, double> dblpair;
typedef std::pair<std::string, double> strdblpair;

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

namespace glmmr {

enum class CovFunc {
    gr = 0,
    ar = 1,
    fexp0 = 2,
    fexp = 3,
    sqexp0 = 4,
    sqexp = 5,
    bessel = 6,
    matern = 7,
    truncpow2 = 8,
    truncpow3 = 9,
    truncpow4 = 10,
    cauchy = 11,
    cauchy3 = 12,
    truncpow20 = 13,
    truncpow30 = 14,
    truncpow40 = 15,
    cauchy0 = 16,
    cauchy30 = 17,
    ar0 = 18,
    ar1 = 19,
    dist = 20
};

enum class Fam {
    gaussian = 0,
    bernoulli = 1,
    poisson = 2,
    gamma = 3,
    beta = 4,
    binomial = 5,
    quantile = 6, // quantile is the asymmetric Laplacian distribution
    quantile_scaled = 7
};

enum class Link {
   logit = 0,
    loglink = 1, // to avoid conflicting with log() function
    probit = 2,
    identity = 3,
    inverse = 4
};

const std::map<str, Fam> str_to_family = {
  {"gaussian",Fam::gaussian},
  {"bernoulli",Fam::bernoulli},
  {"poisson",Fam::poisson},
  {"gamma",Fam::gamma},
  {"Gamma",Fam::gamma},
  {"beta",Fam::beta},
  {"binomial",Fam::binomial},
  {"quantile",Fam::quantile},
  {"quantile_scaled",Fam::quantile_scaled}
};

const std::map<str, Link> str_to_link = {
  {"logit",Link::logit},
  {"log",Link::loglink},
  {"probit",Link::probit},
  {"identity",Link::identity},
  {"inverse",Link::inverse}
};

const std::map<str, CovFunc> str_to_covfunc = {
  {"gr", CovFunc::gr},
  {"ar", CovFunc::ar},
  {"fexp0", CovFunc::fexp0},
  {"fexp", CovFunc::fexp},
  {"sqexp0",CovFunc::sqexp0},
  {"sqexp",CovFunc::sqexp},
  {"bessel",CovFunc::bessel},
  {"matern",CovFunc::matern},
  {"truncpow2",CovFunc::truncpow2},
  {"truncpow3",CovFunc::truncpow3},
  {"truncpow4",CovFunc::truncpow4},
  {"cauchy",CovFunc::cauchy},
  {"cauchy3",CovFunc::cauchy3},
  {"truncpow20",CovFunc::truncpow20},
  {"truncpow30",CovFunc::truncpow30},
  {"truncpow40",CovFunc::truncpow40},
  {"cauchy0",CovFunc::cauchy0},
  {"cauchy30",CovFunc::cauchy30},
  {"ar0", CovFunc::ar0},
  {"ar1", CovFunc::ar1},
  {"dist",CovFunc::dist}
};

// unfortunately need bidirectional map so need to duplicate this unless there's
// a better way??
const std::map<CovFunc, str> covfunc_to_str = {
  {CovFunc::gr, "gr"},
  {CovFunc::ar, "ar"},
  {CovFunc::fexp0, "fexp0"},
  {CovFunc::fexp, "fexp"},
  {CovFunc::sqexp0, "sqexp0"},
  {CovFunc::sqexp, "sqexp"},
  {CovFunc::bessel, "bessel"},
  {CovFunc::matern, "matern"},
  {CovFunc::truncpow2, "truncpow2"},
  {CovFunc::truncpow3, "truncpow3"},
  {CovFunc::truncpow4, "truncpow4"},
  {CovFunc::cauchy, "cauchy"},
  {CovFunc::cauchy3, "cauchy3"},
  {CovFunc::truncpow20, "truncpow20"},
  {CovFunc::truncpow30, "truncpow30"},
  {CovFunc::truncpow40, "truncpow40"},
  {CovFunc::cauchy0, "cauchy0"},
  {CovFunc::cauchy30, "cauchy30"},
  {CovFunc::ar0, "ar0"},
  {CovFunc::ar1, "ar1"},
  {CovFunc::dist, "dist"}
};

const std::map<CovFunc, int> covfunc_to_nvar = {
  {CovFunc::gr, 1},
  {CovFunc::ar, 2},
  {CovFunc::fexp0, 1},
  {CovFunc::fexp, 2},
  {CovFunc::sqexp0, 1},
  {CovFunc::sqexp, 2},
  {CovFunc::bessel, 1},
  {CovFunc::matern, 2},
  {CovFunc::truncpow2, 2},
  {CovFunc::truncpow3, 2},
  {CovFunc::truncpow4, 2},
  {CovFunc::cauchy, 3},
  {CovFunc::cauchy3, 2},
  {CovFunc::truncpow20, 1},
  {CovFunc::truncpow30, 1},
  {CovFunc::truncpow40, 1},
  {CovFunc::cauchy0, 2},
  {CovFunc::cauchy30, 1},
  {CovFunc::ar0, 1},
  {CovFunc::ar1, 1},
  {CovFunc::dist, 0}
};

inline bool validate_fn(const str& fn){
  bool not_fn = str_to_covfunc.find(fn) == str_to_covfunc.end();
  return not_fn;
}

//const static intvec xvar_rpn = {0,1,4,17};
#ifdef R_BUILD

template<typename T>
inline void print_vec_1d(const T& vec){
  Rcpp::Rcout << "\n[1]: ";
  for(auto j: vec) Rcpp::Rcout << j << " ";
}

template<typename T>
inline void print_vec_2d(const T& vec){
  for(int i = 0; i < vec.size(); i++){
    Rcpp::Rcout << "\n[" << i << "]: ";
    for(auto j: vec[i]) Rcpp::Rcout << j << " ";
  }
}

template<typename T>
inline void print_vec_3d(const T& vec){
  for(int i = 0; i < vec.size(); i++){
    Rcpp::Rcout << "\n[" << i << "]:";
    for(int k = 0; k < vec[i].size(); k++){
      Rcpp::Rcout << "\n   [" << i <<","<< k << "]: ";
      for(auto j: vec[i][k]) Rcpp::Rcout << j << " ";
    }
  }
}

inline void print_sparse(const sparse& A){
  Rcpp::Rcout << "\nmatL Ap: ";
  for(auto i: A.Ap)Rcpp::Rcout << " " << i;
  Rcpp::Rcout << "\nmatL Ai: ";
  for(auto i: A.Ai)Rcpp::Rcout << " " << i;
  Rcpp::Rcout << "\nmatL Ax: ";
  for(auto i: A.Ax)Rcpp::Rcout << " " << i;
}

#endif

inline bool is_number(const std::string& s)
{
  bool isnum = true;
  try {
    float a = std::stod(s);
  }
  catch (std::invalid_argument const& ex)
  {
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
    Rcpp::Rcout << " Not double: " << ex.what() << '\n';
#endif
    isnum = false;
  }
  return isnum;
}

inline bool isalnum_or_uscore(const char& s)
{
  return (isalnum(s) || s=='_');
}

template<typename T>
inline bool expect_number_of_unique_elements(const std::vector<T> vec,
                                             int n){
  int vec_size = std::set<T>(vec.begin(),vec.end()).size();
  return vec_size==n;
}

inline bool is_compact_fn(const CovFunc& fn){
  bool compact = false;
  if(fn == CovFunc::truncpow2 || fn == CovFunc::truncpow3 || fn == CovFunc::truncpow4 || fn == CovFunc::truncpow20 || fn == CovFunc::truncpow30
       || fn == CovFunc::truncpow40 || fn == CovFunc::cauchy || fn == CovFunc::cauchy3 || fn == CovFunc::cauchy0  || fn == CovFunc::cauchy30) compact = true;
  return compact;
}

enum class MarginType {
  DyDx = 0,
    Diff = 1,
    Ratio = 2
};

enum class SE {
  GLS = 0,
  KR = 1,
  Robust = 2,
  BW = 3,
  KR2 = 4,
  Sat = 5,
  KRBoth = 6 // used for when two types of correction are required
};

}

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

struct CorrectionDataBase {
public:
  MatrixXd vcov_beta;
  MatrixXd vcov_theta;
  VectorXd dof;
  VectorXd lambda;
  CorrectionDataBase(int n1, int m1, int n2, int m2): vcov_beta(n1,m1), vcov_theta(n2,m2), dof(n1), lambda(n1) {};
  CorrectionDataBase(const MatrixXd& vcov_beta_, const MatrixXd& vcov_theta_, const MatrixXd& dof_, const MatrixXd& lambda_) : 
    vcov_beta(vcov_beta_), vcov_theta(vcov_theta_), dof(dof_), lambda(lambda_)  {};
  CorrectionDataBase(const CorrectionDataBase& x) : vcov_beta(x.vcov_beta), vcov_theta(x.vcov_theta), dof(x.dof), lambda(x.lambda) {};
  CorrectionDataBase& operator=(const CorrectionDataBase& x) = default;
};

template<glmmr::SE corr>
struct CorrectionData : public CorrectionDataBase {
public:
  CorrectionData(int n1, int m1, int n2, int m2): CorrectionDataBase(n1,m1,n2,m2) {};
  CorrectionData(const CorrectionData& x) : CorrectionDataBase(x.vcov_beta, x.vcov_theta, x.dof, x.lambda) {};
  CorrectionData& operator=(const CorrectionData& x){
    CorrectionDataBase::operator=(x);
    return *this;
  };
};

template<>
struct CorrectionData<glmmr::SE::KRBoth> : public CorrectionDataBase {
public:
  MatrixXd vcov_beta_second;
  CorrectionData(int n1, int m1, int n2, int m2): CorrectionDataBase(n1,m1,n2,m2), vcov_beta_second(n1,m1) {};
  CorrectionData(const CorrectionData& x) : CorrectionDataBase(x.vcov_beta, x.vcov_theta, x.dof, x.lambda), vcov_beta_second(x.vcov_beta_second) {};
  CorrectionData& operator=(const CorrectionData& x){
    CorrectionDataBase::operator=(x);
    vcov_beta_second = x.vcov_beta_second;
    return *this;
  };
};

struct BoxResults {
  dblvec dof;
  dblvec scale;
  dblvec test_stat;
  dblvec p_value;
  BoxResults(const int r) : dof(r), scale(r), test_stat(r), p_value(r) {};
};
