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

//glmmrbase version 
#define GLMMR10
#define GLMMR11

// includes

#include <vector>
#include <array>
#include <string>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <set>
#include <map>
#include <unordered_map>
#include <random>
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

// inline bool isalnum_or_uscore(const char& s)
// {
//   return (isalnum(s) || s=='_');
// }
// 
// template<typename T>
// inline bool expect_number_of_unique_elements(const std::vector<T> vec,
//                                              int n){
//   int vec_size = std::set<T>(vec.begin(),vec.end()).size();
//   return vec_size==n;
// }

}

