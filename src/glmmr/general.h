#ifndef GENERAL_H
#define GENERAL_H

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS 
#define _USE_MATH_DEFINES

// includes
#include <Eigen/Core>
#include <vector>
#include <string>
#include <cstring>
#include <sstream>
#include <iostream>
#include <regex>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <stack>
#include <set>
#include <unordered_map>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <random>
#include <boost/math/special_functions/digamma.hpp>
#include "rbobyqa.h"
#include "sparsechol.h"


using namespace Eigen;

typedef std::string str;
typedef std::vector<str> strvec;
typedef std::vector<int> intvec;
typedef std::vector<double> dblvec;
typedef std::vector<strvec> strvec2d;
typedef std::vector<dblvec> dblvec2d;
typedef std::vector<intvec> intvec2d;
typedef std::vector<dblvec2d> dblvec3d;
typedef std::vector<intvec2d> intvec3d;

namespace glmmr {

//useful things used in a few places
const static std::unordered_map<str, double> nvars = {  
  {"gr", 1},
  {"ar", 2},
  {"fexp0", 1},
  {"fexp", 2},
  {"sqexp0",1},
  {"sqexp",2},
  {"bessel",1},
  {"matern",2},
  {"wend0",2},
  {"wend1",2},
  {"wend2",2},
  {"prodwm",2},
  {"prodcb",2},
  {"prodek",2},
  {"ar0", 1},
  {"ar1", 1}
};

const static std::unordered_map<str,int> string_to_case{
  {"gr",1},
  {"ar",2},
  {"fexp0", 3},
  {"fexp", 4},
  {"sqexp0",5},
  {"sqexp",6},
  {"bessel",7},
  {"matern",8},
  {"wend0",9},
  {"wend1",10},
  {"wend2",11},
  {"prodwm",12},
  {"prodcb",13},
  {"prodek",14},
  {"ar0",15},
  {"ar1",16}
};

inline bool validate_fn(const str& fn){
  bool not_fn = string_to_case.find(fn) == string_to_case.end();
  return not_fn;
}

const static intvec xvar_rpn = {0,1,4,17};

inline bool is_number(const std::string& s)
{
  std::string::const_iterator it = s.begin();
  while (it != s.end() && std::isdigit(*it)) ++it;
  return !s.empty() && it == s.end();
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
}

struct vector_matrix{
public:
  VectorXd vec;
  MatrixXd mat;
  vector_matrix(int n): vec(n), mat(n,n) {};
  vector_matrix(const vector_matrix& x) : vec(x.vec), mat(x.mat) {};
  vector_matrix& operator=(vector_matrix x){
    vec = x.vec;
    mat = x.mat;
    return *this;
  };
};

struct matrix_matrix{
public:
  MatrixXd mat1;
  MatrixXd mat2;
  matrix_matrix(int n1, int m1, int n2, int m2): mat1(n1,m1), mat2(n2,m2) {};
  matrix_matrix(const matrix_matrix& x) : mat1(x.mat1), mat2(x.mat2) {};
  matrix_matrix& operator=(matrix_matrix x){
    mat1 = x.mat1;
    mat2 = x.mat2;
    return *this;
  };
};

struct kenward_data {
public:
    MatrixXd vcov_beta;
    MatrixXd vcov_theta;
    VectorXd dof;
    VectorXd lambda;
    kenward_data(int n1, int m1, int n2, int m2) : vcov_beta(n1, m1), vcov_theta(n2, m2), dof(n1), lambda(n1) {};
    kenward_data(const kenward_data& x) : vcov_beta(x.vcov_beta), vcov_theta(x.vcov_theta), dof(x.dof), lambda(x.lambda) {};
    kenward_data& operator=(kenward_data x) {
        vcov_beta = x.vcov_beta;
        vcov_theta = x.vcov_theta;
        dof = x.dof;
        lambda = x.lambda;
        return *this;
    };
};

#endif
