#pragma once

#include <stack>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/erf.hpp>
#include "maths.h"
#include "instructions.h"

enum class CalcDyDx {
  None,
  BetaFirst,
  BetaSecond,
  XBeta,
  Zu
};

namespace glmmr {

class calculator {
public:
  std::vector<Do>       instructions; // vector of insructions to execute
  intvec                indexes; // indexes of data or parameter vectors
  dblvec                y;  // outcome data
  std::array<double,20> numbers = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  strvec                parameter_names; //vector of parameter names
  strvec                data_names; // vector of data names
  ArrayXd               variance = ArrayXd::Constant(1,1.0); // variance values
  int                   data_count = 0; // number of data items
  int                   parameter_count = 0; // number of parameters
  int                   user_number_count = 0; // number of numbers in the function
  int                   data_size = 0; // number of data items in the calculation
  bool                  any_nonlinear = false; // for linear predictor - any non-linear functions?
  MatrixXd              data = MatrixXd::Zero(1,1); // the data for the calculation
  dblvec                parameters;
  intvec                parameter_indexes;
  calculator() {};
  
  template<CalcDyDx dydx>
  dblvec calculate(const int i, 
                   const int j = 0,
                   const int parameterIndex = 0,
                   const double extraData = 0.0) const;
  
  calculator& operator= (const glmmr::calculator& calc);
  VectorXd      linear_predictor();
  MatrixXd      jacobian();
  MatrixXd      jacobian(const VectorXd& extraData);
  MatrixXd      jacobian(const MatrixXd& extraData);
  MatrixMatrix  jacobian_and_hessian(const MatrixXd& extraData);
  VectorMatrix  jacobian_and_hessian();
  void          update_parameters(const dblvec& parameters_in);
  double        get_covariance_data(const int i, const int j, const int fn);
  void          print_instructions() const;
  void          print_names(bool print_data = true, bool print_parameters = false) const;
};

}

inline void glmmr::calculator::update_parameters(const dblvec& parameters_in){
  if(static_cast<int>(parameters_in.size()) < parameter_count)throw std::runtime_error("Expecting "+std::to_string(parameter_count)+" parameters in calculator but got "+std::to_string(parameters_in.size()));
  for(int i = 0; i < parameter_indexes.size(); i++)parameters[i] = parameters_in[parameter_indexes[i]];
}

inline VectorXd glmmr::calculator::linear_predictor(){
  int n = data.rows();
  VectorXd x(n);
#pragma omp parallel for
  for(int i = 0; i < n; i++){
    x(i) = calculate<CalcDyDx::None>(i)[0];
  }
  return x;
};

inline glmmr::calculator& glmmr::calculator::operator= (const glmmr::calculator& calc){
  instructions = calc.instructions;
  indexes = calc.indexes;
  parameter_names = calc.parameter_names;
  data_names = calc.data_names;
  variance.conservativeResize(calc.variance.size());
  variance = calc.variance;
  data_count = calc.data_count;
  parameter_count = calc.parameter_count;
  any_nonlinear = calc.any_nonlinear;
  data.conservativeResize(calc.data.rows(),calc.data.cols());
  data = calc.data;
  parameters.resize(calc.parameters.size());
  parameters = calc.parameters;
  return *this;
};

inline MatrixXd glmmr::calculator::jacobian(){
  int n = data.rows();
#ifdef ENABLE_DEBUG
  if(n==0)throw std::runtime_error("No data initialised in calculator");
#endif
  MatrixXd J(n,parameter_count);
#pragma omp parallel for
  for(int i = 0; i<n ; i++){
    dblvec out = calculate<CalcDyDx::BetaFirst>(i);
    for(int j = 0; j<parameter_count; j++){
      J(i,j) = out[j+1];
    }
  }
  return J;
};

inline MatrixXd glmmr::calculator::jacobian(const VectorXd& extraData){
  int n = data.rows();
#ifdef ENABLE_DEBUG
  if(n==0)throw std::runtime_error("No data initialised in calculator");
#endif 
  MatrixXd J(n,parameter_count);
#pragma omp parallel for
  for(int i = 0; i<n ; i++){
    dblvec out = calculate<CalcDyDx::BetaFirst>(i,0,0,extraData(i));
    for(int j = 0; j<parameter_count; j++){
      J(i,j) = out[j+1];
    }
  }
  return J;
};

inline MatrixXd glmmr::calculator::jacobian(const MatrixXd& extraData){
  int n = data.rows();
  
#ifdef ENABLE_DEBUG
  if(n==0)throw std::runtime_error("No data initialised in calculator");
  if(extraData.rows()!=n)throw std::runtime_error("Extra data not of length n");
#endif
  
  int iter = extraData.cols();
  MatrixXd J = MatrixXd::Zero(parameter_count,n);
#pragma omp parallel for
  for(int i = 0; i<n ; i++){
    dblvec out;
    for(int k = 0; k < iter; k++){
      out = calculate<CalcDyDx::BetaFirst>(i,0,0,extraData(i,k));
      for(int j = 0; j < parameter_count; j++){
        J(j,i) += out[1+j]/iter;
      }
    }
  }
  return J;
};


inline MatrixMatrix glmmr::calculator::jacobian_and_hessian(const MatrixXd& extraData){
  int n = data.rows();
  MatrixMatrix result(parameter_count,parameter_count,parameter_count,n);
  
#ifdef ENABLE_DEBUG
  if(n==0)throw std::runtime_error("No data initialised in calculator");
  if(extraData.rows()!=n)throw std::runtime_error("Extra data not of length n");
#endif
  
  int iter = extraData.cols();
  int n2d = parameter_count*(parameter_count + 1)/2;
  MatrixXd H = MatrixXd::Zero(n2d,n);
  MatrixXd J = MatrixXd::Zero(parameter_count,n);
#pragma omp parallel for collapse(2)
  for(int i = 0; i<n ; i++){
    for(int k = 0; k < iter; k++){
      dblvec out = calculate<CalcDyDx::BetaSecond>(i,0,0,extraData(i,k));
      for(int j = 0; j < parameter_count; j++){
        J(j,i) += out[1+j]/iter;
      }
      for(int j = 0; j < n2d; j++){
        H(j,i) += out[parameter_count + 1 + j]/iter;
      }
    }
  }
  VectorXd Hmean = H.rowwise().sum();
  MatrixXd H0 = MatrixXd::Zero(parameter_count, parameter_count);
  int index_count = 0;
  for(int j = 0; j < parameter_count; j++){
    for(int k = j; k < parameter_count; k++){
      H0(k,j) = Hmean[index_count];
      if(j != k) H0(j,k) = H0(k,j);
      index_count++;
    }
  }
  result.mat1 = H0;
  result.mat2 = J;
  return result;
};

inline VectorMatrix glmmr::calculator::jacobian_and_hessian(){
  VectorMatrix result(parameter_count);
  int n2d = parameter_count*(parameter_count + 1)/2;
  VectorXd H = VectorXd::Zero(n2d);
  VectorXd J = VectorXd::Zero(parameter_count);
  MatrixXd dat = MatrixXd::Zero(1,1);
  dblvec out = calculate<CalcDyDx::BetaSecond>(0,0,2,0);
  for(int j = 0; j < parameter_count; j++){
    J(j,0) += out[1+j];
  }
  for(int j = 0; j < n2d; j++){
    H(j) += out[parameter_count + 1 + j];
  }
  MatrixXd H0 = MatrixXd::Zero(parameter_count, parameter_count);
  int index_count = 0;
  for(int j = 0; j < parameter_count; j++){
    for(int k = j; k < parameter_count; k++){
      H0(k,j) = H(index_count);
      if(j != k) H0(j,k) = H0(k,j);
      index_count++;
    }
  }
  result.mat = H0;
  result.vec = J;
  return result;
};

inline double glmmr::calculator::get_covariance_data(const int i, const int j, const int fn){
  int i1 = i < j ? (data_size-1)*i - ((i-1)*i/2) + (j-i-1) : (data_size-1)*j - ((j-1)*j/2) + (i-j-1);
#ifdef ENABLE_DEBUG
  if(i1 >= data.rows())throw std::runtime_error("PushCovData: Index out of range: "+std::to_string(i1)+" versus "+std::to_string(data.rows()));
#endif
  return data(i1,fn);
}

template<CalcDyDx dydx>
inline dblvec glmmr::calculator::calculate(const int i, 
                                           const int j,
                                           const int parameterIndex,
                                           const double extraData) const {
  int idx_iter = 0;
  double a,b;
  std::stack<double> stack;
  std::vector<std::stack<double> > first_dx;
  std::vector<std::stack<double> > second_dx;
  
  if constexpr(dydx != CalcDyDx::None){
    if constexpr(dydx == CalcDyDx::BetaFirst || dydx == CalcDyDx::BetaSecond){
      first_dx.resize(parameter_count);
    } else if constexpr(dydx == CalcDyDx::XBeta){
      first_dx.resize(1+parameter_count);
    } else if constexpr(dydx == CalcDyDx::Zu){
      first_dx.resize(1);
    }
  }
  
  if constexpr(dydx == CalcDyDx::XBeta || dydx == CalcDyDx::BetaSecond){
    if constexpr(dydx == CalcDyDx::BetaSecond){
      second_dx.resize(parameter_count*(parameter_count + 1)/2);
    } else if constexpr(dydx ==  CalcDyDx::XBeta){
      second_dx.resize(parameter_count);
    }
  }
  
  auto addZeroDx = [&] (){
    for(auto& fstack: first_dx){
      fstack.push(0.0);
    }
  };
  
  auto addZeroDx2 = [&] (){
    for(auto& sstack: second_dx){
      sstack.push(0.0);
    }
  };
  
  auto allDyDxZero = [&](){
    if constexpr (dydx != CalcDyDx::None)addZeroDx();
    if constexpr (dydx == CalcDyDx::BetaSecond || dydx == CalcDyDx::XBeta)addZeroDx2();
  };
  
  for(const auto& k: instructions){
    switch(k){
    case Do::PushData:
    {
      // debugging statements to find possible errors
#ifdef ENABLE_DEBUG
      if(idx_iter >= indexes.size())throw std::runtime_error("Index out of range: case 0 idx iter: "+std::to_string(idx_iter)+" versus "+std::to_string(indexes.size()));
      if(indexes[idx_iter] >= data.cols())throw std::runtime_error("Index out of range: case 0 indexes: "+std::to_string(indexes[idx_iter])+" versus "+std::to_string(data.cols()));
      if(i >= data.rows())throw std::runtime_error("Row index out of range: case 0: "+std::to_string(i)+" versus "+std::to_string(data.rows()));
#endif
      stack.push(data(i,indexes[idx_iter]));
      if constexpr (dydx == CalcDyDx::BetaFirst || dydx == CalcDyDx::BetaSecond)addZeroDx();
      if constexpr (dydx == CalcDyDx::BetaSecond)addZeroDx2();
      if constexpr (dydx == CalcDyDx::XBeta){
        if(parameterIndex == indexes[idx_iter]){
          first_dx[0].push(1.0);
        } else {
          first_dx[0].push(0.0);
        }
        for(int i = 0; i < parameter_count; i++)first_dx[i+1].push(0.0);
        addZeroDx2();
      }
      if constexpr (dydx == CalcDyDx::Zu)first_dx[0].push(0.0);
      idx_iter++;
      break;
    }
    case Do::PushCovData:
    {
#ifdef ENABLE_DEBUG
      if(idx_iter >= indexes.size())throw std::runtime_error("Index out of range: case 1 idx iter: "+std::to_string(idx_iter)+" versus "+std::to_string(indexes.size()));
      if(indexes[idx_iter] >= data.cols())throw std::runtime_error("Index out of range: case 1 indexes: "+std::to_string(indexes[idx_iter])+" versus "+std::to_string(data.cols()));
#endif
      
      if(i==j){
        stack.push(0.0);
      } else {
        int i1 = i < j ? (data_size-1)*i - ((i-1)*i/2) + (j-i-1) : (data_size-1)*j - ((j-1)*j/2) + (i-j-1);
#ifdef ENABLE_DEBUG
        if(i1 >= data.rows())throw std::runtime_error("PushCovData: Index out of range: "+std::to_string(i1)+" versus "+std::to_string(data.rows()));
#endif
        stack.push(data(i1,indexes[idx_iter]));
      }
      if constexpr (dydx != CalcDyDx::None)addZeroDx();
      if constexpr (dydx == CalcDyDx::BetaSecond || dydx == CalcDyDx::XBeta)addZeroDx2();
      idx_iter++;
      break;
    }
    case Do::PushParameter:
    {
#ifdef ENABLE_DEBUG
      if(idx_iter >= indexes.size())throw std::runtime_error("Index out of range: case 2 idx iter: "+std::to_string(idx_iter)+" versus "+std::to_string(indexes.size()));
      if(indexes[idx_iter] >= parameter_count)throw std::runtime_error("Index out of range: case 2 indexes: "+std::to_string(indexes[idx_iter])+" versus "+std::to_string(parameter_count));
      if(indexes[idx_iter] >= parameters.size())throw std::runtime_error("Index out of range (pars): case 2 indexes: "+std::to_string(indexes[idx_iter])+" versus "+std::to_string(parameters.size()));
#endif
      
      stack.push(parameters[indexes[idx_iter]]);
      if constexpr (dydx == CalcDyDx::BetaFirst || dydx == CalcDyDx::BetaSecond){
        for(int idx = 0; idx < parameter_count; idx++){
          if(idx == indexes[idx_iter]){
            first_dx[idx].push(1.0);
          } else {
            first_dx[idx].push(0.0);
          }
        }
      }
      if constexpr (dydx == CalcDyDx::BetaSecond)addZeroDx2();
      if constexpr (dydx == CalcDyDx::XBeta){
        first_dx[0].push(0.0);
        for(int idx = 0; idx < parameter_count; idx++){
          if(idx == indexes[idx_iter]){
            first_dx[idx+1].push(1.0);
          } else {
            first_dx[idx+1].push(0.0);
          }
        }
        addZeroDx2();
      }
      if constexpr (dydx == CalcDyDx::Zu)addZeroDx();
      idx_iter++;
      break;
    }
    case Do::Add:
    {
#ifdef ENABLE_DEBUG
      if(stack.size()<2)throw std::runtime_error("Stack too small (3)");
#endif
      a = stack.top();
      stack.pop();
      b = stack.top();
      stack.pop();
      stack.push(a+b);
      if constexpr (dydx != CalcDyDx::None){
        for(auto& fstack: first_dx){
          a = fstack.top();
          fstack.pop();
          b = fstack.top();
          fstack.pop();
          fstack.push(a+b);
        }
      }
      if constexpr (dydx == CalcDyDx::BetaSecond || dydx == CalcDyDx::XBeta){
        for(auto& sstack: second_dx){
          a = sstack.top();
          sstack.pop();
          b = sstack.top();
          sstack.pop();
          sstack.push(a+b);
        }
      }
      break;
    }
    case Do::Subtract:
    {
#ifdef ENABLE_DEBUG
      if(stack.size()<2)throw std::runtime_error("Stack too small (4)");
#endif
      a = stack.top();
      stack.pop();
      b = stack.top();
      stack.pop();
      stack.push(a-b);
      if constexpr (dydx != CalcDyDx::None){
        for(auto& fstack: first_dx){
          a = fstack.top();
          fstack.pop();
          b = fstack.top();
          fstack.pop();
          fstack.push(a-b);
        }
      }
      if constexpr (dydx == CalcDyDx::BetaSecond || dydx == CalcDyDx::XBeta){
        for(auto& sstack: second_dx){
          a = sstack.top();
          sstack.pop();
          b = sstack.top();
          sstack.pop();
          sstack.push(a-b);
        }
      }
      break;
    }
    case Do::Multiply:
    {
#ifdef ENABLE_DEBUG
      if(stack.size()<2)throw std::runtime_error("Stack too small (5)");
#endif
      
      a = stack.top();
      stack.pop();
      b = stack.top();
      stack.pop();
      stack.push(a*b);
      if constexpr (dydx != CalcDyDx::None){
        dblvec a_top_dx;
        dblvec b_top_dx;
        for(auto& fstack: first_dx){
          a_top_dx.push_back(fstack.top());
          fstack.pop();
          b_top_dx.push_back(fstack.top());
          fstack.pop();
          fstack.push(a*b_top_dx.back() + b*a_top_dx.back());
        }
        if constexpr (dydx == CalcDyDx::BetaSecond){
          int index_count = 0;
          for(int idx = 0; idx < parameter_count; idx++){
            for(int jdx = idx; jdx < parameter_count; jdx++){
              double adx2 = second_dx[index_count].top();
              second_dx[index_count].pop();
              double bdx2 = second_dx[index_count].top();
              second_dx[index_count].pop();
              double result = a*bdx2 + b*adx2 + a_top_dx[idx]*b_top_dx[jdx] + a_top_dx[jdx]*b_top_dx[idx];
              second_dx[index_count].push(result);
              index_count++;
            }
          }
        }
        if constexpr (dydx == CalcDyDx::XBeta){
          for(int jdx = 0; jdx < parameter_count; jdx++){
            double adx2 = second_dx[jdx].top();
            second_dx[jdx].pop();
            double bdx2 = second_dx[jdx].top();
            second_dx[jdx].pop();
            double result = a*bdx2 + b*adx2 + a_top_dx[0]*b_top_dx[jdx+1] + a_top_dx[jdx+1]*b_top_dx[0];
            second_dx[jdx].push(result);
          }
        }
      }
      break;
    }
    case Do::Divide:
    {
#ifdef ENABLE_DEBUG
      if(stack.size()<2)throw std::runtime_error("Stack too small (6)");
#endif
      
      a = stack.top();
      stack.pop();
      b = stack.top();
      
#ifdef ENABLE_DEBUG
      if(b == 0.0)throw std::runtime_error("Divide by zero (6)");
      if(isinf(a/b))throw std::runtime_error("Infinity in division (6b)");
#endif
      
      stack.pop();
      stack.push(a/b);
      if constexpr (dydx != CalcDyDx::None){
        dblvec a_top_dx;
        dblvec b_top_dx;
        for(auto& fstack: first_dx){
          a_top_dx.push_back(fstack.top());
          fstack.pop();
          b_top_dx.push_back(fstack.top());
          fstack.pop();
          double result = (b*a_top_dx.back() - a*b_top_dx.back())/(b*b);
          fstack.push(result);
        }
        if constexpr (dydx == CalcDyDx::BetaSecond){
          int index_count = 0;
          for(int idx = 0; idx < parameter_count; idx++){
            for(int jdx = idx; jdx < parameter_count; jdx++){
              double adx2 = second_dx[index_count].top();
              second_dx[index_count].pop();
              double bdx2 = second_dx[index_count].top();
              second_dx[index_count].pop();
              double result = (adx2*b - a_top_dx[idx]*b_top_dx[jdx]- a_top_dx[jdx]*b_top_dx[idx])/(b*b) - (a*bdx2*b - 2*b_top_dx[idx]*b_top_dx[jdx])/(b*b*b);
              second_dx[index_count].push(result);
              index_count++;
            }
          }
        }
        if constexpr (dydx == CalcDyDx::XBeta){
          for(int jdx = 0; jdx < parameter_count; jdx++){
            double adx2 = second_dx[jdx].top();
            second_dx[jdx].pop();
            double bdx2 = second_dx[jdx].top();
            second_dx[jdx].pop();
            double result = (adx2*b - a_top_dx[0]*b_top_dx[jdx+1]- a_top_dx[jdx+1]*b_top_dx[0])/(b*b) - (a*bdx2*b - 2*b_top_dx[0]*b_top_dx[jdx])/(b*b*b);
            second_dx[jdx].push(result);
          }
        }
      }
      break;
    }
    case Do::Sqrt:
      
#ifdef ENABLE_DEBUG
      if(stack.size()==0)throw std::runtime_error("Stack too small (7)");
#endif
      
      a = stack.top();
      stack.pop();
      stack.push(sqrt(a));
      if constexpr (dydx != CalcDyDx::None){
        dblvec a_top_dx;
        for(auto& fstack: first_dx){
          a_top_dx.push_back(fstack.top());
          fstack.pop();
          double result = a==0 ? 0 : 0.5*pow(a,-0.5)*a_top_dx.back();
          fstack.push(result);
        }
        if constexpr (dydx == CalcDyDx::BetaSecond){
          int index_count = 0;
          for(int idx = 0; idx < parameter_count; idx++){
            for(int jdx = idx; jdx < parameter_count; jdx++){
              double adx2 = second_dx[index_count].top();
              second_dx[index_count].pop();
              double result = a==0? 0 : 0.5*pow(a,-0.5)*adx2 - 0.25*a_top_dx[idx]*a_top_dx[jdx]*pow(a,-3/2);
              second_dx[index_count].push(result);
              index_count++;
            }
          }
        }
        if constexpr (dydx == CalcDyDx::XBeta){
          for(int jdx = 0; jdx < parameter_count; jdx++){
            double adx2 = second_dx[jdx].top();
            second_dx[jdx].pop();
            double result = a==0? 0 : 0.5*pow(a,-0.5)*adx2 - 0.25*a_top_dx[0]*a_top_dx[jdx+1]*pow(a,-3/2);
            second_dx[jdx].push(result);
          }
        }
      }
      break;
    case Do::Power:
    {
#ifdef ENABLE_DEBUG
      if(stack.size()<2)throw std::runtime_error("Stack too small (8)");
#endif
      
      a = stack.top();
      stack.pop();
      b = stack.top();
      stack.pop();
      double out = pow(a,b);
      
#ifdef R_BUILD
      // try and catch these possible failures as it seems to cause crash if nan allowed to propogate
      if(out != out)throw std::runtime_error("Exponent fail: "+std::to_string(a)+"^"+std::to_string(b));
#endif
      
      stack.push(out);
      
      if constexpr (dydx != CalcDyDx::None){
        dblvec a_top_dx;
        dblvec b_top_dx;
        for(auto& fstack: first_dx){
          a_top_dx.push_back(fstack.top());
          fstack.pop();
          b_top_dx.push_back(fstack.top());
          fstack.pop();
          double result = 0;
          if(a > 0) result += pow(a,b)*b_top_dx.back()*log(a) + pow(a,b-1)*b*a_top_dx.back(); 
#ifdef R_BUILD
          // this can sometimes result in a crash if the values of the parameters aren't right
          if(result != result)throw std::runtime_error("Exponent dydx fail: "+std::to_string(a)+"^"+std::to_string(b-1));
#endif
          fstack.push(result);
        }
        if constexpr (dydx == CalcDyDx::BetaSecond){
          int index_count = 0;
          for(int idx = 0; idx < parameter_count; idx++){
            for(int jdx = idx; jdx < parameter_count; jdx++){
              double adx2 = second_dx[index_count].top();
              second_dx[index_count].pop();
              double bdx2 = second_dx[index_count].top();
              second_dx[index_count].pop();
              double result1 = first_dx[jdx].top()*b_top_dx[idx]*log(a) + stack.top()*(bdx2*log(a) + b_top_dx[idx]*(1/a)*a_top_dx[jdx]);
              double result2 = pow(a,b-1)*b_top_dx[jdx]*log(a) + pow(a,b-2)*(b-1)*a_top_dx[jdx];
              double result3 = result2*b*a_top_dx[idx] + pow(a,b-1)*(b*adx2+b_top_dx[jdx]*a_top_dx[idx]);
              double result4 = 0;
              if(a > 0) result4 = result1 + result3;
              #ifdef R_BUILD
              // this can sometimes result in a crash if the values of the parameters aren't right
              if(result4 != result4)throw std::runtime_error("Exponent d2ydx2 fail: "+std::to_string(a)+"^"+std::to_string(b-2)+" (i,j) = ("+std::to_string(idx)+","+std::to_string(jdx)+")");
              #endif
              second_dx[index_count].push(result4);
              index_count++;
            }
          }
        }
        if constexpr (dydx == CalcDyDx::XBeta){
          for(int jdx = 0; jdx < parameter_count; jdx++){
            double adx2 = second_dx[jdx].top();
            second_dx[jdx].pop();
            double bdx2 = second_dx[jdx].top();
            second_dx[jdx].pop();
            double result1 = first_dx[jdx+1].top()*b_top_dx[0]*log(a) + stack.top()*(bdx2*log(a) + b_top_dx[0]*(1/a)*a_top_dx[jdx+1]);
            double result2 = pow(a,b-1)*b_top_dx[jdx+1]*log(a) + pow(a,b-2)*(b-1)*a_top_dx[jdx+1];
            double result3 = result2*b*a_top_dx[0] + pow(a,b-1)*(b*adx2+b_top_dx[jdx+1]*a_top_dx[0]);
            second_dx[jdx].push(result1 + result3);
          }
        }
      }
      break;
    }
    case Do::Exp:
      
#ifdef ENABLE_DEBUG
      if(stack.size()==0)throw std::runtime_error("Stack too small (9)");
#endif
      
      a = stack.top();
      stack.pop();
      stack.push(exp(a));
      if constexpr (dydx != CalcDyDx::None){
        dblvec a_top_dx;
        for(auto& fstack: first_dx){
          a_top_dx.push_back(fstack.top());
          fstack.pop();
          double result = stack.top()*a_top_dx.back();
          fstack.push(result);
        }
        if constexpr (dydx == CalcDyDx::BetaSecond){
          int index_count = 0;
          for(int idx = 0; idx < parameter_count; idx++){
            for(int jdx = idx; jdx < parameter_count; jdx++){
              double adx2 = second_dx[index_count].top();
              second_dx[index_count].pop();
              double result = stack.top()*(a_top_dx[idx]*a_top_dx[jdx] + adx2);
              second_dx[index_count].push(result);
              index_count++;
            }
          }
        }
        if constexpr (dydx == CalcDyDx::XBeta){
          for(int jdx = 0; jdx < parameter_count; jdx++){
            double adx2 = second_dx[jdx].top();
            second_dx[jdx].pop();
            double result = stack.top()*(a_top_dx[0]*a_top_dx[jdx+1] + adx2);
            second_dx[jdx].push(result);
          }
        }
      }
      break;
    case Do::Negate:
      
#ifdef ENABLE_DEBUG
      if(stack.size()==0)throw std::runtime_error("Stack too small (10)");
#endif
      
      a = stack.top();
      stack.pop();
      stack.push(-1*a);
      if constexpr (dydx != CalcDyDx::None){
        for(auto& fstack: first_dx){
          double ftop = fstack.top();
          fstack.pop();
          fstack.push(-1.0*ftop);
        }
      }
      if constexpr (dydx == CalcDyDx::BetaSecond || dydx == CalcDyDx::XBeta){
        for(auto& sstack: second_dx){
          double adx2 = sstack.top();
          sstack.pop();
          sstack.push(-1*adx2);
        }
      }
      break;
    case Do::Bessel:
      
#ifdef ENABLE_DEBUG
      if(stack.size()==0)throw std::runtime_error("Stack too small (11)");
#endif
      
      a = stack.top();
      stack.pop();
      b = boost::math::cyl_bessel_k(1,a);
      stack.push(b);
      if constexpr (dydx != CalcDyDx::None){
        dblvec a_top_dx;
        for(auto& fstack: first_dx){
          a_top_dx.push_back(fstack.top());
          fstack.pop();
          double result = -0.5*boost::math::cyl_bessel_k(0,a)-0.5*boost::math::cyl_bessel_k(2,a);
          fstack.push(result*a_top_dx.back());
        }
        if constexpr (dydx == CalcDyDx::BetaSecond){
          int index_count = 0;
          for(int idx = 0; idx < parameter_count; idx++){
            for(int jdx = idx; jdx < parameter_count; jdx++){
              double adx2 = second_dx[index_count].top();
              second_dx[index_count].pop();
              double result = 0.25*boost::math::cyl_bessel_k(-1,a)+0.5*boost::math::cyl_bessel_k(1,a)+0.25*boost::math::cyl_bessel_k(3,a);
              double result1 = -0.5*boost::math::cyl_bessel_k(0,a)-0.5*boost::math::cyl_bessel_k(2,a);
              double result2 = result*a_top_dx[idx]*a_top_dx[jdx]+result1*adx2;
              second_dx[index_count].push(result2);
              index_count++;
            }
          }
        }
        if constexpr (dydx == CalcDyDx::XBeta){
          for(int jdx = 0; jdx < parameter_count; jdx++){
            double adx2 = second_dx[jdx].top();
            second_dx[jdx].pop();
            double result = 0.25*boost::math::cyl_bessel_k(-1,a)+0.5*boost::math::cyl_bessel_k(1,a)+0.25*boost::math::cyl_bessel_k(3,a);
            double result1 = -0.5*boost::math::cyl_bessel_k(0,a)-0.5*boost::math::cyl_bessel_k(2,a);
            double result2 = result*a_top_dx[0]*a_top_dx[jdx+1]+result1*adx2;
            second_dx[jdx].push(result2);
          }
        }
      }
      break;
    case Do::Gamma:
#ifdef ENABLE_DEBUG
      if(stack.size()==0)throw std::runtime_error("Stack too small (12)");
#endif
      a = stack.top();
      stack.pop();
      stack.push(tgamma(a));
      if constexpr (dydx != CalcDyDx::None){
        dblvec a_top_dx;
        for(auto& fstack: first_dx){
          a_top_dx.push_back(fstack.top());
          fstack.pop();
          double result = stack.top()*boost::math::polygamma(0,a);
          fstack.push(result*a_top_dx.back());
        }
        if constexpr (dydx == CalcDyDx::BetaSecond){
          int index_count = 0;
          for(int idx = 0; idx < parameter_count; idx++){
            for(int jdx = idx; jdx < parameter_count; jdx++){
              double adx2 = second_dx[index_count].top();
              second_dx[index_count].pop();
              double result = stack.top()*boost::math::polygamma(0,a)*boost::math::polygamma(0,a) + stack.top()*boost::math::polygamma(1,a);
              double result1 = stack.top()*boost::math::polygamma(0,a);
              double result2 = result*a_top_dx[idx]*a_top_dx[jdx]+result1*adx2;
              second_dx[index_count].push(result2);
              index_count++;
            }
          }
        }
        if constexpr (dydx == CalcDyDx::XBeta){
          for(int jdx = 0; jdx < parameter_count; jdx++){
            double adx2 = second_dx[jdx].top();
            second_dx[jdx].pop();
            double result = stack.top()*boost::math::polygamma(0,a)*boost::math::polygamma(0,a) + stack.top()*boost::math::polygamma(1,a);
            double result1 = stack.top()*boost::math::polygamma(0,a);
            double result2 = result*a_top_dx[0]*a_top_dx[jdx+1]+result1*adx2;
            second_dx[jdx].push(result2);
          }
        }
      }
      break;
    case Do::Sin:
      
#ifdef ENABLE_DEBUG
      if(stack.size()==0)throw std::runtime_error("Stack too small (13)");
#endif
      
      a = stack.top();
      stack.pop();
      stack.push(sin(a));
      if constexpr (dydx != CalcDyDx::None){
        dblvec a_top_dx;
        for(auto& fstack: first_dx){
          a_top_dx.push_back(fstack.top());
          fstack.pop();
          double result = cos(a)*a_top_dx.back();
          fstack.push(result);
        }
        if constexpr (dydx == CalcDyDx::BetaSecond){
          int index_count = 0;
          for(int idx = 0; idx < parameter_count; idx++){
            for(int jdx = idx; jdx < parameter_count; jdx++){
              double adx2 = second_dx[index_count].top();
              second_dx[index_count].pop();
              double result = -1.0*sin(a);
              double result1 = cos(a);
              double result2 = result*a_top_dx[idx]*a_top_dx[jdx]+result1*adx2;
              second_dx[index_count].push(result2);
              index_count++;
            }
          }
        }
        if constexpr (dydx == CalcDyDx::XBeta){
          for(int jdx = 0; jdx < parameter_count; jdx++){
            double adx2 = second_dx[jdx].top();
            second_dx[jdx].pop();
            double result = -1.0*sin(a);
            double result1 = cos(a);
            double result2 = result*a_top_dx[0]*a_top_dx[jdx+1]+result1*adx2;
            second_dx[jdx].push(result2);
          }
        }
      }
      break;
    case Do::Cos:
#ifdef ENABLE_DEBUG
      if(stack.size()==0)throw std::runtime_error("Stack too small (14)");
#endif
      
      a = stack.top();
      stack.pop();
      stack.push(cos(a));
      if constexpr (dydx != CalcDyDx::None){
        dblvec a_top_dx;
        for(auto& fstack: first_dx){
          a_top_dx.push_back(fstack.top());
          fstack.pop();
          double result = -1.0*sin(a)*a_top_dx.back();
          fstack.push(result);
        }
        if constexpr (dydx == CalcDyDx::BetaSecond){
          int index_count = 0;
          for(int idx = 0; idx < parameter_count; idx++){
            for(int jdx = idx; jdx < parameter_count; jdx++){
              double adx2 = second_dx[index_count].top();
              second_dx[index_count].pop();
              double result = -1.0*cos(a);
              double result1 = -1.0*sin(a);
              double result2 = result*a_top_dx[idx]*a_top_dx[jdx]+result1*adx2;
              second_dx[index_count].push(result2);
              index_count++;
            }
          }
        }
        if constexpr (dydx == CalcDyDx::XBeta){
          for(int jdx = 0; jdx < parameter_count; jdx++){
            double adx2 = second_dx[jdx].top();
            second_dx[jdx].pop();
            double result = -1.0*cos(a);
            double result1 = -1.0*sin(a);
            double result2 = result*a_top_dx[0]*a_top_dx[jdx+1]+result1*adx2;
            second_dx[jdx].push(result2);
          }
        }
      }
      break;
    case Do::BesselK:
      
#ifdef ENABLE_DEBUG
      if(stack.size()==0)throw std::runtime_error("Stack too small (15)");
#endif
      
      a = stack.top();
      stack.pop();
      b = stack.top();
      stack.pop();
      stack.push(boost::math::cyl_bessel_k(b,a));
      if constexpr (dydx != CalcDyDx::None){
        dblvec a_top_dx;
        for(auto& fstack: first_dx){
          a_top_dx.push_back(fstack.top());
          fstack.pop();
          double result = -0.5*boost::math::cyl_bessel_k(b-1,a)-0.5*boost::math::cyl_bessel_k(b+1,a);
          fstack.push(result*a_top_dx.back());
        }
        if constexpr (dydx == CalcDyDx::BetaSecond){
          int index_count = 0;
          for(int idx = 0; idx < parameter_count; idx++){
            for(int jdx = idx; jdx < parameter_count; jdx++){
              double adx2 = second_dx[index_count].top();
              second_dx[index_count].pop();
              double result = 0.25*boost::math::cyl_bessel_k(b-2,a)+0.5*boost::math::cyl_bessel_k(b,a)+0.25*boost::math::cyl_bessel_k(b+2,a);
              double result1 = -0.5*boost::math::cyl_bessel_k(b-1,a)-0.5*boost::math::cyl_bessel_k(b+1,a);
              double result2 = result*a_top_dx[idx]*a_top_dx[jdx]+result1*adx2;
              second_dx[index_count].push(result2);
              index_count++;
            }
          }
        }
        if constexpr (dydx == CalcDyDx::XBeta){
          for(int jdx = 0; jdx < parameter_count; jdx++){
            double adx2 = second_dx[jdx].top();
            second_dx[jdx].pop();
            double result = 0.25*boost::math::cyl_bessel_k(b-2,a)+0.5*boost::math::cyl_bessel_k(b,a)+0.25*boost::math::cyl_bessel_k(b+2,a);
            double result1 = -0.5*boost::math::cyl_bessel_k(b-1,a)-0.5*boost::math::cyl_bessel_k(b+1,a);
            double result2 = result*a_top_dx[0]*a_top_dx[jdx+1]+result1*adx2;
            second_dx[jdx].push(result2);
          }
        }
      }
      break;
    case Do::ErrorFunc:
        a = stack.top();
        stack.pop();
        stack.push(boost::math::erf(a));
        if constexpr (dydx != CalcDyDx::None)
        {
          dblvec a_top_dx;
          for(auto& fstack: first_dx)
          {
            a_top_dx.push_back(fstack.top());
            fstack.pop();
            double result = 1.12837916709551257390 * exp(-1.0 * a * a);
            result *= a_top_dx.back();
            fstack.push(result);
          }
          if constexpr (dydx == CalcDyDx::BetaSecond)
          {
            int index_count = 0;
            for(int idx = 0; idx < parameter_count; idx++){
              for(int jdx = idx; jdx < parameter_count; jdx++){
                double adx2 = second_dx[index_count].top();
                second_dx[index_count].pop();
                double result = -2.0 * a * 1.12837916709551257390 * exp(-1.0 * a * a);
                double result1 = 1.12837916709551257390 * exp(-1.0 * a * a);
                double result2 = result*a_top_dx[idx]*a_top_dx[jdx]+result1*adx2;
                second_dx[index_count].push(result2);
                index_count++;
              }
            }
          }
          if constexpr (dydx == CalcDyDx::XBeta){
          for(int jdx = 0; jdx < parameter_count; jdx++){
            double adx2 = second_dx[jdx].top();
            second_dx[jdx].pop();
            double result = -2.0 * a * 1.12837916709551257390 * exp(-1.0 * a * a);
            double result1 = 1.12837916709551257390 * exp(-1.0 * a * a);
            double result2 = result*a_top_dx[0]*a_top_dx[jdx+1]+result1*adx2;
            second_dx[jdx].push(result2);
          }
        }
        }        
        break;
    case Do::Log:
#ifdef ENABLE_DEBUG
      if(stack.size()==0)throw std::runtime_error("Stack too small (16)");
#endif
      
      a = stack.top();
      stack.pop();
      stack.push(log(a));
      if constexpr (dydx != CalcDyDx::None){
        dblvec a_top_dx;
        for(auto& fstack: first_dx){
          a_top_dx.push_back(fstack.top());
          fstack.pop();
          double result = (1/a)*a_top_dx.back();
          fstack.push(result);
        }
        if constexpr (dydx == CalcDyDx::BetaSecond){
          int index_count = 0;
          for(int idx = 0; idx < parameter_count; idx++){
            for(int jdx = idx; jdx < parameter_count; jdx++){
              double adx2 = second_dx[index_count].top();
              second_dx[index_count].pop();
              double result = -1.0/(a*a);
              double result1 = 1/a;
              double result2 = result*a_top_dx[idx]*a_top_dx[jdx]+result1*adx2;
              second_dx[index_count].push(result2);
              index_count++;
            }
          }
        }
        if constexpr (dydx == CalcDyDx::XBeta){
          for(int jdx = 0; jdx < parameter_count; jdx++){
            double adx2 = second_dx[jdx].top();
            second_dx[jdx].pop();
            double result = -1.0/(a*a);
            double result1 = 1/a;
            double result2 = result*a_top_dx[0]*a_top_dx[jdx+1]+result1*adx2;
            second_dx[jdx].push(result2);
          }
        }
      }
      break;
    case Do::Square:
#ifdef ENABLE_DEBUG
      if(stack.size()==0)throw std::runtime_error("Stack too small (17)");
#endif
      a = stack.top();
      stack.pop();
      stack.push(a*a);
      if constexpr (dydx != CalcDyDx::None){
        dblvec a_top_dx;
        for(auto& fstack: first_dx){
          a_top_dx.push_back(fstack.top());
          fstack.pop();
          double result = 2*a*a_top_dx.back();
          fstack.push(result);
        }
        if constexpr (dydx == CalcDyDx::BetaSecond){
          int index_count = 0;
          for(int idx = 0; idx < parameter_count; idx++){
            for(int jdx = idx; jdx < parameter_count; jdx++){
              double adx2 = second_dx[index_count].top();
              second_dx[index_count].pop();
              double result2 = 2*a_top_dx[idx]*a_top_dx[jdx]+2*a*adx2;
              second_dx[index_count].push(result2);
              index_count++;
            }
          }
        }
        if constexpr (dydx == CalcDyDx::XBeta){
          for(int jdx = 0; jdx < parameter_count; jdx++){
            double adx2 = second_dx[jdx].top();
            second_dx[jdx].pop();
            double result2 = 2*a_top_dx[0]*a_top_dx[jdx+1]+2*a*adx2;
            second_dx[jdx].push(result2);
          }
        }
      }
      break;
    case Do::PushExtraData:
    {
      stack.push(extraData);
      allDyDxZero();
      break;
    }
    case Do::Sign:
    {
      a = data(i,indexes[idx_iter]);
      if(a > 0){
        stack.push(1.0);
      } else if(a< 0){
        stack.push(-1.0);
      } else {
        stack.push(0.0);
      }     
      allDyDxZero();
      idx_iter++;
      break;
    }
    case Do::SignNoZero:
    {
      a = data(i,indexes[idx_iter]);
      if(a >= 0){
        stack.push(1.0);
      } else {
        stack.push(-1.0);
      }    
      allDyDxZero();
      idx_iter++;
      break;
    }  
    case Do::PushY:
    {
      stack.push(y[i]);
      allDyDxZero();
      break;
    }
    case Do::Int10:
      stack.push(10);
      allDyDxZero();
      break;
    case Do::Int1:
      stack.push(1);
      allDyDxZero();
      break;
    case Do::Int2:
      stack.push(2);
      allDyDxZero();
      break;
    case Do::Int3:
      stack.push(3);
      allDyDxZero();
      break;
    case Do::Int4:
      stack.push(4);
      allDyDxZero();
      break;
    case Do::Int5:
      stack.push(5);
      allDyDxZero();
      break;
    case Do::Int6:
      stack.push(6);
      allDyDxZero();
      break;
    case Do::Int7:
      stack.push(7);
      allDyDxZero();
      break;
    case Do::Int8:
      stack.push(8);
      allDyDxZero();
      break;
    case Do::Int9:
      stack.push(9);
      allDyDxZero();
      break;
    case Do::Pi:
      stack.push(M_PI);
      allDyDxZero();
      break;
    case Do::Constant1:
      stack.push(0.3275911);
      allDyDxZero();
      break;
    case Do::Constant2:
      stack.push(0.254829592);
      allDyDxZero();
      break;
    case Do::Constant3:
      stack.push(-0.284496736);
      allDyDxZero();
      break;
    case Do::Constant4:
      stack.push(1.421413741);
      allDyDxZero();
      break;
    case Do::Constant5:
      stack.push(-1.453152027);
      allDyDxZero();
      break;
    case Do::Constant6:
      stack.push(1.061405429);
      allDyDxZero();
      break;
    case Do::LogFactorialApprox:
    {
      //log factorial approximation
#ifdef ENABLE_DEBUG
      if(stack.size()==0)throw std::runtime_error("Stack too small (40)");
#endif
      
      a = stack.top();
      stack.pop();
      // Ramanujan approximation
      if(a == 0){
        stack.push(0);
      } else {
        double result = a*log(a) - a + log(a*(1+4*a*(1+2*a)))/6 + log(3.141593)/2;
        stack.push(result);
      }
      // NOTE: this function is only ever used in Poisson/binom log likelihood and so the top of the derivative
      // stacks should be 0, so we don't need to do anything. However, this should be updated if ever this
      // function is used more broadly.
      //if(order > 0){
      //  addZeroDx();
      //}
      //if(order == 2){
      //  addZeroDx2();
      //}
      break;
    }
    case Do::PushVariance:
    {
      stack.push(variance(i));
      allDyDxZero();
      break;
    }
    case Do::PushUserNumber0:
      stack.push(numbers[0]);
      allDyDxZero();
      break;
    case Do::PushUserNumber1:
      stack.push(numbers[1]);
      allDyDxZero();
      break;
    case Do::PushUserNumber2:
      stack.push(numbers[2]);
      allDyDxZero();
      break;
    case Do::PushUserNumber3:
      stack.push(numbers[3]);
      allDyDxZero();
      break;
    case Do::PushUserNumber4:
      stack.push(numbers[4]);
      allDyDxZero();
      break;
    case Do::PushUserNumber5:
      stack.push(numbers[5]);
      allDyDxZero();
      break;
    case Do::PushUserNumber6:
      stack.push(numbers[6]);
      allDyDxZero();
      break;
    case Do::PushUserNumber7:
      stack.push(numbers[7]);
      allDyDxZero();
      break;
    case Do::PushUserNumber8:
      stack.push(numbers[8]);
      allDyDxZero();
      break;
    case Do::PushUserNumber9:
      stack.push(numbers[9]);
      allDyDxZero();
      break;
    case Do::PushUserNumber10:
      stack.push(numbers[10]);
      allDyDxZero();
      break;
    case Do::PushUserNumber11:
      stack.push(numbers[11]);
      allDyDxZero();
      break;
    case Do::PushUserNumber12:
      stack.push(numbers[12]);
      allDyDxZero();
      break;
    case Do::PushUserNumber13:
      stack.push(numbers[13]);
      allDyDxZero();
      break;
    case Do::PushUserNumber14:
      stack.push(numbers[14]);
      allDyDxZero();
      break;
    case Do::PushUserNumber15:
      stack.push(numbers[15]);
      allDyDxZero();
      break;
    case Do::PushUserNumber16:
      stack.push(numbers[16]);
      allDyDxZero();
      break;
    case Do::PushUserNumber17:
      stack.push(numbers[17]);
      allDyDxZero();
      break;
    case Do::PushUserNumber18:
      stack.push(numbers[18]);
      allDyDxZero();
      break;
    case Do::PushUserNumber19:
      stack.push(numbers[19]);
      allDyDxZero();
      break;
    case Do::SqrtTwo:
      stack.push(0.7071068);
      allDyDxZero();
      break;
    case Do::Half:
      stack.push(0.5);
      allDyDxZero();
      break;
    case Do::Pi2:
      stack.push(2*M_PI);
      allDyDxZero();
      break;  
    case Do::HalfLog2Pi:
      stack.push(0.9189385);
      allDyDxZero();
      break;    
    }
    
#ifdef ENABLE_DEBUG
    if(stack.size() == 0)throw std::runtime_error("Error stack empty");
    if(stack.top() != stack.top() || isnan(stack.top()))throw std::runtime_error("Calculation evaluates to NaN");
#endif
  }
  
#ifdef ENABLE_DEBUG
  if(stack.size()>1)Rcpp::warning("More than one element on the stack at end of calculation");
#endif
  
  dblvec result;
  result.push_back(stack.top());
  
  if constexpr (dydx != CalcDyDx::None){
    for(const auto& fstack: first_dx){
#ifdef ENABLE_DEBUG
      if(fstack.size()==0)throw std::runtime_error("Error derivative stack empty");
#endif
      result.push_back(fstack.top());
    }
  }
  if constexpr (dydx == CalcDyDx::BetaSecond || dydx == CalcDyDx::XBeta){
    for(const auto& sstack: second_dx){
#ifdef ENABLE_DEBUG
      if(sstack.size()==0)throw std::runtime_error("Error second derivative stack empty");
#endif
      result.push_back(sstack.top());
    }
  }
  
  return result;
}

inline void glmmr::calculator::print_instructions() const {
  //currently only setup for R
#ifdef R_BUILD
  int counter = 1;
  int idx_iter = 0;
  Rcpp::Rcout << "\nInstructions:\n";
  for(const auto& i: instructions){
    Rcpp::Rcout << counter << ". " << instruction_str.at(i);
    switch(i){
    case Do::PushUserNumber0:
      Rcpp::Rcout << " = " << numbers[0] << "\n";
      break;
    case Do::PushUserNumber1:
      Rcpp::Rcout << " = " << numbers[1] << "\n";
      break;
    case Do::PushUserNumber2:
      Rcpp::Rcout << " = " << numbers[2] << "\n";
      break;
    case Do::PushUserNumber3:
      Rcpp::Rcout << " = " << numbers[3] << "\n";
      break;
    case Do::PushUserNumber4:
      Rcpp::Rcout << " = " << numbers[4] << "\n";
      break;
    case Do::PushUserNumber5:
      Rcpp::Rcout << " = " << numbers[5] << "\n";
      break;
    case Do::PushUserNumber6:
      Rcpp::Rcout << " = " << numbers[6] << "\n";
      break;
    case Do::PushUserNumber7:
      Rcpp::Rcout << " = " << numbers[7] << "\n";
      break;
    case Do::PushUserNumber8:
      Rcpp::Rcout << " = " << numbers[8] << "\n";
      break;
    case Do::PushUserNumber9:
      Rcpp::Rcout << " = " << numbers[9] << "\n";
      break;
    case Do::PushParameter:
      {
        if(indexes[idx_iter] >= parameter_names.size()){
          Rcpp::Rcout << "\nError in instruction set";
          Rcpp::Rcout << "\nIndex " << indexes[idx_iter] << " requested for parameter size " << parameter_names.size();
          Rcpp::Rcout << "\nIndexes: ";
          glmmr::print_vec_1d(indexes);
          Rcpp::Rcout << "\nParameter names: ";
          glmmr::print_vec_1d(parameter_names);
          throw std::runtime_error("Execution halted");
        }
        Rcpp::Rcout << ": " << parameter_names[indexes[idx_iter]] << "; index " << indexes[idx_iter] <<"\n";
        idx_iter++;
        break;
      }
    case Do::PushData: case Do::Sign: case Do::SignNoZero:
      {
        if(indexes[idx_iter] >= data_names.size()){
          Rcpp::Rcout << "\nError in instruction set";
          Rcpp::Rcout << "\nIndex " << indexes[idx_iter] << " requested for data size " << data_names.size();
          Rcpp::Rcout << "\nIndexes: ";
          glmmr::print_vec_1d(indexes);
          Rcpp::Rcout << "\nData names: ";
          glmmr::print_vec_1d(data_names);
          throw std::runtime_error("Execution halted");
        }
        Rcpp::Rcout << " (column " << data_names[indexes[idx_iter]] << "; index " << indexes[idx_iter] <<")\n";
        idx_iter++;
        break;
      }
      
    case Do::PushCovData:
      Rcpp::Rcout << " (column " << indexes[idx_iter] << ")\n";
      idx_iter++;
      break;
    default:
      Rcpp::Rcout << "\n";
    }
    counter++;
  }
  Rcpp::Rcout << "\n";
#endif
}

inline void glmmr::calculator::print_names(bool print_data, bool print_parameters) const {
#ifdef R_BUILD
  Rcpp::Rcout << "\nParameter count " << parameter_count << " vec size: " << parameters.size();
  Rcpp::Rcout << "\nData count " << data_count << " mat size: " << data.rows() << " x " << data.cols();
  Rcpp::Rcout << "\nIndexes: ";
  glmmr::print_vec_1d(indexes);
  Rcpp::Rcout << "\nAny nonlinear? " << any_nonlinear;
  if(print_data)
  {
    Rcpp::Rcout << "\nData names: ";
    glmmr::print_vec_1d<strvec>(data_names);
  }
  if(print_parameters)
  {
    Rcpp::Rcout << "\nParameter names: ";
    glmmr::print_vec_1d<strvec>(parameter_names);
  }
  VectorXd x(10);
  for(int i = 0; i < 10; i++) x(i) = calculate<CalcDyDx::None>(i)[0];
  Rcpp::Rcout << "\nExample data: " << x.transpose() << "\n";
  
#endif
}