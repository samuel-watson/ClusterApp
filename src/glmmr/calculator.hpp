#pragma once

#include "general.h"
#include "openmpheader.h"
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
  std::vector<Do> instructions;
  intvec indexes;
  dblvec y;
  std::array<double,10> numbers = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  strvec parameter_names;
  strvec data_names;
  ArrayXd variance = ArrayXd::Constant(1,1.0);
  int data_count = 0;
  int parameter_count = 0;
  int user_number_count = 0;
  bool any_nonlinear = false;
  calculator() {};
  
  template<CalcDyDx dydx>
  dblvec calculate(const int i, 
                   const dblvec& parameters, 
                   const MatrixXd& data,
                   const int j = 0,
                   const int parameterIndex = 0,
                   const double extraData = 0.0,
                   const int n = 0) const;
  
  calculator& operator= (const glmmr::calculator& calc);
  VectorXd linear_predictor(const dblvec& parameters,const MatrixXd& data);
  MatrixXd jacobian(const dblvec& parameters, const MatrixXd& data);
  MatrixXd jacobian(const dblvec& parameters,const MatrixXd& data,const VectorXd& extraData);
  MatrixXd jacobian(const dblvec& parameters,const MatrixXd& data,const MatrixXd& extraData);
  MatrixMatrix jacobian_and_hessian(const dblvec& parameters,const MatrixXd& data,const MatrixXd& extraData);
  VectorMatrix jacobian_and_hessian(const dblvec& parameters);
  void print_instructions() const;
  void print_names(bool data = true, bool parameters = false) const;
};

}

inline VectorXd glmmr::calculator::linear_predictor(const dblvec& parameters, 
                                                    const MatrixXd& data){
  int n = data.rows();
  VectorXd x(n);
#pragma omp parallel for
  for(int i = 0; i < n; i++){
    x(i) = calculate<CalcDyDx::None>(i,parameters,data)[0];
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
  return *this;
};

inline MatrixXd glmmr::calculator::jacobian(const dblvec& parameters,
                                            const MatrixXd& data){
  int n = data.rows();
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
  if(n==0)Rcpp::stop("No data initialised in calculator");
#endif
  MatrixXd J(n,parameter_count);
#pragma omp parallel for
  for(int i = 0; i<n ; i++){
    dblvec out = calculate<CalcDyDx::BetaFirst>(i,parameters,data);
    for(int j = 0; j<parameter_count; j++){
      J(i,j) = out[j+1];
    }
  }
  return J;
};

inline MatrixXd glmmr::calculator::jacobian(const dblvec& parameters,
                                            const MatrixXd& data,
                                            const VectorXd& extraData){
  int n = data.rows();
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
  if(n==0)Rcpp::stop("No data initialised in calculator");
#endif 
  MatrixXd J(n,parameter_count);
#pragma omp parallel for
  for(int i = 0; i<n ; i++){
    dblvec out = calculate<CalcDyDx::BetaFirst>(i,parameters,data,0,0,extraData(i));
    for(int j = 0; j<parameter_count; j++){
      J(i,j) = out[j+1];
    }
  }
  return J;
};

inline MatrixXd glmmr::calculator::jacobian(const dblvec& parameters,
                                            const MatrixXd& data,
                                            const MatrixXd& extraData){
  int n = data.rows();
  
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
  if(n==0)Rcpp::stop("No data initialised in calculator");
  if(extraData.rows()!=n)Rcpp::stop("Extra data not of length n");
#endif
  
  int iter = extraData.cols();
  MatrixXd J = MatrixXd::Zero(parameter_count,n);
#pragma omp parallel for
  for(int i = 0; i<n ; i++){
    dblvec out;
    for(int k = 0; k < iter; k++){
      out = calculate<CalcDyDx::BetaFirst>(i,parameters,data,0,0,extraData(i,k));
      for(int j = 0; j < parameter_count; j++){
        J(j,i) += out[1+j]/iter;
      }
    }
  }
  return J;
};


inline MatrixMatrix glmmr::calculator::jacobian_and_hessian(const dblvec& parameters,
                                                             const MatrixXd& data,
                                                             const MatrixXd& extraData){
  int n = data.rows();
  MatrixMatrix result(parameter_count,parameter_count,parameter_count,n);
  
#ifdef ENABLE_DEBUG
  if(n==0)Rcpp::stop("No data initialised in calculator");
  if(extraData.rows()!=n)Rcpp::stop("Extra data not of length n");
#endif
  
  int iter = extraData.cols();
  int n2d = parameter_count*(parameter_count + 1)/2;
  MatrixXd H = MatrixXd::Zero(n2d,n);
  MatrixXd J = MatrixXd::Zero(parameter_count,n);
#pragma omp parallel for collapse(2)
  for(int i = 0; i<n ; i++){
    for(int k = 0; k < iter; k++){
      dblvec out = calculate<CalcDyDx::BetaSecond>(i,parameters,data,0,0,extraData(i,k));
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

inline VectorMatrix glmmr::calculator::jacobian_and_hessian(const dblvec& parameters){
  VectorMatrix result(parameter_count);
  int n2d = parameter_count*(parameter_count + 1)/2;
  VectorXd H = VectorXd::Zero(n2d);
  VectorXd J = VectorXd::Zero(parameter_count);
  MatrixXd dat = MatrixXd::Zero(1,1);
  dblvec out = calculate<CalcDyDx::BetaSecond>(0,parameters,dat,0,2,0);
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

template<CalcDyDx dydx>
inline dblvec glmmr::calculator::calculate(const int i, 
                                           const dblvec& parameters, 
                                           const MatrixXd& data,
                                           const int j,
                                           const int parameterIndex,
                                           const double extraData,
                                           const int n) const {
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
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(idx_iter >= indexes.size())Rcpp::stop("Index out of range: case 0 idx iter: "+std::to_string(idx_iter)+" versus "+std::to_string(indexes.size()));
      if(indexes[idx_iter] >= data.cols())Rcpp::stop("Index out of range: case 0 indexes: "+std::to_string(indexes[idx_iter])+" versus "+std::to_string(data.cols()));
      if(i >= data.rows())Rcpp::stop("Row index out of range: case 0: "+std::to_string(i)+" versus "+std::to_string(data.rows()));
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
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(idx_iter >= indexes.size())Rcpp::stop("Index out of range: case 1 idx iter: "+std::to_string(idx_iter)+" versus "+std::to_string(indexes.size()));
      if(indexes[idx_iter] >= data.cols())Rcpp::stop("Index out of range: case 1 indexes: "+std::to_string(indexes[idx_iter])+" versus "+std::to_string(data.cols()));
#endif
      
      if(i==j){
        stack.push(0.0);
      } else {
        int i1 = i < j ? (n-1)*i - ((i-1)*i/2) + (j-i-1) : (n-1)*j - ((j-1)*j/2) + (i-j-1);
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
        if(i1 >= data.rows())Rcpp::stop("PushCovData: Index out of range: "+std::to_string(i1)+" versus "+std::to_string(data.rows()));
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
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if((unsigned)idx_iter >= indexes.size())Rcpp::stop("Index out of range: case 2 idx iter: "+std::to_string(idx_iter)+" versus "+std::to_string(indexes.size()));
      if(indexes[idx_iter] >= parameter_count)Rcpp::stop("Index out of range: case 2 indexes: "+std::to_string(indexes[idx_iter])+" versus "+std::to_string(parameter_count));
      if(indexes[idx_iter] >= parameters.size())Rcpp::stop("Index out of range (pars): case 2 indexes: "+std::to_string(indexes[idx_iter])+" versus "+std::to_string(parameters.size()));
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
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(stack.size()<2)Rcpp::stop("Stack too small (3)");
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
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(stack.size()<2)Rcpp::stop("Stack too small (4)");
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
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(stack.size()<2)Rcpp::stop("Stack too small (5)");
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
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(stack.size()<2)Rcpp::stop("Stack too small (6)");
#endif
      
      a = stack.top();
      stack.pop();
      b = stack.top();
      
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(b == 0.0)Rcpp::stop("Divide by zero (6)");
      if(isinf(a/b))Rcpp::stop("Infinity in division (6b)");
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
      
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(stack.size()==0)Rcpp::stop("Stack too small (7)");
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
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(stack.size()<2)Rcpp::stop("Stack too small (8)");
#endif
      
      a = stack.top();
      stack.pop();
      b = stack.top();
      stack.pop();
      double out = pow(a,b);
      
#ifdef R_BUILD
      // try and catch these possible failures as it seems to cause crash if nan allowed to propogate
      if(out != out)Rcpp::stop("Exponent fail: "+std::to_string(a)+"^"+std::to_string(b));
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
          double result = pow(a,b)*b_top_dx.back()*log(a) + pow(a,b-1)*b*a_top_dx.back();
#ifdef R_BUILD
          // this can sometimes result in a crash if the values of the parameters aren't right
          if(result != result)Rcpp::stop("Exponent dydx fail: "+std::to_string(a)+"^"+std::to_string(b-1));
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
              second_dx[index_count].push(result1 + result3);
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
      
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(stack.size()==0)Rcpp::stop("Stack too small (9)");
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
      
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(stack.size()==0)Rcpp::stop("Stack too small (10)");
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
      
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(stack.size()==0)Rcpp::stop("Stack too small (11)");
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
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(stack.size()==0)Rcpp::stop("Stack too small (12)");
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
      
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(stack.size()==0)Rcpp::stop("Stack too small (13)");
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
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(stack.size()==0)Rcpp::stop("Stack too small (14)");
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
      
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(stack.size()==0)Rcpp::stop("Stack too small (15)");
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
    case Do::Log:
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(stack.size()==0)Rcpp::stop("Stack too small (16)");
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
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(stack.size()==0)Rcpp::stop("Stack too small (17)");
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
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(stack.size()==0)Rcpp::stop("Stack too small (40)");
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
    }
    
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
    if(stack.size() == 0)Rcpp::stop("Error stack empty");
    if(stack.top() != stack.top() || isnan(stack.top()))Rcpp::stop("Calculation evaluates to NaN");
#endif
  }
  
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
  if(stack.size()>1)Rcpp::warning("More than one element on the stack at end of calculation");
#endif
  
  dblvec result;
  result.push_back(stack.top());
  
  if constexpr (dydx != CalcDyDx::None){
    for(const auto& fstack: first_dx){
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(fstack.size()==0)Rcpp::stop("Error derivative stack empty");
#endif
      result.push_back(fstack.top());
    }
  }
  if constexpr (dydx == CalcDyDx::BetaSecond || dydx == CalcDyDx::XBeta){
    for(const auto& sstack: second_dx){
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
      if(sstack.size()==0)Rcpp::stop("Error second derivative stack empty");
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
      Rcpp::Rcout << ": " << parameter_names[indexes[idx_iter]] << "\n";
      idx_iter++;
      break;
    case Do::PushData:
      Rcpp::Rcout << "(column " << data_names[indexes[idx_iter]] << ")\n";
      idx_iter++;
      break;
    case Do::PushCovData:
      Rcpp::Rcout << "(column " << indexes[idx_iter] << ")\n";
      idx_iter++;
      break;
    default:
      Rcpp::Rcout << "\n";
    }
    counter++;
  }
#endif
}

inline void glmmr::calculator::print_names(bool data, bool parameters) const{
#ifdef R_BUILD
  if(data){
    Rcpp::Rcout << "\nData names: \n";
    glmmr::print_vec_1d<strvec>(data_names);
  }
  if(data){
    Rcpp::Rcout << "\nParameter names: \n";
    glmmr::print_vec_1d<strvec>(parameter_names);
  }
#endif
}