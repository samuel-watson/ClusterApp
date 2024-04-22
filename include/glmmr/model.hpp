#pragma once

#include "general.h"
#include "modelbits.hpp"
#include "randomeffects.hpp"
#include "modelmatrix.hpp"
#include "modelmcmc.hpp"
#include "modeloptim.hpp"

namespace glmmr {

using namespace Eigen;

template<class>
struct check_type : std::false_type {};

template<>
struct check_type<glmmr::ModelBits<glmmr::Covariance, glmmr::LinearPredictor> > : std::true_type {};

template<typename modeltype>
class Model {
public:
  // model objects
  modeltype                       model;
  glmmr::RandomEffects<modeltype> re;
  glmmr::ModelMatrix<modeltype>   matrix;
  glmmr::ModelOptim<modeltype>    optim;
  glmmr::ModelMCMC<modeltype>     mcmc;
  // constructor
  Model(const std::string& formula_,const ArrayXXd& data_,const strvec& colnames_,std::string family_,std::string link_);
  //functions
  virtual void    set_offset(const VectorXd& offset_);
  virtual void    set_weights(const ArrayXd& weights_);
  virtual void    set_y(const VectorXd& y_);
  virtual void    update_beta(const dblvec &beta_);
  virtual void    update_theta(const dblvec &theta_);
  virtual void    update_u(const MatrixXd &u_, bool append = false);
  virtual void    reset_u(); // just resets the random effects samples to zero
  virtual void    set_trace(int trace_);
  virtual dblpair marginal(const MarginType type,
                             const std::string& x,
                             const strvec& at,
                             const strvec& atmeans,
                             const strvec& average,
                             const RandomEffectMargin re_type,
                             const SE se_type,
                             const dblpair& xvals,
                             const dblvec& atvals,
                             const dblvec& atrevals);
                                             
};

}

template<typename modeltype>
inline glmmr::Model<modeltype>::Model(const std::string& formula_,
      const ArrayXXd& data_,
      const strvec& colnames_,
      std::string family_, 
      std::string link_) : model(formula_,data_,colnames_,family_,link_), 
      re(model), 
      matrix(model,re,check_type<modeltype>::value,check_type<modeltype>::value),  
      optim(model,matrix,re), mcmc(model,matrix,re) {};

template<typename modeltype>
inline void glmmr::Model<modeltype>::set_offset(const VectorXd& offset_){
  model.data.set_offset(offset_);
}

template<typename modeltype>
inline void glmmr::Model<modeltype>::set_weights(const ArrayXd& weights_){
  model.data.set_weights(weights_);
  if((weights_ != 1.0).any()){
    model.weighted = true;
  }
}

template<typename modeltype>
inline void glmmr::Model<modeltype>::set_y(const VectorXd& y_){
  model.data.update_y(y_);
}

template<typename modeltype>
inline void glmmr::Model<modeltype>::update_beta(const dblvec &beta_){
  optim.update_beta(beta_);
}

template<typename modeltype>
inline void glmmr::Model<modeltype>::update_theta(const dblvec &theta_){
  model.covariance.update_parameters(theta_);
  re.zu_ = model.covariance.ZLu(re.u_);
  // model.vcalc.data = model.covariance.ZL();
}

template<typename modeltype>
inline void glmmr::Model<modeltype>::reset_u(){
  re.u_.resize(model.covariance.Q(),1);
  re.u_.setZero();
  re.zu_.resize(NoChange,1);
  re.zu_.setZero();
}

template<typename modeltype>
inline void glmmr::Model<modeltype>::update_u(const MatrixXd &u_, bool append){
#ifdef R_BUILD
  if(u_.rows()!=model.covariance.Q())Rcpp::stop(std::to_string(u_.rows())+" rows provided, "+std::to_string(model.covariance.Q())+" expected");
#endif
  
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
  // if(optim.ll_previous.rows() != optim.ll_current.rows()) optim.ll_previous.resize(optim.ll_current.rows(),NoChange);
  // optim.ll_previous = optim.ll_current;
  
  if(action_append){
    re.u_.conservativeResize(NoChange,currcolsize + newcolsize);
    re.zu_.conservativeResize(NoChange,currcolsize + newcolsize);
    re.u_.rightCols(newcolsize) = u_;
    optim.ll_current.resize(currcolsize + newcolsize,NoChange);
  } else {
    if(u_.cols()!=re.u_.cols()){
      #if defined(ENABLE_DEBUG) && defined(R_BUILD)
        Rcpp::Rcout << "\nResize u: " << model.covariance.Q() << "x" << u_.cols();
      #endif
      re.u_.resize(NoChange,newcolsize);
      re.zu_.resize(NoChange,newcolsize);
    }
    re.u_ = u_;
    if(re.u_.cols() != optim.ll_current.rows()) optim.ll_current.resize(newcolsize,NoChange);
  }
  re.zu_ = model.covariance.ZLu(re.u_);
}

template<typename modeltype>
inline void glmmr::Model<modeltype>::set_trace(int trace_){
  optim.trace = trace_;
  mcmc.trace = trace_;
  model.trace = trace_;
  if(trace_ > 0){
    mcmc.verbose = true;
  } else {
    mcmc.verbose = false;
  }
}

// marginal effects:
// type is margin type 
// x is name of variable to calculate margin
// at is fixed effects at a set value specified in atvals
// atmeans specifies the fixed effects to set at their mean value
// average specifies the fixed effects to average over
// re_type is random effects margin type
// se_type is the standard error type
// xvals values of the x variable at which to evaluate marginal effect
// atvals is the values for at argument
// atrevals is the random effects values if re_type is At
template<typename modeltype>
inline dblpair glmmr::Model<modeltype>::marginal(const MarginType type,
                                                 const std::string& x,
                                                 const strvec& at,
                                                 const strvec& atmeans,
                                                 const strvec& average,
                                                 const RandomEffectMargin re_type,
                                                 const SE se_type,
                                                 const dblpair& xvals,
                                                 const dblvec& atvals,
                                                 const dblvec& atrevals){
#pragma omp declare reduction(vec_dbl_plus : std::vector<double> :                                    \
  std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
  initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
  
  int total_p = at.size() + atmeans.size() + average.size() + 1;
  int intercept = 1- (int)model.linear_predictor.form.RM_INT;
  
#ifdef R_BUILD
  if(total_p != (model.linear_predictor.P() - intercept))Rcpp::stop("All variables must be named");
  if(at.size() != atvals.size())Rcpp::stop("Not enough values specified for at");
  if(re_type == RandomEffectMargin::Average && re.zu_.cols()<=1)Rcpp::warning("No MCMC samples of random effects. Random effects will be set at estimated values.");
#endif
    
  bool single_row = true;
  MatrixXd newXdata(1,model.linear_predictor.calc.data.cols());
  int P = model.linear_predictor.P();
  int N = 1;
  auto xidx = std::find(model.linear_predictor.calc.data_names.begin(),model.linear_predictor.calc.data_names.end(),x);
  int xcol = xidx - model.linear_predictor.calc.data_names.begin();
  
  if(average.size() > 0 || (total_p == 1 && (re_type == RandomEffectMargin::Average || re_type == RandomEffectMargin::AtEstimated))){
    single_row = false;
    N = model.n();
    newXdata.conservativeResize(model.n(),NoChange);
    newXdata.col(xcol) = model.linear_predictor.calc.data.col(xcol);
    
#ifdef R_BUILD
   if(re_type == RandomEffectMargin::At && atrevals.size() != model.covariance.Q())Rcpp::stop("Need to provide values for u vector");
#endif
   
   for(const auto& p: average){
     auto colidx = std::find(model.linear_predictor.calc.data_names.begin(),model.linear_predictor.calc.data_names.end(),p);
     if(colidx != model.linear_predictor.calc.data_names.end()){
       int pcol = colidx - model.linear_predictor.calc.data_names.begin();
       newXdata.col(pcol) = model.linear_predictor.calc.data.col(pcol);
     } else {
#ifdef R_BUILD
       Rcpp::stop("Variable "+p+" not in data names");  
#endif
     }
   }
  } else {
    newXdata(0,xcol) = xvals.first;
#ifdef R_BUILD
    if(re_type == RandomEffectMargin::At && atrevals.size() != 1)Rcpp::stop("Need to provide single value for Zu");
    if(re_type == RandomEffectMargin::AtEstimated)Rcpp::stop("All covariates are at fixed values, cannot used estimated random effects.");
#endif
  }
  
  if(at.size() > 0){
    for(int p = 0; p < at.size(); p++){
      auto colidx = std::find(model.linear_predictor.calc.data_names.begin(),model.linear_predictor.calc.data_names.end(),at[p]);
      if(colidx != model.linear_predictor.calc.data_names.end()){
        int pcol = colidx - model.linear_predictor.calc.data_names.begin();
        for(int i = 0; i < newXdata.rows(); i++){
          newXdata(i,pcol) = atvals[p];
        }
      } else {
#ifdef R_BUILD
        Rcpp::stop("Variable "+at[p]+" not in data names");  
#endif
      }
    }
  }
  if(atmeans.size() > 0){
    for(int p = 0; p < atmeans.size(); p++){
      auto colidx = std::find(model.linear_predictor.calc.data_names.begin(),model.linear_predictor.calc.data_names.end(),atmeans[p]);
      if(colidx != model.linear_predictor.calc.data_names.end()){
        int pcol = colidx - model.linear_predictor.calc.data_names.begin();
        double xmean = 0;
        for(int i = 0; i < model.n(); i++) xmean += model.linear_predictor.calc.data(i,pcol);
        xmean *= (1.0/model.n());
        for(int i = 0; i < newXdata.rows(); i++){
          newXdata(i,pcol) = xmean;
        }
      } else {
#ifdef R_BUILD
        Rcpp::stop("Variable "+atmeans[p]+" not in data names");  
#endif
      }
    }
  }
  
  // now create the new calculator object
  glmmr::calculator mcalc = model.linear_predictor.calc;
  mcalc.instructions.push_back(Do::PushExtraData);
  mcalc.instructions.push_back(Do::Add);
  glmmr::linear_predictor_to_link(mcalc,model.family.link);
  mcalc.data.conservativeResize(newXdata.rows(),NoChange);
  mcalc.data = newXdata;
  
  dblpair result;
  VectorXd delta = VectorXd::Zero(P);
  MatrixXd M(P,P);
  
#if defined(R_BUILD) && defined(ENABLE_DEBUG)
  Rcpp::Rcout << "\nMARGINS\nN: " << N << " xcol: " << xcol;
  Rcpp::Rcout << "\nTest calculator\nUsing X: " << newXdata.row(0);
  dblvec test_result = mcalc.calculate<CalcDyDx::XBeta>(0,0,xcol,0);
  Rcpp::Rcout << "\nValues: ";
  for(const auto& i: test_result) Rcpp::Rcout << i << " ";
#endif
  
  switch(se_type){
    case SE::KR:
      {
      CorrectionData<SE::KR> kdata = matrix.template small_sample_correction<SE::KR>();
      M = kdata.vcov_beta;
      break;
      }
    case SE::KR2:
    {
      CorrectionData<SE::KR2> kdata = matrix.template small_sample_correction<SE::KR2>();
      M = kdata.vcov_beta;
      break;
    }
    case SE::Robust:
      M = matrix.sandwich_matrix();
      break;
    default:
      M = matrix.information_matrix();
      M = M.llt().solve(MatrixXd::Identity(M.rows(),M.cols()));
      break;
  }
  
#if defined(R_BUILD) && defined(ENABLE_DEBUG)
  Rcpp::Rcout << "\nDone calculating SE, calculating margin...";
#endif
  
  switch(re_type){
  case RandomEffectMargin::At: case RandomEffectMargin::AtEstimated: case RandomEffectMargin::AtZero:
  {
    VectorXd zu(N);
    if(re_type == RandomEffectMargin::At){
      if(single_row){
        zu(0) = atrevals[0];
      } else {
        VectorXd u(model.covariance.Q());
        for(int i = 0; i < model.covariance.Q(); i++)u(i) = atrevals[i];
        zu = model.covariance.Z()*u;
      }
    } else if(re_type == RandomEffectMargin::AtEstimated) {
      zu = re.zu_.rowwise().mean();
    } else if(re_type == RandomEffectMargin::AtZero){
      zu.setZero();
    }
    
    switch(type){
    case MarginType::DyDx:
    {
      double d_result = 0;
      dblvec m_result(2+2*P);
      dblvec delta_vec(P,0.0);
      for(int i = 0; i < newXdata.rows(); i++)newXdata(i,xcol) = xvals.first;
#pragma omp parallel for reduction(+:d_result) reduction(vec_dbl_plus:delta_vec) private(m_result)
      for(int i = 0; i < N; i++){
        newXdata(i,xcol) = xvals.first;
        m_result = mcalc.calculate<CalcDyDx::XBeta>(i,0,xcol,zu(i));
        d_result += m_result[1];
        for(int p = 0; p < P; p++)delta_vec[p] += m_result[p+2+P];
      }
      result.first = d_result/N;
      for(int p = 0; p < P; p++)delta(p) = delta_vec[p]/N;
      result.second = sqrt((delta.transpose()*M*delta)(0));
      break;
    }
    case MarginType::Diff:
      {
      double d_result = 0;
      dblvec delta_vec(P,0.0);
      dblvec m_result(1+P);
      for(int i = 0; i < newXdata.rows(); i++)mcalc.data(i,xcol) = xvals.first;
#pragma omp parallel for reduction(+:d_result) reduction(vec_dbl_plus:delta_vec) private(m_result)
      for(int i = 0; i < N; i++)
      {
        m_result = mcalc.calculate<CalcDyDx::BetaFirst>(i,0,0,zu(i));
        d_result += m_result[0];
        for(int p = 0; p < P; p++)delta_vec[p] += m_result[p+1];
      }
      
      for(int j = 0; j < newXdata.rows(); j++)mcalc.data(j,xcol) = xvals.second;
      
#pragma omp parallel for reduction(+:d_result) reduction(vec_dbl_plus:delta_vec) private(m_result)      
      for(int i = 0; i < N; i++)
      {
        m_result = mcalc.calculate<CalcDyDx::BetaFirst>(i,0,0,zu(i));
        d_result += -1.0*m_result[0];
        for(int p = 0; p < P; p++)delta_vec[p] += -1.0*m_result[p+1];
      }
      result.first = d_result/N;
      for(int p = 0; p < P; p++)delta(p) = delta_vec[p]/N;
      result.second = sqrt((delta.transpose()*M*delta)(0));
      break;
      }
    case MarginType::Ratio:
    {
      double d_result0 = 0;
      double d_result1 = 0;
      dblvec delta0(P,0);
      dblvec delta1(P,0);
      dblvec m_result0(P+1);
      dblvec m_result1(P+1);
      for(int i = 0; i < newXdata.rows(); i++) newXdata(i,xcol) = xvals.first;
#pragma omp parallel for private(m_result0) reduction(+:d_result0) reduction(vec_dbl_plus:delta0)  
      for(int i = 0; i < N; i++)
      {
        m_result0 = mcalc.calculate<CalcDyDx::BetaFirst>(i,0,0,zu(i));
        d_result0 += m_result0[0];
        for(int p = 0; p < P; p++) delta0[p] += m_result0[p+1];
      }
      
      for(int i = 0; i < newXdata.rows(); i++) mcalc.data(i,xcol) = xvals.second;
#pragma omp parallel for private(m_result1) reduction(+:d_result1) reduction(vec_dbl_plus:delta1) 
      for(int i = 0; i < N; i++)
      {
        m_result1 = mcalc.calculate<CalcDyDx::BetaFirst>(i,0,0,zu(i));
        d_result1 += m_result1[0];
        for(int p = 0; p < P; p++) delta1[p] += m_result1[p+1];
      }
      
      result.first = log(d_result0/N) - log(d_result1/N);
      for(int p = 0; p < P; p++)delta(p) = delta0[p]/d_result0 - delta1[p]/d_result1;
      result.second = sqrt((delta.transpose()*M*delta)(0));
      break;
    }
    
    }
    break;
  }
    case RandomEffectMargin::Average:
      {
      int iter = re.zu_.cols();
      switch(type){
      case MarginType::DyDx:
      {
        double d_result = 0;
        dblvec m_result(2+2*P);
        dblvec delta_vec(P,0.0);
        for(int i = 0; i < newXdata.rows(); i++) newXdata(i,xcol) = xvals.first;
#pragma omp parallel for private(m_result) reduction(+:d_result) reduction(vec_dbl_plus:delta_vec) collapse(2)
        for(int i = 0; i < model.n(); i++){
          for(int j = 0; j < iter; j++){
            if(N==1){
              m_result = mcalc.calculate<CalcDyDx::XBeta>(0,0,xcol,re.zu_(i,j));
            } else {
              m_result = mcalc.calculate<CalcDyDx::XBeta>(i,0,xcol,re.zu_(i,j));
            }
            d_result += m_result[1];
            for(int p = 0; p < P; p++)delta_vec[p] += m_result[p+2+P];
          }
        }
        result.first = d_result/(N*iter);
        for(int p = 0; p < P; p++)delta(p) = delta_vec[p];
        delta.array() *= (1.0/(N*iter));
        result.second = sqrt((delta.transpose()*M*delta)(0));
        break;
      }
      case MarginType::Diff:
      {
        double d_result = 0;
        dblvec m_result(P+1);
        dblvec delta_vec(P,0.0);
        for(int i = 0; i < newXdata.rows(); i++)newXdata(i,xcol) = xvals.first;
#pragma omp parallel for private(m_result) reduction(+:d_result) reduction(vec_dbl_plus:delta_vec) collapse(2)
        for(int i = 0; i < model.n(); i++)
        {
          for(int j = 0; j < iter; j++)
          {
            if(N==1){
              m_result = mcalc.calculate<CalcDyDx::BetaFirst>(0,0,0,re.zu_(i,j));
            } else {
              m_result = mcalc.calculate<CalcDyDx::BetaFirst>(i,0,0,re.zu_(i,j));
            }
            d_result += m_result[0];
            for(int p = 0; p < P; p++)delta_vec[p] += m_result[p+1];
           }
        }
        
        for(int i = 0; i < newXdata.rows(); i++)mcalc.data(i,xcol) = xvals.second;
#pragma omp parallel for private(m_result) reduction(+:d_result) reduction(vec_dbl_plus:delta_vec) collapse(2)
        for(int i = 0; i < model.n(); i++)
        {
          for(int j = 0; j < iter; j++)
          {
            if(N==1){
              m_result = mcalc.calculate<CalcDyDx::BetaFirst>(0,0,0,re.zu_(i,j));
            } else {
              m_result = mcalc.calculate<CalcDyDx::BetaFirst>(i,0,0,re.zu_(i,j));
            }
            d_result += -1.0*m_result[0];
            for(int p = 0; p < P; p++)delta_vec[p] += -1.0*m_result[p+1];
          }
        }
        result.first = d_result/(N*iter);
        for(int p = 0; p < P; p++)delta(p) = delta_vec[p];
        delta.array() *= (1.0/(N*iter));
        result.second = sqrt((delta.transpose()*M*delta)(0));
        break;
      }
      case MarginType::Ratio:
      {
        double d_result0 = 0;
        double d_result1 = 0;
        dblvec delta0(P,0);
        dblvec delta1(P,0);
        dblvec m_result0(1+P);
        dblvec m_result1(1+P);
        for(int i = 0; i < newXdata.rows(); i++)newXdata(i,xcol) = xvals.first;
#pragma omp parallel for private(m_result0) reduction(+:d_result0) reduction(vec_dbl_plus:delta0) collapse(2)
        for(int i = 0; i < model.n(); i++)
        {
          for(int j = 0; j < iter; j++)
          {
            if(N==1){
              m_result0 = mcalc.calculate<CalcDyDx::BetaFirst>(1,0,0,re.zu_(i,j));
            } else {
              m_result0 = mcalc.calculate<CalcDyDx::BetaFirst>(i,0,0,re.zu_(i,j));
            }
            d_result0 += m_result0[0];
            for(int p = 0; p < P; p++) delta0[p] += m_result0[p+1];
          }
        }
        
        for(int i = 0; i < newXdata.rows(); i++)mcalc.data(i,xcol) = xvals.second;
#pragma omp parallel for private(m_result1) reduction(+:d_result1) reduction(vec_dbl_plus:delta1) collapse(2)
        for(int i = 0; i < model.n(); i++)
        {
          for(int j = 0; j < iter; j++)
          {
            if(N==1){
              m_result1 = mcalc.calculate<CalcDyDx::BetaFirst>(1,0,0,re.zu_(i,j));
            } else {
              m_result1 = mcalc.calculate<CalcDyDx::BetaFirst>(i,0,0,re.zu_(i,j));
            }
            d_result1 += m_result1[0];
            for(int p = 0; p < P; p++) delta1[p] += m_result1[p+1];
          }
        }
        
        result.first = log(d_result0/(N*iter)) - log(d_result1/(N*iter));
        for(int p = 0; p < P; p++){
          delta(p) = delta0[p]/d_result0 - delta1[p]/d_result1;
        }
        result.second = sqrt((delta.transpose()*M*delta)(0));
        break;
      }
      }
      break;
      }
  }
  
  return result;
  
}

