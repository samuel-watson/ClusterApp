#pragma once

#include "general.h"
#include "covariance.hpp"
#include "hsgpcovariance.hpp"
#include "nngpcovariance.hpp"
#include "linearpredictor.hpp"
#include "family.hpp"
#include "modelextradata.hpp"
#include "formula.hpp"

namespace glmmr {

using namespace Eigen;

template<typename cov, typename linpred>
class ModelBits{
public:
  glmmr::Formula        formula;
  linpred               linear_predictor;
  cov                   covariance;
  glmmr::ModelExtraData data;
  glmmr::Family         family;
  bool                  weighted = false;
  int                   trace = 0;
    
  ModelBits(const std::string& formula_,const ArrayXXd& data_,const strvec& colnames_,std::string family_,std::string link_);
  //functions
  virtual int       n() const {return linear_predictor.n();};
  virtual ArrayXd   xb() {return linear_predictor.xb().array() + data.offset.array();};
  virtual void      make_covariance_sparse(bool amd = true);
  virtual void      make_covariance_dense();
};

}

template<>
inline glmmr::ModelBits<glmmr::Covariance, glmmr::LinearPredictor>::ModelBits(const std::string& formula_,
          const ArrayXXd& data_,
          const strvec& colnames_,
          std::string family_, 
          std::string link_) : 
  formula(formula_), 
  linear_predictor(formula,data_,colnames_),
  covariance(formula,data_,colnames_),
  data(data_.rows()),
  family(family_,link_) {
  covariance.linear_predictor_ptr(&linear_predictor);
};

template<>
inline glmmr::ModelBits<glmmr::nngpCovariance, glmmr::LinearPredictor>::ModelBits(const std::string& formula_,
                                                 const ArrayXXd& data_,
                                                 const strvec& colnames_,
                                                 std::string family_, 
                                                 std::string link_) : 
  formula(formula_), 
  linear_predictor(formula,data_,colnames_),
  covariance(formula,data_,colnames_),
  data(data_.rows()),
  family(family_,link_) {};

template<>
inline glmmr::ModelBits<glmmr::hsgpCovariance, glmmr::LinearPredictor>::ModelBits(const std::string& formula_,
                                                                                  const ArrayXXd& data_,
                                                                                  const strvec& colnames_,
                                                                                  std::string family_, 
                                                                                  std::string link_) : 
  formula(formula_), 
  linear_predictor(formula,data_,colnames_),
  covariance(formula,data_,colnames_),
  data(data_.rows()),
  family(family_,link_) {};

template<typename cov, typename linpred>
inline void glmmr::ModelBits<cov, linpred>::make_covariance_sparse(bool amd){
  covariance.set_sparse(true,amd);
}

template<typename cov, typename linpred>
inline void glmmr::ModelBits<cov, linpred>::make_covariance_dense(){
  covariance.set_sparse(false);
}

typedef glmmr::Covariance covariance;
typedef glmmr::nngpCovariance nngp;
typedef glmmr::hsgpCovariance hsgp;
typedef glmmr::LinearPredictor xb;
typedef glmmr::ModelBits<covariance, xb> bits;
typedef glmmr::ModelBits<nngp, xb> bits_nngp;
typedef glmmr::ModelBits<hsgp, xb> bits_hsgp;